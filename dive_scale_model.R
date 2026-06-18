# Code for Johansson et al. 
# dive_scale_model

# 0. Load Libraries and Data
# 1. Multiple Imputation Loop
## 1a. set parameters
## 1b. aggregate and filter movement data
## 1c. calculate spatiotemporally-paired chla metric
## 1d. effort function
## 1e. efficiency function
## 1f. fit model
# 2. Pool Models

#############################
# 0. Load Libraries and Data
#############################

library(tidyverse)
library(sf)
library(data.table)
library(glmmTMB)
library(broom.mixed)
library(ggeffects)

load("byTrip_df.RData")
load("sim_dive_locs.RData")
load("chla.RData")
load("internal_states_df.RData")
load("byDive_PP_df.RData")

#############################
# 1. Multiple Imputation Loop
#############################
MI_list <- list()

# fit model 20 times using 20 different estimated dive locations
for(m in 1:20){ 
  
  #### 1a. set parameters ####
  thinning_threshold <- 4 # thinning dives to account for autocorrelation
  curr_succ_cutoff <- 4 # number of prior dives used to define current success
  effort_cutoff <- 4 # number of subsequent minutes used to define foraging effort
  PP_col <- "cluster16_PP95_TF" # metric of prey pursuit used
  imputation_id <- m 
  
  #### 1b. aggregate and filter movement data #### 
  penguin_data <- byDive_PP_df |>
    rename(dive_energy = energy_J) |> 
    left_join(byTrip_df,
              by = c('InstrumentSeq','TripNumber')) |>
    rename(trip_energy = energy_J) |>
    left_join(sim_dive_locs |> filter(sim_id == imputation_id) |> mutate(TripNumber=as.numeric(TripNumber)),
              by = c("InstrumentSeq",'TripNumber',"DiveNumber")) |>
    left_join(internal_states_df,
              by = c("InstrumentSeq")) |>
    filter(Duration >= 24,
           Duration < 144,
           !is.na(TripNumber),
           divetim > 10,
           maxdep > 5,
           MCP_Area_km2 < 25,
           !is.na(sim_id),
           divetim < 1000)
  
  #### 1c. calculate spatiotemporally-paired chla metric ####
  chla_res  <- 0.041666668
  
  # integer-index binning to avoid floating point drift
  snap_to_grid <- function(x, res) {
    round(round(x / res) * res, digits = 6)
  }
  
  # aggregate to gridded bins
  chla_reg <- chla |>
    mutate(
      date    = as.Date(time),
      lon_bin = snap_to_grid(longitude, chla_res),
      lat_bin = snap_to_grid(latitude,  chla_res),
      chlor_a_log = log(chlor_a + 0.0001)
    ) %>%
    group_by(date, lon_bin, lat_bin) %>%
    summarize(chlor_a = mean(chlor_a_log, na.rm = TRUE), .groups = "drop")

  # project  locations to WGS84
  penguin_data_utm <- st_as_sf(penguin_data, coords = c("mu.x", "mu.y"), crs = 32720)
  penguin_data_ll  <- st_transform(penguin_data_utm, crs = 4326)
  coords <- st_coordinates(penguin_data_ll)
  
  penguin_env_paired <- penguin_data_ll %>%
    mutate(
      Longitude    = coords[, 1],
      Latitude     = coords[, 2],
      date         = as.Date(GMTDateTime),
      lon_bin_chla = snap_to_grid(Longitude, chla_res),
      lat_bin_chla = snap_to_grid(Latitude,  chla_res)
    )
  
  # helper function
  add_broad_metrics <- function(penguin_df, env_reg, res, lon_col, lat_col, env_var, prefix, n = 3) {
    offsets <- seq(-(n %/% 2), n %/% 2) * res
    neighbors <- tidyr::crossing(
      lon_offset = offsets,
      lat_offset = offsets
    )
    
    broad_df <- penguin_df %>%
      mutate(row_id = row_number()) %>%
      crossing(neighbors) %>%
      mutate(
        lon_bin_n = snap_to_grid(Longitude, res) + round(lon_offset, 6),
        lat_bin_n = snap_to_grid(Latitude,  res) + round(lat_offset, 6)
      ) %>%
      left_join(env_reg, by = c("date", "lon_bin_n" = "lon_bin", "lat_bin_n" = "lat_bin")) %>%
      group_by(row_id) %>%
      summarize(
        "{prefix}_broad" := mean(.data[[env_var]], na.rm = TRUE),
        "{prefix}_sd"    := sd(.data[[env_var]],   na.rm = TRUE),
        "{prefix}_count" := sum(!is.na(.data[[env_var]])),
        .groups = "drop"
      )
    
    penguin_df %>%
      mutate(row_id = row_number()) %>%
      left_join(broad_df, by = "row_id") %>%
      select(-row_id)
  }
  
  # 3x3 for chla (~14 km)
  penguin_env_paired <- penguin_env_paired %>%
    add_broad_metrics(chla_reg, chla_res, "lon_bin_chla", "lat_bin_chla","chlor_a","chla",n = 3) 
  
  #### 1d. effort function ####
  calc_effort <- function(dives, time_cutoff) {
    time_cutoff <- time_cutoff*60 # convert min to sec
    
    dives <- dives %>%
      arrange(InstrumentSeq, begdesc) %>%
      mutate(
        end_window = begdesc + time_cutoff,
        row_id = row_number()
      )
    
    # Use non-equi join
    dt <- as.data.table(dives)
    result <- dt[dt, 
                 .(effort_timediving = sum(i.divetim),
                   effort_dives = .N,
                   effort_energy = sum(i.dive_energy)),
                 on = .(InstrumentSeq, begdesc > begdesc, begdesc <= end_window),
                 by = .EACHI
    ]
    
    return(as_tibble(result[, .(row_id = dt$row_id, effort_timediving, effort_dives, effort_energy)]))
  }
  
  #### 1e. efficiency function  ####
  calc_efficiency <- function(dives, dive_cutoff, pp_col) {
    dt <- as.data.table(dives)
    setorder(dt, InstrTrip, row_id)
    
    # date column for grouping (PP_in windows are constrained within trip-days)
    dt[, date := as.Date(begdesc)]
    
    dt[, `:=`(
      # PP_in: sum of prey captures in the dive_cutoff-wide window ending at the
      # previous dive
      PP_in = data.table::shift(
        frollapply(get(pp_col),
                   dive_cutoff,
                   sum,
                   na.rm = TRUE,
                   align = "right"),
        n = 1,
        fill = 0
      ),
      # timestamps bracketing the PP_in window, used to compute time_elapsed_in
      # after the shift, .end is the dive just before the current one
      # .start is dive_cutoff steps before that.
      .end   = data.table::shift(begdesc, n = 1),
      .start = data.table::shift(begdesc, n = dive_cutoff),
      # Cumulative prey captures since start of trip-day (for PP_out)
      cum_day = cumsum(replace(get(pp_col),
                               is.na(get(pp_col)),
                               0)),
      dives_elapsed_out = seq_len(.N),
      time_elapsed_out  = as.numeric(
        difftime(begdesc, first(begdesc), units = "mins")
      )
    ), by = .(InstrTrip, date)]
    
    # time_elapsed_in: minutes spanned by the PP_in window
    dt[, time_elapsed_in := as.numeric(
      difftime(.end, .start, units = "mins")
    ), by = .(InstrTrip, date)]
    
    # PP_out: cumulative captures before this dive, minus whats already in PP_in
    # shift(cum_day) gives captures up to (but not including) the current dive.
    dt[, PP_out := data.table::shift(cum_day, fill = 0) -
         fifelse(is.na(PP_in), 0, PP_in),
       by = .(InstrTrip, date)]
    
    # mask the first dive_cutoff rows per trip-day (not enough history for PP_in)
    dt[, c("PP_in", "time_elapsed_in") := lapply(
      .SD, function(x) replace(x, seq_len(.N) <= dive_cutoff, NA_real_)
    ), .SDcols = c("PP_in", "time_elapsed_in"),
    by = .(InstrTrip, date)]
    
    # drop temp columns
    dt[, c(".end", ".start") := NULL]
    
    return(as_tibble(dt[, .(row_id,
                            PP_in,
                            PP_out,
                            dives_elapsed_out,
                            time_elapsed_out,
                            time_elapsed_in)]))
  }
  
  #####
  
  # add row labels
  penguin_env_paired$row_id <- seq_len(nrow(penguin_env_paired))
  penguin_env_paired$InstrTrip <- paste(penguin_env_paired$InstrumentSeq,penguin_env_paired$TripNumber,sep='_')

  # calculate effort and efficiency
  effort <- calc_effort(penguin_env_paired,effort_cutoff) 
  efficiency <- calc_efficiency(penguin_env_paired,curr_succ_cutoff,PP_col)

  
  #### 1f. fit model ####
  model_df_prep <- penguin_env_paired %>%
    left_join(efficiency, by = "row_id") %>%
    left_join(effort, by = "row_id") %>%
    filter(time_elapsed_in<30) |> # exclude dives where current success is defined by >30min
    select(-row_id) |>
    mutate(sex = as.factor(Sex),
           chla_broad_scaled = as.vector(scale(chla_broad)),
           PP_in_scaled = as.vector(scale(PP_in)),
           PP_out_rate = as.vector(scale(PP_out/time_elapsed_out))) |>
    st_drop_geometry() |>
    dplyr::select(
      InstrTrip, DiveNumber,
      chla_broad_scaled,
      PP_in_scaled,
      effort_timediving,
      PP_out_rate,
      nChicks,
      sex) |>
    na.omit()

  model_df <- model_df_prep |>
    arrange(InstrTrip, DiveNumber) %>%              # ensure correct order
    group_by(InstrTrip) %>%                         # work within each trip
    filter(if(thinning_threshold == 1) TRUE else row_number() %% thinning_threshold == 1) %>%
    ungroup()
  
  mod <- glmmTMB(
    effort_timediving ~ PP_in_scaled*PP_out_rate + sex + nChicks + chla_broad_scaled + 
      (1|InstrTrip),
    family = Gamma(link = "log"),
    data = model_df
  )

  MI_list[[paste0("MI_",m)]] <- mod
}


#############################
# 2. Pool Models
#############################
ests <- lapply(MI_list, tidy, effects = "fixed", conf.int = FALSE)

# Combine into one df
ests_df <- bind_rows(ests, .id = "imp")

# apply Rubin's Rules
pooled <- ests_df %>%
  group_by(term) %>%
  summarise(
    m = n(),
    pooled_estimate = mean(estimate),
    
    # within-imputation variance (mean of squared SEs)
    Ubar = mean(std.error^2),
    
    # between-imputation variance
    B = var(estimate),
    
    # total variance
    total_var = Ubar + (1 + 1/m) * B,
    pooled_se = sqrt(total_var),
    
    # z / p-value
    z = pooled_estimate / pooled_se,
    p.value = 2 * pnorm(abs(z), lower.tail = FALSE)
  )


