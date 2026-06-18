# Code for Johansson et al. 
# trip_scale_model

# 0. Load Libraries and Data
# 1. Multiple Imputation Loop
## 1a. set parameters
## 1b. aggregate and filter trip data
## 1c. generate available trips
## 1d. generate KDEs for used and available trips
## 1e. calculate overlap
## 1f. calculate trip-specific chla
## 1g. join data
## 1h. fit model
# 2. Pool Models

#############################
# 0. Load Libraries and Data
#############################
library(tidyverse)
library(rnaturalearth)
library(sf)
library(terra)
library(amt)
library(glmmTMB)
library(survival)
library(broom.mixed)
library(profvis)
library(ggthemes)
library(patchwork)
library(broom)

setwd("O:/PERSONAL/Erik J/Chapter3") 

load("byTrip_df.RData")
load("byDive_PP_df.RData")
load("sim_dive_locs.RData")
load("chla.RData")
load("internal_states_df.RData")

setwd("O:/PERSONAL/Erik J") 

#############################
# 1. Multiple Imputation Loop
#############################
MI_list <- list()
model_df_list <- list()

for(m in 1:10){
  #### 1a. set parameters ####
  n_avail <- 20
  imputation_id <- m
  
  #### 1b. aggregate and filter trip data ####
  # add row identifiers to byTrip and byDive
  byTrip_df <- byTrip_df |>
    mutate(InstrTrip = paste(InstrumentSeq,TripNumber,sep='_'))
  
  byDive_PP_df <- byDive_PP_df |>
    mutate(InstrTrip = paste(InstrumentSeq,TripNumber,sep='_'),
           InstrTripDive = paste(InstrumentSeq,TripNumber,DiveNumber,sep='_'))
  
  # filter sim_dive_locs based on trip duration, dive depth, and dive MCP area
  sim_dive_locs_f <- sim_dive_locs |>
    mutate(InstrTrip = paste(InstrumentSeq,TripNumber,sep='_'),
           InstrTripDive = paste(InstrTrip,DiveNumber,sep='_')) |>
    filter(InstrTrip %in% 
             byTrip_df$InstrTrip[byTrip_df$Duration>24 & byTrip_df$Duration<96 & !is.na(byTrip_df$TripNumber)],
           InstrTripDive %in%
             byDive_PP_df$InstrTripDive[byDive_PP_df$maxdep>5 & byDive_PP_df$MCP_Area_km2<25])
  
  #### 1c. generate available trips ####
  # function to choose a random trip and rotate it randomly
  rotate_points <- function(all_locs,
                            land,
                            colonyLat = -44.04527,
                            colonyLon = -65.22437,
                            utm_crs = 32720){  # Add UTM zone as parameter
    # get land polygons
    trys <- 0
    repeat {
      trys <- trys + 1
      random_trip <- sample(unique(all_locs$InstrTrip), 1)
      dive_points <- st_as_sf(all_locs |>
                                filter(InstrTrip == random_trip,
                                       sim_id == imputation_id),
                              coords = c("mu.x", "mu.y"),
                              crs = utm_crs)  
      
      # transform to lat/long first, then to azimuthal equidistant
      dive_points_ll <- st_transform(dive_points, crs = 4326)
      
      # project to planar coordinates centered on colony
      proj <- sprintf("+proj=aeqd +lat_0=%f +lon_0=%f", colonyLat, colonyLon)
      dive_points_proj <- st_transform(dive_points_ll, crs = proj)
      
      # calculate offsets from colony (which is at origin in this projection)
      coords <- st_coordinates(dive_points_proj)
      theta <- runif(1, 0, 2*pi)
      
      # 2D rotation matrix
      R <- matrix(c(cos(theta), -sin(theta),
                    sin(theta),  cos(theta)), nrow = 2)
      
      # rotate coords
      new_coords <- t(R %*% t(coords))
      df_rotated <- data.frame(
        X = new_coords[,1],
        Y = new_coords[,2]
      )
      
      # make sf MULTIPOINT
      trip_rotated_proj <- st_as_sf(df_rotated, coords=c("X","Y"), crs=proj)
      trip_rotated_ll <- st_transform(trip_rotated_proj, crs=4326)
      
      # check overlap with land
      intersections <- st_intersects(trip_rotated_ll, land)
      if(all(lengths(intersections) == 0)) {
        trip_rotated_utm <- st_transform(trip_rotated_ll, utm_crs)  # Back to original UTM
        return(trip_rotated_utm)
      }
    }
  }
  
  sa_land <- ne_countries(scale='large', continent = 'south america', returnclass='sf')
  sa_land <- st_crop(sa_land, xmin = -68, xmax = -61, ymin = -46, ymax = -42)
  
  # chunk to exclude first trip and singular trips (dont need avail trips for these)
  first_trips <-
    sim_dive_locs_f |> 
    group_by(InstrumentSeq) |> 
    summarise(first_trip = min(InstrTrip))
  
  trips_to_simulate <- unique(sim_dive_locs_f$InstrTrip)[!unique(sim_dive_locs_f$InstrTrip) %in% first_trips$first_trip]
  
  avail_list <- vector("list",length(trips_to_simulate) * n_avail)
  k <- 1
  for(i in trips_to_simulate){
    # control for trip distance when generating random uds
    trips_within_10 <- byTrip_df |> filter(abs(maxDistance - byTrip_df$maxDistance[byTrip_df$InstrTrip==i]) <=10)
    points_from_trips_within_10 <- sim_dive_locs_f |>
      filter(InstrTrip%in%trips_within_10$InstrTrip)
    
    for(j in c(1:n_avail)){
      avail_list[[k]] <- rotate_points(all_locs=points_from_trips_within_10,land=sa_land) |>
        mutate(InstrTrip = i,
               avail_id = j)
      k <- k+1
    }
  }
  avail_utm <- bind_rows(avail_list)
  
  
  
  used_utm <- st_as_sf(sim_dive_locs_f|>filter(sim_id == imputation_id),
                       coords = c("mu.x", "mu.y"),
                       crs = 32720)
  used_utm <- used_utm |> 
    mutate(avail_id = NA) |>
    dplyr::select(InstrTrip,avail_id,geometry)
  
  
  #### 1d. generate KDEs for used and available trips #### 
  all_pts <- rbind(avail_utm, used_utm)  
  bb <- st_bbox(all_pts)
  
  rast_template <- rast(
    xmin = bb["xmin"]-10000,
    xmax = bb["xmax"]+10000,
    ymin = bb["ymin"]-10000,
    ymax = bb["ymax"]+10000,
    resolution = 2000,     
    crs = "EPSG:32720"
  )
  
  all_pts <- all_pts %>%
    mutate(trip_uid = ifelse(is.na(avail_id),
                             paste0(InstrTrip, "_used_0"),
                             paste0(InstrTrip, "_avail_", avail_id)))
  ud_list <- all_pts %>%
    split(.$trip_uid) %>%
    map(~{
      coords <- st_coordinates(.x)
      df <- data.frame(x = coords[,1], y = coords[,2])
      track <- make_track(df, x, y, crs = 32720)
      hr_kde(track, trast = rast_template, h=4000, levels = c(0.95))
    })
  
  
  #### 1e. calculate overlap ####
  # Parse trip names
  trip_info <- tibble(trip_uid = names(ud_list)) %>%
    tidyr::extract(trip_uid, 
                   into = c("InstrSeq", "TripNum", "Type", "AvailID"),
                   regex = "([0-9]+)_([0-9]+)_(used|avail)_([0-9]*)",
                   remove = FALSE) %>%
    mutate(TripNum = as.integer(TripNum))
  
  # Get all trips with TripNum > 1 (exclude first trip of each InstrSeq)
  current_trips <- trip_info %>%
    filter(TripNum > 1)
  
  # For each current trip, find the previous USED trip and calculate overlap
  overlap_results <- current_trips %>%
    rowwise() %>%
    mutate(
      # Find previous used trip (allows for non-sequential trip numbers)
      prev_used_trip = {
        prev_trips <- trip_info %>%
          filter(InstrSeq == .env$InstrSeq, 
                 TripNum < .env$TripNum,  # Get all previous trips
                 Type == "used") %>%
          arrange(desc(TripNum)) %>%      # Sort descending to get most recent first
          pull(trip_uid)
        
        if(length(prev_trips) > 0) prev_trips[1] else NA_character_
      }
    ) %>%
    ungroup() %>%
    filter(!is.na(prev_used_trip)) %>%  # Only keep if previous used trip exists in ud_list
    rowwise() %>%
    mutate(
      used = ifelse(Type == "used", 1, 0),
      ba_overlap = hr_overlap(ud_list[[trip_uid]], 
                              ud_list[[prev_used_trip]], 
                              type = "ba")$overlap
    ) %>%
    ungroup() %>%
    select(trip = trip_uid, used, compared_trip = prev_used_trip, ba_overlap)
  
  #### 1f. calculate trip-specific chla ####
  chla_sf <- chla %>%
    mutate(time = as_date(time),
           chlor_a = log(chlor_a + 0.0001)) %>%
    drop_na(longitude, latitude, chlor_a) %>%  # Remove NAs upfront
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
  
  # extract all isopleths
  all_isos <- map_dfr(names(ud_list), function(uid) {
    tryCatch({
      iso <- hr_isopleths(ud_list[[uid]])
      st_crs(iso) <- 32720
      iso_fixed <- st_buffer(iso, dist = 0)
      iso_latlon <- st_transform(iso_fixed, crs = 4326)
      iso_latlon$trip_uid <- uid
      iso_latlon
    })
  })
  
  # join trip dates
  ud_dates <- tibble(trip_uid = names(ud_list)) %>%
    mutate(InstrTrip = str_extract(trip_uid, "^[0-9]+_[0-9]+")) %>%
    left_join(
      byTrip_df %>% select(InstrTrip, StartTime, EndTime), 
      by = "InstrTrip"
    ) %>%
    mutate(
      startDate = as_date(StartTime),
      endDate = as_date(EndTime)
    )
  
  # spatial join and calculate means
  byTrip_envvar_df <- all_isos %>%
    left_join(ud_dates, by = "trip_uid") %>%
    st_join(chla_sf, join = st_intersects) %>%
    st_drop_geometry() %>%
    filter(time >= startDate & time <= endDate) %>%
    group_by(trip_uid) %>%
    summarise(mean_chla = mean(chlor_a, na.rm = TRUE),
              chla_count = sum(!is.na(chlor_a)), .groups = "drop") |>
    # add back trips with no chlorophyll data
    right_join(
      ud_dates %>% select(trip_uid),
      by = c("trip_uid")
    ) %>%
    replace_na(list(mean_chla = 0))

  #### 1g. join data ####
  model_df_prep <- overlap_results |>
    mutate(InstrumentSeq = as.numeric(substr(trip,1,4)),
           InstrTrip = substr(trip,1,6),
           prev_InstrTrip = substr(compared_trip,1,6)) |> # assumes no more than 9 trips per deployment
    rowwise() |>
    mutate(sex = internal_states_df$Sex[internal_states_df$InstrumentSeq==InstrumentSeq],
           n_chicks = internal_states_df$nChicks[internal_states_df$InstrumentSeq==InstrumentSeq],
           prev_tripEnergy = byTrip_df$energy_J[byTrip_df$InstrTrip==prev_InstrTrip],
           chla = byTrip_envvar_df$mean_chla[byTrip_envvar_df$trip_uid==trip]
           ) |>
    ungroup() |>
    rename(prev_trip = compared_trip,
           stratum = InstrTrip,
           overlap = ba_overlap)
  
  model_df <- model_df_prep |>
    mutate(
      overlap_scaled = scale(overlap)[,1],
      prev_tripEnergy_scaled = scale(prev_tripEnergy)[,1],
      chla_scaled = scale(chla)[,1],
      sex = factor(sex)) 
  
  
  #### 1h. fit model ####
  mod <- clogit(
    used ~ overlap_scaled + chla_scaled + 
      overlap_scaled:prev_tripEnergy_scaled + 
      overlap_scaled:n_chicks + 
      overlap_scaled:sex + 
      strata(stratum),
    data = model_df
  )

  MI_list[[paste0("MI_",m)]] <- mod
  model_df_list[[paste0("MI_",m)]] <- model_df
  

  save_these <- c("byDive_PP_df","byTrip_df","chla","internal_states_df",
                  "sim_dive_locs","MI_list",
                  "model_df_list")
  rm(list = setdiff(ls(), save_these))
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
    
    # total variance (Rubin's Rules)
    total_var = Ubar + (1 + 1/m) * B, # 
    pooled_se = sqrt(total_var),
    
    # Z / p-value
    z = pooled_estimate / pooled_se,
    p.value = 2 * pnorm(abs(z), lower.tail = FALSE)
  )

