# Data and code for: 'Outcomes of Previous Experiences Shape Behavioral Decision-Making Across Scales'

**Corresponding author:** Erik Johansson, ejohanss@uw.edu

**Repository DOI:**  10.5281/zenodo.20751054

## Overview

This repository contains the data and R code used to examine how prior foraging success shapes broad- and fine-scale foraging decisions in *Spheniscus magellanicus*, using GPS, TDR, and accelerometry biologging data collected at Punta Tombo during 2024-2025. 

## Repository structure

| File | Type | Description |
|---|---|---|
| `chla.RData` | Data | Satellite-derived chlorophyll-a concentration data, extracted from the NOAA ERDDAP data server [ID: pmlEsaCCI60OceanColorDaily] |
| `internal_states_df.RData` | Data | Individual-level records of sex and number of chicks |
| `byDive_df.RData` | Data | Dive-level dataset, one row per dive. |
| `byTrip_df.RData` | Data | Trip-level dataset, one row per foraging trip. |
| `sim_dive_locs.RData` | Data | Simulated dive locations, generated from crawl model fit to GPS data |
| `trip_scale_model.R` | Code | Fits the trip-scale conditional logistic regression model (with multiple imputation) |
| `dive_scale_model.R` | Code | Fits the dive-scale mixed effects linear model (with multiple imputation) |

## Data file details

### `byTrip_df.RData`

One row per foraging trip. n = 254 rows.

| Column | Type | Description |
|---|---|---|
| `InstrumentSeq` | numeric | unique identifier of tag deployment |
| `TripNumber` | numeric | identifier of foraging trip, unique within each tag deployment |
| `maxDistance` | numeric | maximum distance from colony reached during forgaging trip (km) |
| `Duration` | numeric | duration of foraging trip (hrs) |
| `StartTime` | POSIXct | timestamp of first at-sea location after leaving colony |
| `EndTime` | POSIXct | timestamp of last at-sea location before returning to colony |
| `energy_J` | numeric | energetic expenditure summed across foraging trip (J) - calculated from accelerometry data |
| `InstrTrip` | character | concatenated 'InstrumentSeq_TripNumber' identifier |

### `byDive_df.RData`

One row per dive. n = 260903 rows.

| Column | Type | Description |
|---|---|---|
| `begdesc` | POSIXct | timestamp of first underwater record |
| `divetim` | numeric | duration of underwater period (s) |
| `maxdep` | numeric | maximum depth below surface reached (m) |
| `DiveNumber` | numeric | identifier of dive, unique within each tag deployment |
| `energy_J` | numeric | energetic expenditure summed across dive period (J) - calculated from accelerometry data |
| `TripNumber` | numeric | identifier of foraging trip, unique within each tag deployment |
| `MCP_Area_km2` | numeric | uncertainty in dive location, measured as the area of the 90% MCP of simulated dive locations (km2) |
| `InstrumentSeq` | numeric | unique identifier of tag deployment |
| `cluster16_PP95_TF` | logical | whether the dive included prey pursuit, based on Gaussian mixture model (T/F) |
| `InstrTrip` | character | concatenated 'InstrumentSeq_TripNumber' identifier |
| `InstrTripDive` | character | concatenated 'InstrumentSeq_TripNumber_DiveNumber' identifier |

### `chla.RData`

One row per gridcell per day. n = 282492 rows.

| Column | Type | Description |
|---|---|---|
| `longitude` | numeric | longitude of gridcell |
| `latitude` | numeric | latitude of gridcell |
| `time` | character | day of chla record |
| `chlor_a` | numeric | extracted chla_a value |
| `chlor_a_log10_bias` | numeric | extracted chla_a bias |
| `chla_corrected` | numeric | corrected chlor_a value, using chla_corrected = 10^(log10(chla)+chla_bias) |

### `internal_states_df.RData`

One row per individual. n = 70 rows.

| Column | Type | Description |
|---|---|---|
| `PenguinSeq` | numeric | unique identifier of individual |
| `InstrumentSeq` | numeric | unique identifier of tag deployment |
| `nChicks` | numeric | number of chicks being provisioned by individual (1 or 2) |
| `Sex` | character | M/F |

### `sim_dive_locs.RData`

One row per simulated dive location (20 per dive). n = 5218060 rows. 

| Column | Type | Description |
|---|---|---|
| `InstrumentSeq` | numeric | unique identifier of tag deployment |
| `TripNumber` | character | identifier of foraging trip, unique within each tag deployment |
| `DiveNumber` | numeric | identifier of dive, unique within each tag deployment |
| `GMTDateTime` | POSIXct | timestamp of first underwater record |
| `mu.x` | numeric | simulated latitude (UTM) |
| `mu.y` | numeric | simulated longitude (UTM) |
| `sim_id` | numeric | simulation identifier (1-20) |

