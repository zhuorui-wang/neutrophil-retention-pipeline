# neutrophil_retention_pipeline.R
# --------------------------------------------
#
# This script calculates neutrophil retention time and interaction dynamics
# with pre-neoplastic cells (PNCs) from 3D time-lapse imaging data.

# ---- Load Required Libraries ----
library(tidyverse)
library(dplyr)
library(data.table)


# ----Define Parameters ----
# Distance threshold for assigning neutrophil to PNC (in µm)
threshold <- 15

# ---- Step 1: Calculate 3D Distance Between Neutrophils and PNCs ----
calculate_distances <- function(neutrophil_df, pnc_df) {
  pnc_fixed <- pnc_df %>% filter(Time == 1)  # Assume PNCs are fixed in space
  all_distances <- list()
  
  time_points <- unique(neutrophil_df$Time)
  
  for (t in time_points) {
    neutrophils_t <- neutrophil_df %>% filter(Time == t)
    
    if (nrow(neutrophils_t) == 0) next
    
    distance_df <- inner_join(neutrophils_t, pnc_fixed, by = character(), suffix = c("_neu", "_pnc")) %>%
      mutate(Distance = sqrt((`Position X_neu` - `Position X_pnc`)^2 +
                               (`Position Y_neu` - `Position Y_pnc`)^2 +
                               (`Position Z_neu` - `Position Z_pnc`)^2),
             TimeFrame = t)
    
    all_distances[[as.character(t)]] <- distance_df
  }
  
  final_df <- bind_rows(all_distances)
  return(final_df)
}

# ---- Step 2: Assign Closest PNC (Within Threshold) to Each Neutrophil ----
assign_closest_pnc <- function(distance_df) {
  closest_pnc <- distance_df %>%
    filter(Distance <= threshold) %>%
    group_by(TimeFrame, TrackID) %>%
    slice_min(Distance, with_ties = FALSE) %>%
    ungroup()
  
  return(closest_pnc)
}

# ---- Step 3: Calculate Retention Time for Neutrophil-PNC Interactions ----
calculate_retention_time <- function(closest_pnc_df) {
  retention_df <- closest_pnc_df %>%
    arrange(TrackID, TimeFrame) %>%
    group_by(TrackID) %>%
    mutate(Interaction_Group = rleid(ID_pnc)) %>%
    ungroup() %>%
    group_by(TrackID, ID_pnc, Interaction_Group) %>%
    summarise(Start_Time = min(TimeFrame),
              End_Time = max(TimeFrame),
              Duration_Frames = End_Time - Start_Time + 1,
              .groups = "drop") %>%
    arrange(TrackID, Start_Time)
  
  return(retention_df)
}

# ---- Example Usage ----
# distance_df <- calculate_distances(neu_df, pnc_df)
# closest_pnc_df <- assign_closest_pnc(distance_df)
# retention_time_df <- calculate_retention_time(closest_pnc_df)

# ---- Optional: Convert Frame Count to Seconds ----
# frame_rate <- 20  # e.g., seconds per frame
# retention_time_df <- retention_time_df %>%
#   mutate(Retention_Time_sec = Duration_Frames * frame_rate)
