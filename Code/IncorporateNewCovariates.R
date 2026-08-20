#This code incorporates updated metadata for SE Alaska sites into my current df for modelling. 
#Need to run the ode twice, once for the KSP resutls and nother for the Manual vetted results

library(dplyr)

df <- read.csv("Data/All Year Activity by Night 2016-2025 AK-Aug13-KSP-temp.csv")
cov <- read.csv("Data/NABat_SEAK_UpdatedSiteCovariates-temp.csv")

cov_update <- cov %>%
  dplyr::transmute(
    join_grts = trimws(as.character(grtsID)),
    join_site = tolower(trimws(as.character(site_code))),
    Elevation_new       = as.numeric(site_elev),
    Percent.Clutter_new = as.numeric(clutter_pct),
    Habitat.Type_new    = as.character(habitat_category),
    WaterType_new       = as.character(water_type),
    strm_dist           = strm_dist,
    lake_dist           = lake_dist,
    harv_dist           = harv_dist,
    harv_year           = harv_year,
    road_dist           = road_dist
  )


# 2. Add matching normalized keys to df, join, then overwrite / append
df1 <- df %>%
  dplyr::mutate(
    join_grts = trimws(as.character(GRTS.Cell.ID)),
    join_site = tolower(trimws(sub("-.*$", "", as.character(SiteName)))) #gets rid fo the "-SM" for the duplicate sites
  ) %>%
  dplyr::left_join(cov_update, by = c("join_grts", "join_site")) %>%
  dplyr::mutate(
    Elevation       = dplyr::coalesce(Elevation_new,       as.numeric(Elevation)),
    Percent.Clutter = dplyr::coalesce(Percent.Clutter_new, as.numeric(Percent.Clutter)),
    Habitat.Type    = dplyr::coalesce(Habitat.Type_new,    as.character(Habitat.Type)),
    WaterType       = dplyr::coalesce(WaterType_new,       as.character(WaterType))
  ) %>%
  dplyr::select(-Elevation_new, -Percent.Clutter_new, -Habitat.Type_new, -WaterType_new,
                -join_grts, -join_site)

# 3. Update the distance to water based on the smallest number between current distance to lake vs distance to stream
df1 <- df1 %>%
  dplyr::mutate(
    WaterDist = pmin(strm_dist, lake_dist, na.rm = TRUE)
  )

write.csv(df1,"V:/ARU/SENSR-BAT/NABat/2025/AK/Analyzed/All Year Activity by Night 2016-2025 AK-Aug19-MV.csv")
