##Code to clean up and check data before running models on it

set.seed(234324)

library(glmmTMB)
library(ggplot2)
library(GGally)
library(ggmap)
library(insight)
library(kableExtra)
library(lmerTest)
library(tidyverse)
library(ggforce)
library(huxtable)
library(AICcmodavg)
library(glmpathcr)
library(glmmTMB)
library(dplyr)
library(interactions)

#load in the data using the code that runs preliminary stats on it
source(file.path("c:","Users","cami","Documents","SE-Alaska_NABat","Code","read.detect.data.R"), local=TRUE, chdir=TRUE)
source(file.path("c:","Users","cami","Documents","SE-Alaska_NABat","Code","MV.read.detect.data.R"), local=TRUE, chdir=TRUE)

##------------------------------------------Formatting the dataframe for the modeling-----------------------------------------------------#
#change columns that should be character from number
bat.data$WeatherSource <- as.character(bat.data$WeatherSource)
#remove duplicate entries for sites (-dup and -sm)
bat.data <- bat.data[!grepl("-dup|-SM", bat.data$Quadrant, ignore.case = TRUE), ]
#separate transect from stationary data
bat.data <- bat.data[ bat.data$Quadrant != "Transects",]

###NEEED TO DO THISSS!!! In the Kaleidoscope outputs it looks like there are a lot of LACI, but they are all noise.
### Need to remove all of the LACI  outputs form 2019-2023 and only include one confirmed call from Yakitat (33486) in 2022. 

#Clean up transect data (comes from manual verified data rather tahn autoID)
#change columns that should be character from number
bat.data.MV$WeatherSource <- as.character(bat.data.MV$WeatherSource)
#remove duplicate entries for sites (-dup and -sm)
bat.data.MV <- bat.data.MV[!grepl("-dup|-SM", bat.data.MV$Quadrant, ignore.case = TRUE), ]
#separate transect from stationary data
bat.transect.data <- bat.data.MV[ bat.data.MV$Quadrant == "Transects",]
#add transect length (this should be done when putting otgether the raw data. Doign it here now becuase of timing, but try to be better next time)
# Grid-cell -> transect length lookup
tlength_lookup <- tibble(
  GRTS.Cell.ID = c(206, 5006, 18318, 34338, 47310, 14222,24782,30606, 31438,34254,50638,315563,33486),
  TLength      = c(48.88809967, 49.30189896, 38.09939957, 25.78700066, 47.06299973,71.96230316,27.53079987,
                   49.19760132,31.82169914,44.92910004,27.29159927,54.00579834,22.46789932)
)
# ---- pre-check: confirm the exact label used for transects ----
bat.transect.data %>% count(Quadrant)
# ---------------------------------------------------------------
bat.transect.data <- bat.transect.data %>%
  left_join(tlength_lookup, by = "GRTS.Cell.ID") %>%
  mutate(TLength = if_else(
    grepl("^transect", tolower(trimws(Quadrant))),  # transect rows only
    TLength,
    NA_real_
  ))

# ---- Check: transect rows that did NOT get a length ----
bat.transect.data %>%
  filter(grepl("^transect", tolower(trimws(Quadrant))), is.na(TLength)) %>%
  count(GRTS.Cell.ID, Quadrant)
# ------------------------------------------------------------



#create summarized detection counts for all nights for each species. 
#This code has different columns for single, all couplets included, and couplets include only when there is a single

##For Stationary
bat.data.long <- plyr::ldply(species.id, function(species, bat.data){
  bat.data2 <- bat.data
  bat.data$SpeciesGroup    <- species
  SpeciesGroupCountVars    <- bat.data.species[ grepl(species, bat.data.species)] 
  bat.data$SpeciesFullPool <- apply(bat.data[, SpeciesGroupCountVars],1,sum,na.rm=TRUE)
  bat.data$SpeciesSingleton<- as.vector(bat.data[, species, drop=TRUE])
  bat.data2[,SpeciesGroupCountVars] <- (diag( as.numeric(bat.data2[,species]>0))) %*% as.matrix(bat.data2[,SpeciesGroupCountVars] )
  bat.data$SpeciesPartPool <- apply(bat.data2[, SpeciesGroupCountVars],1,sum,na.rm=TRUE)
  bat.data
}, bat.data=bat.data)


#For transects
bat.transect.data.long <- plyr::ldply(species.id, function(species, bat.transect.data){
  bat.transect.data2 <- bat.transect.data
  bat.transect.data$SpeciesGroup    <- species
  SpeciesGroupCountVars    <- bat.data.species[ grepl(species, bat.data.species)]
  bat.transect.data$SpeciesFullPool <- apply(bat.transect.data[, SpeciesGroupCountVars],1,sum,na.rm=TRUE)
  bat.transect.data$SpeciesSingleton<- as.vector(bat.transect.data[, species, drop=TRUE])
  bat.transect.data2[,SpeciesGroupCountVars] <- (diag( as.numeric(bat.transect.data2[,species]>0))) %*% as.matrix(bat.transect.data2[,SpeciesGroupCountVars] )
  bat.transect.data$SpeciesPartPool <- apply(bat.transect.data2[, SpeciesGroupCountVars],1,sum,na.rm=TRUE)
  bat.transect.data
}, bat.transect.data=bat.transect.data)

#---------------------------------------Check and clean up data------------------------------------------------------------------------------------#
#__only includes grids that have at least 5 years of data___________________________________________________________________________________________
# If you'd rather treat the whole GRTS cell as the site, drop SiteName from the key.
bat.data.long <- bat.data.long %>%
  mutate(Site_key = paste0(GRTS.Cell.ID))
bat.transect.data.long <- bat.transect.data.long %>%
  mutate(Site_key = paste0(GRTS.Cell.ID))

# ---- Sampling effort per site ------------------------------------------------
site_effort <- bat.data.long %>%
  dplyr::group_by(Site_key) %>%
  dplyr::summarise(
    GRTS.Cell.ID = dplyr::first(GRTS.Cell.ID),
    SiteName     = paste(sort(unique(SiteName)), collapse = ", "),
    n_years      = dplyr::n_distinct(Year),
    years        = paste(sort(unique(Year)), collapse = ", "),
    n_nights     = dplyr::n_distinct(Night),
    .groups      = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_years))

# ---- Sampling effort transect ------------------------------------------------
site_effort_t <- bat.transect.data.long %>%
  dplyr::group_by(Site_key) %>%
  dplyr::summarise(
    GRTS.Cell.ID = dplyr::first(GRTS.Cell.ID),
    SiteName     = paste(sort(unique(SiteName)), collapse = ", "),
    n_years      = dplyr::n_distinct(Year),
    years        = paste(sort(unique(Year)), collapse = ", "),
    n_nights     = dplyr::n_distinct(Night),
    .groups      = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_years))

# ---- Keep sites with >= 5 years, then filter the data ------------------------
keep <- site_effort %>% filter(n_years >= 5) %>% pull(Site_key)
keept <- site_effort_t %>% filter(n_years >= 5) %>% pull(Site_key)

bat.data.5yr <- bat.data.long %>% filter(Site_key %in% keep)
bat.transect.data.5yr <- bat.transect.data.long%>% filter(Site_key %in% keept)

#----------For modelling will only use data from 2019 forward
bat.data.5yr <- bat.data.5yr[!(bat.data.5yr$Year %in% c(2015, 2017)), ]
bat.transect.data.5yr <- bat.transect.data.5yr[!(bat.transect.data.5yr$Year %in% c(2015, 2017)), ]

# ---- Sanity summary ----------------------------------------------------------
cat("Stationary Sites total: ", n_distinct(bat.data.long$Site_key), "\n",
    "Stationary Sites kept:  ", length(keep), "\n",
    "Stationary Rows before: ", nrow(bat.data.long), "\n",
    "Stationary Rows after:  ", nrow(bat.data.5yr), "\n", sep = "")

cat("Transect Routes total: ", n_distinct(bat.transect.data.long$Site_key), "\n",
    "Transect Routes kept:  ", length(keept), "\n",
    "Transect Rows before: ", nrow(bat.transect.data.long), "\n",
    "Transect Rows after:  ", nrow(bat.transect.data.5yr), "\n", sep = "")

#____Check to see if an exesive number of nights with 0 at the sites__________________________________________________________________________________
# 1. Collapse species-long frame to one row per site-night --------------------
#    Fullcount is the night's total detections (bats + NoID), constant per night.
night_lvl <- bat.data.5yr.2019 %>%                      
  mutate(Night = as.Date(Night)) %>%
  group_by(Site_key, GRTS.Cell.ID, SiteName, Year, Night) %>%
  summarise(total_det = first(Fullcount),
            n_rows    = n(),                        # sanity: # species-group rows/night
            .groups   = "drop")

night_lvl %>% count(n_rows)                         # expect one repeated value

# 2. Overall zero rate, per site ---------------------------------------------
night_lvl %>%
  group_by(Site_key) %>%
  summarise(n_nights = n(),
            n_zero   = sum(total_det == 0),
            pct_zero = round(100 * mean(total_det == 0), 1),
            .groups  = "drop") %>%
  arrange(desc(pct_zero)) %>%
  print(n = Inf)

# 3. One deployment per site-year; flag trailing zeros -----------------------
deploy <- night_lvl %>%
  arrange(Site_key, Year, Night) %>%
  group_by(Site_key, Year) %>%                       # <-- the whole year IS the deployment
  mutate(
    night_pos    = row_number(),                     # 1st, 2nd, 3rd… night of the year
    is_zero      = total_det == 0,
    last_det_pos = if (any(!is_zero)) max(night_pos[!is_zero]) else 0L,
    trailing0    = night_pos > last_det_pos           # nights after the final detection
  ) %>%
  ungroup()

remove <- deploy %>%
  filter(trailing0) %>%
  count(Site_key, Year, name = "n_trailing_zeros") %>%
  arrange(desc(n_trailing_zeros)) %>%
  print(n = Inf)
#confirm that the remove df has 0 this means there are no trailing 0 in this df. IF they do then you need to remove these

#_____Check for missing covariates_______________________________________________________________________________________________________________________
# ---- Candidate covariates: trim to what you'll actually model ---------------
site_covars  <- c("Distance.to.Clutter..m.","Percent.Clutter",
                  "Water.Nearby","WaterDist")
night_covars <- c("Nightly.Min.Temp","Nightly.Max.Temp",
                  "Nightly.Mean.RH","Nightly.Min.RH","Nightly.Max.RH",
                  "Nightly.Precipitation","Nightly.Mean.Windsp","jNight")

all_covars <- c(site_covars, night_covars)
all_covars <- all_covars[all_covars %in% names(bat.data.5yr)]   # keep only present cols

# ---- Collapse to one row per site-night (covariates constant within night) --
night_lvl <- bat.data.5yr.2019 %>%
  mutate(Night = as.Date(Night)) %>%
  group_by(Site_key, GRTS.Cell.ID, SiteName, Region, Year, Night) %>%
  summarise(across(all_of(all_covars), ~ first(.x)), .groups = "drop")

# ---- 1. Overall missingness per covariate -----------------------------------
miss_overall <- night_lvl %>%
  summarise(across(all_of(all_covars), ~ sum(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "covariate", values_to = "n_missing") %>%
  mutate(n_total     = nrow(night_lvl),
         pct_missing = round(100 * n_missing / n_total, 1)) %>%
  arrange(desc(pct_missing))

print(miss_overall, n = Inf)

# ---- 2. Missingness by Site (the important cut) ---------------------------
night_lvl %>% count(Site_key, name = "n_nights")     # denominators

miss_site <- night_lvl %>%
  group_by(Site_key) %>%
  summarise(across(all_of(all_covars), ~ round(100 * mean(is.na(.x)), 1)),
            .groups = "drop") %>%
  pivot_longer(-Site_key, names_to = "covariate", values_to = "pct_missing") %>%
  pivot_wider(names_from = Site_key, values_from = pct_missing)

print(miss_site, n = Inf)                        # cols = grids, cells = % missing


#write the df that will be used for the models
write.csv(bat.data.5yr, "Data/bat.data.5yr.csv")
write.csv(bat.transect.data.5yr, "Data/bat.transect.data.5yr.Aug20.csv")


#------------------------------Double checking the data (old code -- here fore reference but was not actually used)---------------------------------------------#

###Table of sites and effort by year
summary.report <- plyr::ddply(bat.data.5yr, c("Year"), function(x){
  n.sites = length(unique(x$GRTS.Cell.ID))
  n.site.nights= nrow( unique(x[,c("GRTS.Cell.ID","Night")]))
  data.frame(n.sites, n.site.nights)
})
summary.report

###plot to examine regularity of sampling by grts over years
ggplot(data=bat.data.5yr, aes(x=Year, y=as.character(GRTS.Cell.ID), color=Region))+
  geom_point()+
  ggtitle("Examining Regularity of Sampling by GRTS Over Years")+
  ylab("GRTS Cell ID")

###Examining the proportion of sites monitored in each region by year. 
#Check for broad biases between where the majority of our cells were sampled over time.
temp <- unique(bat.data.5yr[,c("GRTS.Cell.ID","Region","Year")])
round(prop.table( xtabs(~Year+Region, data=temp, exclude=NULL, na.action=na.pass), 1),2)

### find the calls that occur in a singlet and those in a couplet Skews like we see in LABO indicate that a majority of our labels come from couplets or singlets.
#this can help when determining how sure we feel about our data
 single.vs.double <- plyr::ldply(species.id, function(species, bat.data){
   bat.data$SpeciesGroup    <- species
   SpeciesGroupCountVars    <- bat.data.species[ grepl(species, bat.data.species)]
   bat.data$Doublet <- unlist(apply(bat.data[, SpeciesGroupCountVars],1,sum,na.rm=TRUE))
   bat.data$Doublet <- unlist(bat.data$Doublet - bat.data[, species])
   bat.data$Singlet <- bat.data[, species, drop=TRUE]
   bat.data
 }, bat.data=bat.data)
 
 plotdata <- single.vs.double
 nrow=2
 ncol=2
 plots.per.page <- ncol*nrow
 npages <- ceiling(length(species.id)/plots.per.page)
 if( length(species.id) %% plots.per.page != 0){
   # need to deal with a bug in ggforce where last page needs to be filled
   dummy <- data.frame( SpeciesGroup=paste("zzz",1:(plots.per.page-length(species.id) %% plots.per.page)), detect.night=NA)
   plotdata <- plyr::rbind.fill(plotdata, dummy)
 }

  plyr::l_ply(1:npages, function(page){
    myplot <- ggplot(data=single.vs.double, aes(x=Singlet, y=Doublet))+
      ggtitle("Detections in doublets vs detections in singlets")+
      geom_point()+
      facet_wrap_paginate(~SpeciesGroup, ncol=ncol, nrow=nrow, page=page, scales="free")+
      xlab("Detections in singlets")+ylab("Detections in doublets")
    print(myplot)
  })

###Table of raw count of detections by single id species
  #this can help in checking to see if any one species has a drastically different number of detections from year to year
spdetec <- as_data_frame(addmargins(xtabs(SpeciesSingleton~SpeciesGroup+Year, data=bat.data.5yr, exclude=NULL, na.action=na.pass),2))
spdetec$n <- formatC(spdetec$n, digits=0, format="f", width=6)
spdetec <- reshape2::dcast(spdetec, SpeciesGroup~Year, value.var="n")
singlet.total.counts <- spdetec
spdetec
  #plot this so it is easier to interpret
# Convert to long format for plotting
detection_long <- spdetec %>%
  pivot_longer(cols = -SpeciesGroup, 
               names_to = "Year", 
               values_to = "Detections")
detection_long <- detection_long[detection_long$Year!="Sum",] #get rid of the sum values for the quick plotting exrcise
detection_long$Detections <- as.numeric(detection_long$Detections)
detection_long$Year <- as.numeric(detection_long$Year)
# Create a faceted plot for better individual species visualization
facet_plot <- ggplot(detection_long, aes(x = Year, y = Detections)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue", size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed", alpha = 0.7) +
  facet_wrap(~SpeciesGroup, scales = "free_y", ncol = 4) +
  labs(title = "Individual Species Detection Trends with Linear Trend Lines",
       subtitle = "Blue = actual data, Red dashed = linear trend line with confidence interval",
       x = "Year",
       y = "Number of Detections") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(face = "bold"))
print(facet_plot)



