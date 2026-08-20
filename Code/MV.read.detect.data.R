# read the bat detection data base
library(tidyverse)
library(dplyr)
library(purrr)

#speciesdetails <- readxl::read_excel(file.path("..","Data","Grid Cell Comments.xlsx"),
#                               sheet="SpeciesDetails",
#                               .name_repair='universal')

#gridcomments <- readxl::read_excel(file.path("..","Data","Grid Cell Comments.xlsx"),
#                                       sheet="2023",
#                                       .name_repair='universal')

#quaddetails <- readxl::read_excel(file.path("..","Data","Grid Cell Comments.xlsx"),
#                                   sheet="Quadrants",
#                                   .name_repair='universal')

#bec.data <- readxl::read_excel(file.path("..","Data","NABat Master Bec Zone Union.xlsx"),
#                               sheet="Bec Zone Subzone Areas",
#                               .name_repair='universal')

#bec.sp.data <- readxl::read_excel(file.path("..","Data","NABat Master Bec Zone Union.xlsx"),
#                               sheet="Species by ZoneSubzone",
#                               .name_repair='universal')


#bat.data.MV <- readxl::read_excel("C:/Users/cami/Documents/BC-NABat/Data/All Year Activity by Night 2016-2024 BC V2 - TEMP.xlsx",
# bat.data.MV <- readxl::read_excel(file.path("..","Data","All Year Activity by Night 2016-2023 Manual.xlsx"),
                             ###Toggle between these two lines to switch between all files and only files that passed auto ID
                             # sheet="OccupancyNoID",
                            # sheet="Occupancy",
                             #.name_repair='universal')

bat.data.MV <- readxl::read_excel("V:/ARU/SENSR-BAT/NABat/2025/AK/Analyzed/All Year Activity by Night 2016-2025 AK-Aug19.xlsx",
                               sheet="Occupancy",
                               .name_repair='universal') 

#time.data <- readxl::read_excel(file.path("..","Data","All Year Activity Time 2016-2023 AK and BC.xlsx"),
#                               sheet="Activity Time",
#                               .name_repair='universal')
#change Night to an actual recognisable date
bat.data.MV <- bat.data.MV %>%
  mutate(Night = as.Date(Night))

# Define columns that should stay as character/text
keep_as_text <- c("Quadrant", "Night", "Name", "SiteName", "Habitat.Type", 
                  "Feature.Sampled", "WaterType", "Region","srise", "sset",
                  "mrise",	"mset",	"moonr2",	"moons2","WeatherSource"
                  )
# Convert everything else to numeric
for(col in names(bat.data.MV)) {
  if(!col %in% keep_as_text) {
    bat.data.MV[[col]] <- as.numeric(bat.data.MV[[col]])
  }
}

bat.data.MV$rnum <- 1:nrow(bat.data.MV)
names(bat.data.MV)
bat.data.MV$Night <- as.Date(bat.data.MV$Night)
bat.data.MV$Year <- lubridate::year(bat.data.MV$Night)

# check year of data
xtabs(~Region+Year, data=bat.data.MV, exclude=NULL, na.action=na.pass)
length(unique(bat.data.MV$GRTS.Cell.ID))

#-----------------------
#Only keep data that was indetified to species or couplet level for the species that we know to be present in Alaska

targets <- c("LACI", "LANO", "MYCA", "MYEV", "MYLU", "MYVO", "MYYU")

# Split a column name into consecutive 4-char codes.
# Returns the codes only if the whole name is a clean run of 4-letter uppercase codes,
# otherwise character(0) so it contributes to nothing.
split_codes <- function(nm) {
  if (nchar(nm) == 0 || nchar(nm) %% 4 != 0) return(character(0))
  chunks <- substring(nm, seq(1, nchar(nm), by = 4), seq(4, nchar(nm), by = 4))
  if (all(grepl("^[A-Z]{4}$", chunks))) chunks else character(0)
}

# Columns that sit strictly between Night and Year
night_idx  <- which(names(bat.data.MV) == "Night")
year_idx   <- which(names(bat.data.MV) == "Year")
block_cols <- names(bat.data.MV)[(night_idx + 1):(year_idx - 1)]

# Map each target species to every column whose name contains it
contrib <- setNames(
  lapply(targets, function(sp) {
    block_cols[map_lgl(block_cols, ~ sp %in% split_codes(.x))]
  }),
  targets
)

# ---- diagnostic preview: confirm the fold-in before committing ----
print(contrib)
# -------------------------------------------------------------------

# Build the 7 collapsed species columns (na.rm = TRUE so a couplet NA
# does not wipe out a known count in the pure-species column)
new_species <- map(targets, function(sp) {
  cols <- contrib[[sp]]
  if (length(cols) == 0) return(rep(0, nrow(bat.data.MV)))
  rowSums(as.matrix(bat.data.MV[, cols, drop = FALSE]), na.rm = TRUE)
}) |>
  set_names(targets) |>
  as_tibble()

# Reassemble: everything up to & including Night, the 7 species, then Year onward
bat.data.MV <- bind_cols(
  bat.data.MV[, 1:night_idx, drop = FALSE],
  new_species,
  bat.data.MV[, year_idx:ncol(bat.data.MV), drop = FALSE]
)

#-----------------------
# check the dates on the data
xtabs(~GRTS.Cell.ID+Year, data=bat.data.MV, exclude=NULL, na.action=na.pass)



# check dates when data are collected by year
# email of 2021-02-16. Truncate at the end of August which is about Julian Day 240
bat.data.MV$jNight <- lubridate::yday(bat.data.MV$Night)
myplot <- ggplot(data=bat.data.MV[ bat.data.MV$Year > 2010,], aes(x=jNight, y=Year))+
  ggtitle("Julian nights of data collection")+
  geom_point(position=position_jitter(h=.1, w=.1))
myplot
ggsave(myplot, file=file.path("c:","Users","cami", "Documents","SE-Alaska_NABat","Figures","Summary", "jNight-vs-year.png"), h=4,w=6, units="in", dpi=300)

xtabs(~jNight, data=bat.data.MV, exclude=NULL, na.action=na.pass)
#bat.data.MV1 <- bat.data.MV #keep a copy with all of the dates just in case. But create a new df that only has data before julian 240 (End of August)
bat.data.MV <- bat.data.MV[ bat.data.MV$jNight < 240,]
bat.data.MV$jNightM15 <- bat.data.MV$jNight - lubridate::yday("2020-05-14")

#-----------------------
# check lat long
ggplot(data=bat.data.MV, aes(x=Long, y=Lat,color=as.factor(Year)))+
  ggtitle("Where collected")+
  geom_point( position=position_jitter(h=.1, w=.1), shape=1)


#---------------
# check quadrant
xtabs(~Year+Quadrant, data=bat.data.MV, exclude=NULL, na.action=na.pass)
bat.data.MV$Quadrant <- paste0(toupper(substr(bat.data.MV$Quadrant,1,1)), substring(bat.data.MV$Quadrant,2))
xtabs(~Year+Quadrant, data=bat.data.MV, exclude=NULL, na.action=na.pass)

xtabs(~GRTS.Cell.ID+Quadrant, data=bat.data.MV, exclude=NULL, na.action=na.pass)
#-----------------------
# check habitat type
xtabs(~Habitat.Type+Year, data=bat.data.MV, exclude=NULL, na.action=na.pass)

# extract the first (primary) habitat type
bat.data.MV$Habitat2 <- trimws(substr(bat.data.MV$Habitat.Type,1, -1+regexpr(",", paste0(bat.data.MV$Habitat.Type,","), fixed=TRUE)))
xtabs(~Habitat2+Year, data=bat.data.MV, exclude=NULL, na.action=na.pass)


#-----------------------
# email 2021-02-16 Drop from analysis
#xtabs(~Feature.Sampled+Year, data=bat.data.MV, exclude=NULL, na.action=na.pass)

#----------------------
xtabs(~Distance.to.Clutter..m.,data=bat.data.MV, exclude=NULL, na.action=na.pass)
bat.data.MV$Distance.to.Clutter..m. <- as.numeric(bat.data.MV$Distance.to.Clutter..m.)
bat.data.MV$Distance.to.Clutter..m.[bat.data.MV$Distance.to.Clutter..m. > 80] <- 80
xtabs(~Distance.to.Clutter..m.,data=bat.data.MV, exclude=NULL, na.action=na.pass)

#-----------------------
xtabs(~Percent.Clutter,data=bat.data.MV, exclude=NULL, na.action=na.pass)
bat.data.MV$Percent.Clutter <- as.numeric(bat.data.MV$Percent.Clutter)
xtabs(~Quadrant+Percent.Clutter,data=bat.data.MV, exclude=NULL, na.action=na.pass)


#-----------------------
# confirm that 1=yes
# Email 2021-02-16. Yes, 1=Yes
xtabs(~Quadrant+Water.Nearby,data=bat.data.MV, exclude=NULL, na.action=na.pass)

#-----------------------
# Convert second of daylength to decimal hours
# 
# bat.data.MV$dayl <- bat.data.MV$dayl / 3600
# bat.data.MV$dayl.r <- round(bat.data.MV$dayl)
# xtabs(~dayl,data=bat.data.MV, exclude=NULL, na.action=na.pass)
# xtabs(~dayl.r,data=bat.data.MV, exclude=NULL, na.action=na.pass)
# sum(is.na(bat.data.MV$dayl))
# 
# #-----------------------
# # Seem unlikely that prcp is measured to 1/100 mm
# bat.data.MV$prcp.r <- round(bat.data.MV$prcp)
# xtabs(~prcp,data=bat.data.MV, exclude=NULL, na.action=na.pass)
# xtabs(~prcp.r,data=bat.data.MV, exclude=NULL, na.action=na.pass)

# added a column for nightly precipitation???NOt sure what negative values mean. Most of these are also missing so not useful
xtabs(~Nightly.Precipitation, data=bat.data.MV, exclude=NULL, na.action=na.pass)
bat.data.MV$Nightly.Precipitation <- round(bat.data.MV$Nightly.Precipitation,.1)
xtabs(~Nightly.Precipitation+Year, data=bat.data.MV, exclude=NULL, na.action=na.pass)

#create nightly.mean.temp?



# #-----------------------
# xtabs(~srad,data=bat.data.MV, exclude=NULL, na.action=na.pass)
# sum(is.na(bat.data.MV$srad))
# hist(bat.data.MV$srad)
# 
# #-----------------------
# xtabs(~tmax,data=bat.data.MV, exclude=NULL, na.action=na.pass)
# sum(is.na(bat.data.MV$tmax))
# hist(bat.data.MV$tmax)
# bat.data.MV[ bat.data.MV$tmax < 10 & !is.na(bat.data.MV$tmax),c("Year","GRTS.Cell.ID","Lat","Long","Night","tmax","tmin")]

bat.data.MV$Nightly.Max.Temp <- as.numeric(bat.data.MV$Nightly.Max.Temp)
sum(is.na(bat.data.MV$Nightly.Max.Temp))
xtabs(~Nightly.Max.Temp,data=bat.data.MV, exclude=NULL, na.action=na.pass)
hist(bat.data.MV$Nightly.Max.Temp)
bat.data.MV[ bat.data.MV$Nightly.Max.Temp < 10 & !is.na(bat.data.MV$Nightly.Max.Temp),c("Year","GRTS.Cell.ID","Lat","Long","Night","Nightly.Max.Temp","Nightly.Min.Temp")]

# #-----------------------
# xtabs(~tmin,data=bat.data.MV, exclude=NULL, na.action=na.pass)
# sum(is.na(bat.data.MV$tmin))
# hist(bat.data.MV$tmin)

bat.data.MV$Nightly.Min.Temp <- as.numeric(bat.data.MV$Nightly.Min.Temp)
sum(is.na(bat.data.MV$Nightly.Min.Temp))
xtabs(~Nightly.Min.Temp,data=bat.data.MV, exclude=NULL, na.action=na.pass)
hist(bat.data.MV$Nightly.Min.Temp)


# #-----------------------
# xtabs(~dayl,data=bat.data.MV, exclude=NULL, na.action=na.pass)
# sum(is.na(bat.data.MV$dayl))
# hist(bat.data.MV$dayl)

# #-----------------------
# xtabs(~vp,data=bat.data.MV, exclude=NULL, na.action=na.pass)
# sum(is.na(bat.data.MV$vp))
# hist(bat.data.MV$vp)
# 
# #-----------------------
# xtabs(~srad,data=bat.data.MV, exclude=NULL, na.action=na.pass)
# sum(is.na(bat.data.MV$srad))
# hist(bat.data.MV$srad)

#-----------------------
bat.data.MV$Nightly.Min.RH <- as.numeric(bat.data.MV$Nightly.Min.RH )
sum(is.na(bat.data.MV$Nightly.Min.RH ))
xtabs(~Nightly.Min.RH ,data=bat.data.MV, exclude=NULL, na.action=na.pass)
hist(bat.data.MV$Nightly.Min.RH )

#-----------------------
bat.data.MV$Nightly.Max.RH <- as.numeric(bat.data.MV$Nightly.Max.RH)
sum(is.na(bat.data.MV$Nightly.Max.RH))
xtabs(~Nightly.Max.RH,data=bat.data.MV, exclude=NULL, na.action=na.pass)
hist(bat.data.MV$Nightly.Max.RH)

#-----------------------
bat.data.MV$Nightly.Mean.Windsp <- as.numeric(bat.data.MV$Nightly.Mean.Windsp)
sum(is.na(bat.data.MV$Nightly.Mean.Windsp))
xtabs(~Nightly.Mean.Windsp,data=bat.data.MV, exclude=NULL, na.action=na.pass)
hist(bat.data.MV$Nightly.Mean.Windsp)

#-----------------------
# Moon rise is read in a Excel decimal days
# See https://stackoverflow.com/questions/19172632/converting-excel-datetime-serial-number-to-r-datetime

#bat.data.MV$mrise <-  as.POSIXct(as.numeric(bat.data.MV$mrise) * (60*60*24), origin="1899-12-30", tz="UTC")
#bat.data.MV$mset  <-  as.POSIXct(as.numeric(bat.data.MV$mset ) * (60*60*24), origin="1899-12-30", tz="UTC")
#bat.data.MV$srise <-  as.POSIXct(as.numeric(bat.data.MV$srise) * (60*60*24), origin="1899-12-30", tz="UTC")
#bat.data.MV$sset  <-  as.POSIXct(as.numeric(bat.data.MV$sset ) * (60*60*24), origin="1899-12-30", tz="UTC")

bat.data.MV$mtime <- as.numeric(bat.data.MV$mtime)
sum(is.na(bat.data.MV$mtime))
xtabs(~mtime,data=bat.data.MV, exclude=NULL, na.action=na.pass)
hist(bat.data.MV$mtime)

bat.data.MV[ bat.data.MV$mtime==0, c("mrise","mset","mtime")][1:10,]

# moon time 2 hrs around sunset or sunrise when bats most active
bat.data.MV$mtime2h<- as.numeric(bat.data.MV$mtime2h)
sum(is.na(bat.data.MV$mtime2h))
xtabs(~mtime2h,data=bat.data.MV, exclude=NULL, na.action=na.pass)
hist(bat.data.MV$mtime2h)

bat.data.MV[ bat.data.MV$mtime2h==0, c("mrise","mset","mtime")][1:10,]

bat.data.MV$nightlen <- as.numeric(bat.data.MV$nightlen)
sum(is.na(bat.data.MV$nightlen))
hist(bat.data.MV$nightlen)

#-----------------------




