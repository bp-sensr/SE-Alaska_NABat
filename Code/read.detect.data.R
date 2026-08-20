# read the bat detection data base
library(tidyverse)

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


#bat.data <- readxl::read_excel("C:/Users/cami/Documents/BC-NABat/Data/All Year Activity by Night 2016-2024 BC V2 - TEMP.xlsx",
# bat.data <- readxl::read_excel(file.path("..","Data","All Year Activity by Night 2016-2023 Manual.xlsx"),
                             ###Toggle between these two lines to switch between all files and only files that passed auto ID
                             # sheet="OccupancyNoID",
                            # sheet="Occupancy",
                             #.name_repair='universal')

bat.data <- readxl::read_excel("V:/ARU/SENSR-BAT/NABat/2025/AK/Analyzed/All Year Activity by Night 2016-2025 AK-Aug19.xlsx",
                               sheet="OccupancyKSP",
                               .name_repair='universal') 

#time.data <- readxl::read_excel(file.path("..","Data","All Year Activity Time 2016-2023 AK and BC.xlsx"),
#                               sheet="Activity Time",
#                               .name_repair='universal')
#change Night to an actual recognisable date
bat.data <- bat.data %>%
  mutate(Night = as.Date(Night))

# Define columns that should stay as character/text
keep_as_text <- c("Quadrant", "Night", "Name", "SiteName", "Habitat.Type", 
                  "Feature.Sampled", "WaterType", "Region","srise", "sset",
                  "mrise",	"mset",	"moonr2",	"moons2"
                  )
# Convert everything else to numeric
for(col in names(bat.data)) {
  if(!col %in% keep_as_text) {
    bat.data[[col]] <- as.numeric(bat.data[[col]])
  }
}

bat.data$rnum <- 1:nrow(bat.data)
names(bat.data)
bat.data$Night <- as.Date(bat.data$Night)
bat.data$Year <- lubridate::year(bat.data$Night)

# check year of data
xtabs(~Region+Year, data=bat.data, exclude=NULL, na.action=na.pass)
length(unique(bat.data$GRTS.Cell.ID))

# # Email 2021-02-16 Ignore HIGHF and LOWF fields
bat.data$HIGHF <- NULL
bat.data$LOWF  <- NULL


#-----------------------
# identify the bat species and combination of batspecies column names in the bat.dat
# email 2021-02-16 Ignore the TABR species

species.id <- c("LACI", "LANO", "MYCA", 
                "MYEV", "MYLU", "MYVO", "MYYU")
bat.data.species <- names(bat.data)
bat.data.species <- bat.data.species[ (nchar(bat.data.species) %% 4 ) == 0]
bat.data.species <- bat.data.species[ bat.data.species == toupper(bat.data.species)]
bat.data.species

all.bat.data <- append(bat.data.species,c("NoID"),after = length(bat.data.species))

bat.data$Fullcount <- apply(bat.data[, all.bat.data],1,sum)

check.species <- bat.data.species[substr(bat.data.species,1,4) %in% species.id]
check.species

excel.species <- substr(bat.data.species,1,4)
setdiff(excel.species, species.id)
setdiff(species.id, excel.species)

bat.data.species <- bat.data.species[ substr(bat.data.species,1,4) %in% species.id]
bat.data.species

#-----------------------
# check the dates on the data
xtabs(~GRTS.Cell.ID+Year, data=bat.data, exclude=NULL, na.action=na.pass)



# check dates when data are collected by year
# email of 2021-02-16. Truncate at the end of August which is about Julian Day 240
bat.data$jNight <- lubridate::yday(bat.data$Night)
myplot <- ggplot(data=bat.data[ bat.data$Year > 2010,], aes(x=jNight, y=Year))+
  ggtitle("Julian nights of data collection")+
  geom_point(position=position_jitter(h=.1, w=.1))
myplot
ggsave(myplot, file=file.path("c:","Users","cami", "Documents","SE-Alaska_NABat","Figures","Summary", "jNight-vs-year.png"), h=4,w=6, units="in", dpi=300)

xtabs(~jNight, data=bat.data, exclude=NULL, na.action=na.pass)
#bat.data1 <- bat.data #keep a copy with all of the dates just in case. But create a new df that only has data before julian 240 (End of August)
bat.data <- bat.data[ bat.data$jNight < 240,]
bat.data$jNightM15 <- bat.data$jNight - lubridate::yday("2020-05-14")

#-----------------------
# check lat long
ggplot(data=bat.data, aes(x=Long, y=Lat,color=as.factor(Year)))+
  ggtitle("Where collected")+
  geom_point( position=position_jitter(h=.1, w=.1), shape=1)


#---------------
# check quadrant
xtabs(~Year+Quadrant, data=bat.data, exclude=NULL, na.action=na.pass)
bat.data$Quadrant <- paste0(toupper(substr(bat.data$Quadrant,1,1)), substring(bat.data$Quadrant,2))
xtabs(~Year+Quadrant, data=bat.data, exclude=NULL, na.action=na.pass)

xtabs(~GRTS.Cell.ID+Quadrant, data=bat.data, exclude=NULL, na.action=na.pass)
#-----------------------
# check habitat type
xtabs(~Habitat.Type+Year, data=bat.data, exclude=NULL, na.action=na.pass)

# extract the first (primary) habitat type
bat.data$Habitat2 <- trimws(substr(bat.data$Habitat.Type,1, -1+regexpr(",", paste0(bat.data$Habitat.Type,","), fixed=TRUE)))
xtabs(~Habitat2+Year, data=bat.data, exclude=NULL, na.action=na.pass)


#-----------------------
# email 2021-02-16 Drop from analysis
#xtabs(~Feature.Sampled+Year, data=bat.data, exclude=NULL, na.action=na.pass)

#----------------------
xtabs(~Distance.to.Clutter..m.,data=bat.data, exclude=NULL, na.action=na.pass)
bat.data$Distance.to.Clutter..m. <- as.numeric(bat.data$Distance.to.Clutter..m.)
bat.data$Distance.to.Clutter..m.[bat.data$Distance.to.Clutter..m. > 80] <- 80
xtabs(~Distance.to.Clutter..m.,data=bat.data, exclude=NULL, na.action=na.pass)

#-----------------------
xtabs(~Percent.Clutter,data=bat.data, exclude=NULL, na.action=na.pass)
bat.data$Percent.Clutter <- as.numeric(bat.data$Percent.Clutter)
xtabs(~Quadrant+Percent.Clutter,data=bat.data, exclude=NULL, na.action=na.pass)


#-----------------------
# confirm that 1=yes
# Email 2021-02-16. Yes, 1=Yes
xtabs(~Quadrant+Water.Nearby,data=bat.data, exclude=NULL, na.action=na.pass)

#-----------------------
# Convert second of daylength to decimal hours
# 
# bat.data$dayl <- bat.data$dayl / 3600
# bat.data$dayl.r <- round(bat.data$dayl)
# xtabs(~dayl,data=bat.data, exclude=NULL, na.action=na.pass)
# xtabs(~dayl.r,data=bat.data, exclude=NULL, na.action=na.pass)
# sum(is.na(bat.data$dayl))
# 
# #-----------------------
# # Seem unlikely that prcp is measured to 1/100 mm
# bat.data$prcp.r <- round(bat.data$prcp)
# xtabs(~prcp,data=bat.data, exclude=NULL, na.action=na.pass)
# xtabs(~prcp.r,data=bat.data, exclude=NULL, na.action=na.pass)

# added a column for nightly precipitation???NOt sure what negative values mean. Most of these are also missing so not useful
xtabs(~Nightly.Precipitation, data=bat.data, exclude=NULL, na.action=na.pass)
bat.data$Nightly.Precipitation <- round(bat.data$Nightly.Precipitation,.1)
xtabs(~Nightly.Precipitation+Year, data=bat.data, exclude=NULL, na.action=na.pass)

#create nightly.mean.temp?



# #-----------------------
# xtabs(~srad,data=bat.data, exclude=NULL, na.action=na.pass)
# sum(is.na(bat.data$srad))
# hist(bat.data$srad)
# 
# #-----------------------
# xtabs(~tmax,data=bat.data, exclude=NULL, na.action=na.pass)
# sum(is.na(bat.data$tmax))
# hist(bat.data$tmax)
# bat.data[ bat.data$tmax < 10 & !is.na(bat.data$tmax),c("Year","GRTS.Cell.ID","Lat","Long","Night","tmax","tmin")]

bat.data$Nightly.Max.Temp <- as.numeric(bat.data$Nightly.Max.Temp)
sum(is.na(bat.data$Nightly.Max.Temp))
xtabs(~Nightly.Max.Temp,data=bat.data, exclude=NULL, na.action=na.pass)
hist(bat.data$Nightly.Max.Temp)
bat.data[ bat.data$Nightly.Max.Temp < 10 & !is.na(bat.data$Nightly.Max.Temp),c("Year","GRTS.Cell.ID","Lat","Long","Night","Nightly.Max.Temp","Nightly.Min.Temp")]

# #-----------------------
# xtabs(~tmin,data=bat.data, exclude=NULL, na.action=na.pass)
# sum(is.na(bat.data$tmin))
# hist(bat.data$tmin)

bat.data$Nightly.Min.Temp <- as.numeric(bat.data$Nightly.Min.Temp)
sum(is.na(bat.data$Nightly.Min.Temp))
xtabs(~Nightly.Min.Temp,data=bat.data, exclude=NULL, na.action=na.pass)
hist(bat.data$Nightly.Min.Temp)


# #-----------------------
# xtabs(~dayl,data=bat.data, exclude=NULL, na.action=na.pass)
# sum(is.na(bat.data$dayl))
# hist(bat.data$dayl)

# #-----------------------
# xtabs(~vp,data=bat.data, exclude=NULL, na.action=na.pass)
# sum(is.na(bat.data$vp))
# hist(bat.data$vp)
# 
# #-----------------------
# xtabs(~srad,data=bat.data, exclude=NULL, na.action=na.pass)
# sum(is.na(bat.data$srad))
# hist(bat.data$srad)

#-----------------------
bat.data$Nightly.Min.RH <- as.numeric(bat.data$Nightly.Min.RH )
sum(is.na(bat.data$Nightly.Min.RH ))
xtabs(~Nightly.Min.RH ,data=bat.data, exclude=NULL, na.action=na.pass)
hist(bat.data$Nightly.Min.RH )

#-----------------------
bat.data$Nightly.Max.RH <- as.numeric(bat.data$Nightly.Max.RH)
sum(is.na(bat.data$Nightly.Max.RH))
xtabs(~Nightly.Max.RH,data=bat.data, exclude=NULL, na.action=na.pass)
hist(bat.data$Nightly.Max.RH)

#-----------------------
bat.data$Nightly.Mean.Windsp <- as.numeric(bat.data$Nightly.Mean.Windsp)
sum(is.na(bat.data$Nightly.Mean.Windsp))
xtabs(~Nightly.Mean.Windsp,data=bat.data, exclude=NULL, na.action=na.pass)
hist(bat.data$Nightly.Mean.Windsp)

#-----------------------
# Moon rise is read in a Excel decimal days
# See https://stackoverflow.com/questions/19172632/converting-excel-datetime-serial-number-to-r-datetime

#bat.data$mrise <-  as.POSIXct(as.numeric(bat.data$mrise) * (60*60*24), origin="1899-12-30", tz="UTC")
#bat.data$mset  <-  as.POSIXct(as.numeric(bat.data$mset ) * (60*60*24), origin="1899-12-30", tz="UTC")
#bat.data$srise <-  as.POSIXct(as.numeric(bat.data$srise) * (60*60*24), origin="1899-12-30", tz="UTC")
#bat.data$sset  <-  as.POSIXct(as.numeric(bat.data$sset ) * (60*60*24), origin="1899-12-30", tz="UTC")

bat.data$mtime <- as.numeric(bat.data$mtime)
sum(is.na(bat.data$mtime))
xtabs(~mtime,data=bat.data, exclude=NULL, na.action=na.pass)
hist(bat.data$mtime)

bat.data[ bat.data$mtime==0, c("mrise","mset","mtime")][1:10,]

# moon time 2 hrs around sunset or sunrise when bats most active
bat.data$mtime2h<- as.numeric(bat.data$mtime2h)
sum(is.na(bat.data$mtime2h))
xtabs(~mtime2h,data=bat.data, exclude=NULL, na.action=na.pass)
hist(bat.data$mtime2h)

bat.data[ bat.data$mtime2h==0, c("mrise","mset","mtime")][1:10,]

bat.data$nightlen <- as.numeric(bat.data$nightlen)
sum(is.na(bat.data$nightlen))
hist(bat.data$nightlen)

#-----------------------
# Look at the counts
head(bat.data[, c("Year","GRTS.Cell.ID","Night", bat.data.species)])
bat.data$Count <- apply(bat.data[, bat.data.species],1,sum, na.rm=TRUE)
xtabs(~Count, data=bat.data, exclude=NULL, na.action=na.pass)




