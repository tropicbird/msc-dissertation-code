##########################
# Most of the code in this file is the modification of the code from
# "Best Practices for Using eBirdData. Version 1.0." (Strimas-Mackey et al., 2020)
# (https://cornelllabofornithology.github.io/ebird-best-practices/)
##########################

rm(list=ls(all=T))
gc(reset = TRUE)
library(auk)
library(lubridate)
library(sf)
library(gridExtra)
library(tidyverse)
library(rnaturalearth)
getwd()
#auk::auk_set_ebd_path("D:/eBird_data",overwrite = TRUE)
#?auk_set_ebd_path

PATH='C:/Users/DAISUKE/Dropbox/Leeds/Modules/dissertation/eBirdR/'

# select one species
species='nutwoo'
species='recwoo'
species='lewwoo'

# select one species
common_name="Nuttall's Woodpecker"
common_name='Red-cockaded Woodpecker'
common_name="Lewis's Woodpecker"

# select one bcrNumber
bcrNumber='32'
bcrNumber='27'
bcrNumber='9n10'

DIRECTORY=paste('data_',species,sep="")
eBirdData=paste('ebd_US_',species,'_relMay-2020.txt',sep="")

file1=paste("ebd_",species,"_bcr",bcrNumber,".txt",sep="")
file2=paste("ebd_checklists_bcr",bcrNumber,".txt",sep="")
file3=paste("ebd_",species,"_bcr",bcrNumber,"_zf.csv",sep="")

select<-dplyr::select #dplyr is a grammar of data manipulation

ebd<-auk_ebd(eBirdData,
             file_sampling="ebd_sampling_relMay-2020.txt")

ebd_filters<-ebd %>%
  auk_species(common_name)%>% 
  auk_bcr(bcr=c(32))%>% #bcr=27 #bcr=c(9,10)
  auk_date(date=c("2010-01-01","2019-12-31"))%>%
  auk_protocol(protocol=c("Stationary","Traveling"))%>%
  auk_complete()

f_ebd<-file.path(paste(PATH,DIRECTORY,sep=""),file1)
f_sampling<-file.path(paste(PATH,DIRECTORY,sep=""),file2)

if(!file.exists(f_ebd)){
  auk_filter(ebd_filters,file=f_ebd,file_sampling = f_sampling)
}

ebd_zf<-auk_zerofill(f_ebd, f_sampling, collapse=TRUE)
ebd_zf
names(ebd_zf)

time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}

# clean up variables
ebd_zf <- ebd_zf %>% 
  mutate(
    # convert X to NA
    observation_count = if_else(observation_count == "X", 
                                NA_character_, observation_count),
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for non-travelling counts
    effort_distance_km = if_else(protocol_type != "Traveling", 
                                 0, effort_distance_km),
    # convert time to decimal hours since midnight
    time_observations_started = time_to_decimal(time_observations_started),
    # split date into year and day of year
    year = year(observation_date),
    day_of_year = yday(observation_date)
  )

#Accounting for variation in detectability
# additional filtering
ebd_zf_filtered <- ebd_zf %>% 
  filter(
    # effort filters
    duration_minutes <= 5 * 60,
    effort_distance_km <= 5,
    # last XX years of data
    year >= 2010,
    year <= 2019,
    # 10 or fewer observers
    number_observers <= 10)


ebird <- ebd_zf_filtered %>% 
  select(checklist_id, observer_id, sampling_event_identifier,
         scientific_name,
         observation_count, species_observed, 
         state_code, locality_id, latitude, longitude,
         protocol_type, all_species_reported,
         observation_date, year, day_of_year,
         time_observations_started, 
         duration_minutes, effort_distance_km,
         number_observers)

write_csv(ebird, paste(PATH,DIRECTORY,'/',file3,sep=""), na = "")

ebird=read.csv(paste(PATH,DIRECTORY,'/',file3,sep=""))

#----Exploratory analysis and visualization
#Coordinate Reference System
#map_proj <- st_crs(102003)#<-----102003 doesn't work!!!!
map_proj <- st_crs(3857)#<----alternative

paste(PATH,DIRECTORY,'/data/gis-data.gpkg',sep="")

ne_land <- read_sf(paste(PATH,DIRECTORY,'/data/gis-data.gpkg',sep=""), "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()

bcr <- read_sf(paste(PATH,DIRECTORY,'/data/gis-data.gpkg',sep=""), "bcr") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()

ne_country_lines <- read_sf(paste(PATH,DIRECTORY,'/data/gis-data.gpkg',sep=""), "ne_country_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()

ne_state_lines <- read_sf(paste(PATH,DIRECTORY,'/data/gis-data.gpkg',sep=""), "ne_state_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
     
# prepare ebird data for mapping
ebird_sf <- ebird %>% 
  # convert to spatial points
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = map_proj) %>% 
  select(species_observed) #TRUE or FALSE

# map
par(mar = c(0.25, 0.25, 0.25, 0.25))
# set up plot area
plot(st_geometry(ebird_sf), col = NA)

# contextual gis data
plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)
plot(bcr, col = "#cccccc", border = NA, add = TRUE)
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)

# ebird observations
# not observed
plot(st_geometry(ebird_sf), #<---there are no filter, plot every thing??
     pch = 19, cex = 0.1, col = alpha("#555555", 0.25),
     add = TRUE)

# observed
plot(filter(ebird_sf, species_observed) %>% st_geometry(),
     pch = 19, cex = 0.1, col = alpha("#4daf4a", 0.25),
     add = TRUE)

legend("bottomright", bty = "n",
       col = c("#555555", "#4daf4a"),
       legend = c("Absence checklists", "Presence checklists"),
       pch = 19)

box()
par(new = TRUE, mar = c(0, 0, 3, 0))
title(paste(common_name," eBird Observations\n 2010-2019, BCR ","32",sep=""))