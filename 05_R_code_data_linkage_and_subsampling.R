##########################
# Most of the code in this file is the modification of the code from
# "Best Practices for Using eBirdData. Version 1.0." (Strimas-Mackey et al., 2020)
# (https://cornelllabofornithology.github.io/ebird-best-practices/)
##########################

rm(list=ls(all=T))
gc(reset = TRUE)
library(sf)
library(raster)
library(dggridR)
library(lubridate)
library(ranger)
library(scam)
library(PresenceAbsence)
library(verification)
library(ebirdst)
library(fields)
library(gridExtra)
library(tidyverse)

# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map
projection <- raster::projection

set.seed(1)
getwd()
if(!dir.exists("output")){
  dir.create("output")
}

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

PATH2='D:/eBird_data/'
DIRECTORY=paste('data_',species,sep="")
file4=paste(species,"_with_prep_tmp_clm_add30yMonth.csv",sep="")

ebird<-read_csv(paste(PATH2,file4,sep='')) %>%
  mutate(year=year(observation_date))

habitat<-read_csv(paste(PATH,DIRECTORY,'/pland-elev_location-year.csv',sep='')) %>%
  mutate(year=as.integer(year))

# combine ebird and habitat data
ebird_habitat <- inner_join(ebird, habitat, by = c("locality_id", "year"))

# prediction surface
pred_surface <- read_csv(paste(PATH,DIRECTORY,"/pland-elev_prediction-surface.csv",sep=""))
# latest year of landcover data
max_lc_year <- max(pred_surface$year)
r <- raster(paste(PATH,DIRECTORY,"/prediction-surface.tif",sep=""))

# load gis data for making maps
map_proj <- st_crs(3857)#-----------102003->3857
ne_land <- read_sf(paste(PATH,DIRECTORY,"/gis-data.gpkg",sep=""), "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
bcr <- read_sf(paste(PATH,DIRECTORY,"/gis-data.gpkg",sep=""), "bcr") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_country_lines <- read_sf(paste(PATH,DIRECTORY,"/gis-data.gpkg",sep=""), "ne_country_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_state_lines <- read_sf(paste(PATH,DIRECTORY,"/gis-data.gpkg",sep=""), "ne_state_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()

#Spatiotemporal subsampling
# bounding box to generate points from
bb <- st_bbox(c(xmin = -0.1, xmax = 0.1, ymin = -0.1, ymax = 0.1), 
              crs = 4326) %>% 
  st_as_sfc() %>% 
  st_sf()
# random points
pts <- st_sample(bb, 500) %>% 
  st_sf(as.data.frame(st_coordinates(.)), geometry = .) %>% 
  rename(lat = Y, lon = X)

# contruct a hexagonal grid with ~ 5 km between cells
dggs <- dgconstruct(spacing = 5)
# for each point, get the grid cell
pts$cell <- dgGEO_to_SEQNUM(dggs, pts$lon, pts$lat)$seqnum

# sample one checklist per grid cell
pts_ss <- pts %>% 
  group_by(cell) %>% 
  sample_n(size = 1) %>% 
  ungroup()

# generate polygons for the grid cells
hexagons <- dgcellstogrid(dggs, unique(pts$cell), frame = FALSE) %>% 
  st_as_sf()
ggplot() +
  geom_sf(data = hexagons) +
  geom_sf(data = pts, size = 0.5) +
  geom_sf(data = pts_ss, col = "red") +
  theme_bw()

# generate hexagonal grid with ~ 5 km betweeen cells
dggs <- dgconstruct(spacing = 5)
dggs
# get hexagonal cell id and week number for each checklist
checklist_cell <- ebird_habitat %>% 
  mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum,
         year = year(observation_date),
         week = week(observation_date))

# sample one checklist per grid cell per week
# sample detection/non-detection independently 
ebird_ss <- checklist_cell %>% 
  group_by(species_observed, year, week, cell) %>% 
  sample_n(size = 1) %>% 
  ungroup()

# original data
nrow(ebird_habitat)
count(ebird_habitat, species_observed) %>% 
  mutate(percent = n / sum(n))

# after sampling
nrow(ebird_ss)
count(ebird_ss, species_observed) %>% 
  mutate(percent = n / sum(n))

# convert checklists to spatial features
all_pts <- ebird_habitat %>%  
  st_as_sf(coords = c("longitude","latitude"), crs = 4326) %>%
  st_transform(crs = map_proj) %>% 
  select(species_observed)
ss_pts <- ebird_ss %>%  
  st_as_sf(coords = c("longitude","latitude"), crs = 4326) %>%
  st_transform(crs = map_proj) %>% 
  select(species_observed)
both_pts <- list(before_ss = all_pts, after_ss = ss_pts)

write_csv(ebird_habitat, paste(PATH2,"ebird_habitat_",species,"_add30yMonth.csv",sep=''))
write_csv(ebird_ss, paste(PATH2,"ebird_ss_",species,"_add30yMonth.csv",sep=''))