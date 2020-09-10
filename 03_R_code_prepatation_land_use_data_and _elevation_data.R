##########################
# Most of the code in this file is the modification of the code from
# "Best Practices for Using eBirdData. Version 1.0." (Strimas-Mackey et al., 2020)
# (https://cornelllabofornithology.github.io/ebird-best-practices/)
##########################

rm(list=ls(all=T))
getwd()
library(sf)
library(raster)
library(MODIS)
library(exactextractr)
library(viridis)
library(tidyverse)

# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map
projection <- raster::projection

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
file3=paste("ebd_",species,"_bcr",bcrNumber,"_zf.csv",sep="")


bcr <- read_sf(paste(PATH,DIRECTORY,"/data/gis-data.gpkg",sep=""),"bcr") %>% 
  #filter(bcr_code == as.integer(bcrNumber)) %>% #For nutwoo and recwoo
  filter(bcr_code == 9|bcr_code == 10) %>% #For lewwoo
  # project to the native modis projection
  st_transform(crs = paste("+proj=sinu +lon_0=0 +x_0=0 +y_0=0",
                           "+a=6371007.181 +b=6371007.181 +units=m +no_defs"))

ebird=read.csv(paste(PATH,DIRECTORY,'/',file3,sep=""))
tiles<-getTile(bcr)
tiles@tile

#MODIS setup
MODIS::EarthdataLogin(usr = "XXXX", pwd = "XXXX") #<-usr and pwd are hidden
MODIS:::checkTools("GDAL")
MODIS::MODISoptions(gdalPath = "c:/OSGeo4W64/bin")
MODIS:::checkTools("GDAL")

#Download using R
# earliest year of ebird data
begin_year <- format(as.Date(min(ebird$observation_date)), "%Y.01.01") #as.Date added!!

# end date for ebird data
end_year <- format(as.Date(max(ebird$observation_date)), "%Y.12.31") #as.Date added!!

# download tiles and combine into a single raster for each year
tifs <- runGdal(product = "MCD12Q1", collection = "006", SDSstring = "01", 
                extent = bcr %>% st_buffer(dist = 10000), 
                begin = begin_year, end = end_year, 
                outDirPath = paste(PATH,DIRECTORY,sep=""), job = "modis",
                MODISserverOrder = "LPDAAC") %>% 
  pluck("MCD12Q1.006") %>% 
  unlist()

# rename tifs to have more descriptive names
new_names <- format(as.Date(names(tifs)), "%Y") %>% 
  sprintf("modis_mcd12q1_umd_%s.tif", .) %>% 
  file.path(dirname(tifs), .)
file.rename(tifs, new_names)

# load the landcover data
landcover <- list.files(paste(PATH,DIRECTORY,"/modis",sep=""), "^modis_mcd12q1_umd", 
                        full.names = TRUE) %>% stack()

# label layers with year
landcover <- names(landcover) %>% 
  str_extract("(?<=modis_mcd12q1_umd_)[0-9]{4}") %>% 
  paste0("y", .) %>% 
  setNames(landcover, .)
landcover

#which year is the most recent?->2018
max_lc_year <- names(landcover) %>% 
  str_extract("[0-9]{4}") %>% 
  as.integer() %>% 
  max()

#Landscape metrics
neighborhood_radius <- 5 * ceiling(max(res(landcover))) / 2
neighborhood_radius 

ebird_buff <- ebird %>% 
  distinct(year = format(as.Date(observation_date), "%Y"), 
           locality_id, latitude=round(latitude,4), longitude=round(longitude,4)) %>% #as.Date() added!!!, #round Added!!
  mutate(year_lc = if_else(as.integer(year) > max_lc_year, 
                           as.character(max_lc_year), year),
         year_lc = paste0("y", year_lc)) %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  # transform to modis projection
  st_transform(crs = projection(landcover)) %>% 
  # buffer to create neighborhood around each point
  st_buffer(dist = neighborhood_radius) %>% 
  # nest by year
  nest(data = c(year, locality_id, geometry))

# function to summarize landcover data for all checklists in a given year
calculate_pland <- function(yr, regions, lc) {
  locs <- st_set_geometry(regions, NULL)
  exact_extract(lc[[yr]], regions, progress = FALSE) %>% 
    map(~ count(., landcover = value)) %>% 
    tibble(locs, data = .) %>% 
    unnest(data)
}

# iterate over all years extracting landcover for all checklists in each
lc_extract <- ebird_buff %>% 
  mutate(pland = map2(year_lc, data, calculate_pland, lc = landcover)) %>% 
    select(pland) %>% 
      unnest(cols = pland)

lc_extract

#PLAND: the proportion of the neighborhood within each land cover class.
pland <- lc_extract %>% 
  # calculate proporiton
  group_by(locality_id, year) %>% 
  mutate(pland = n / sum(n)) %>% 
  ungroup() %>% 
  select(-n) %>% 
  # remove NAs after tallying so pland is relative to total number of cells
  filter(!is.na(landcover))

lc_names <- tibble(landcover = 0:15,
                   lc_name = c("pland_00_water", 
                               "pland_01_evergreen_needleleaf", 
                               "pland_02_evergreen_broadleaf", 
                               "pland_03_deciduous_needleleaf", 
                               "pland_04_deciduous_broadleaf", 
                               "pland_05_mixed_forest",
                               "pland_06_closed_shrubland", 
                               "pland_07_open_shrubland", 
                               "pland_08_woody_savanna", 
                               "pland_09_savanna", 
                               "pland_10_grassland", 
                               "pland_11_wetland", 
                               "pland_12_cropland", 
                               "pland_13_urban", 
                               "pland_14_mosiac", 
                               "pland_15_barren"))
pland <- pland %>% 
  inner_join(lc_names, by = "landcover") %>% 
  arrange(landcover) %>% 
  select(-landcover)

# transform to wide format, filling in implicit missing values with 0s%>% 
pland <- pland %>%
  pivot_wider(names_from = lc_name,
              values_from = pland,
              values_fill = list(pland = 0))

write_csv(pland, paste(PATH,DIRECTORY,"/modis_pland_location-year.csv",sep=""))


#Prediction surface
agg_factor <- round(2 * neighborhood_radius / res(landcover))
r <- raster(landcover) %>% 
  aggregate(agg_factor) 
r <- bcr %>% 
  st_transform(crs = projection(r)) %>% 
  rasterize(r, field = 1) %>% 
  # remove any empty cells at edges
  trim()

r <- writeRaster(r, filename = paste(PATH,DIRECTORY,"/prediction-surface.tif",sep=""), overwrite = TRUE)


# get cell centers and create neighborhoods
r_centers <- rasterToPoints(r, spatial = TRUE) %>% 
  st_as_sf() %>% 
  transmute(id = row_number())
r_cells <- st_buffer(r_centers, dist = neighborhood_radius)

# extract landcover values within neighborhoods, only needed most recent year
lc_extract_pred <- landcover[[paste0("y", max_lc_year)]] %>% 
  exact_extract(r_cells, progress = FALSE) %>% 
  map(~ count(., landcover = value)) %>% 
  tibble(id = r_cells$id, data = .) %>% 
  unnest(data)

save(lc_extract_pred,file=paste(PATH,DIRECTORY,"/lc_extract_pred_tmp.RData",sep=""))

# calculate the percent for each landcover class
pland_pred <- lc_extract_pred %>% 
  count(id, landcover) %>% 
  group_by(id) %>% 
  mutate(pland = n / sum(n)) %>% 
  ungroup() %>% 
  select(-n) %>% 
  # remove NAs after tallying so pland is relative to total number of cells
  filter(!is.na(landcover))

# convert names to be more descriptive
pland_pred <- pland_pred %>% 
  inner_join(lc_names, by = "landcover") %>% 
  arrange(landcover) %>% 
  select(-landcover)

# tranform to wide format, filling in implicit missing values with 0s
pland_pred <- pland_pred %>% 
  pivot_wider(names_from = lc_name, 
              values_from = pland, 
              values_fill = list(pland = 0)) %>% 
  mutate(year = max_lc_year) %>% 
  select(id, year, everything())

# join in coordinates
pland_coords <- st_transform(r_centers, crs = 4326) %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  cbind(id = r_centers$id, .) %>% 
  rename(longitude = X, latitude = Y) %>% 
  inner_join(pland_pred, by = "id")


#Elevation
elev <- raster(paste(PATH,DIRECTORY,"/elevation_1KMmd_GMTEDmd.tif",sep=""))

# crop, buffer bcr by 10 km to provide a little wiggly room
elev <- bcr %>% 
  st_buffer(dist = 10000) %>% 
  st_transform(crs = projection(elev)) %>% 
  crop(elev, .) %>% 
  projectRaster(crs = projection(landcover))

# buffer each checklist location
ebird_buff_noyear <- ebird %>% 
  distinct(locality_id, latitude=round(latitude,4), longitude=round(longitude,4)) %>% #Added round!!!!
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(elev)) %>% 
  st_buffer(dist = neighborhood_radius)

# extract elevation values and calculate median and sd
locs <- st_set_geometry(ebird_buff_noyear, NULL) %>% 
  mutate(id = row_number())

#This takes time!!!
elev_checklists <- exact_extract(elev, ebird_buff_noyear, progress = FALSE) %>% 
  map_dfr(~ tibble(elevation_median = mean(.$value, na.rm = TRUE),
                   elevation_sd = sd(.$value, na.rm = TRUE))) %>% 
  # join to lookup table to get locality_id
  bind_cols(locs, .)

save(elev_checklists,file=paste(PATH,DIRECTORY,"/elev_checklist_tmp.RData",sep=""))

# extract and calculate median and sd
#This takes time!!!
elev_pred <- exact_extract(elev, r_cells, progress = FALSE) %>% 
  map_dfr(~ tibble(elevation_median = mean(.$value, na.rm = TRUE),
                   elevation_sd = sd(.$value, na.rm = TRUE))) %>% 
  # join to lookup table to get locality_id
  bind_cols(st_drop_geometry(r_cells), .)

save(elev_pred,file=paste(PATH,DIRECTORY,"/elev_pred_tmp.RData",sep=""))

# checklist covariates
pland=read_csv(paste(PATH,DIRECTORY,"/modis_pland_location-year.csv",sep=""))#<-If you start from elevation
pland_elev_checklist <- inner_join(pland, elev_checklists, by = "locality_id")

write_csv(pland_elev_checklist, paste(PATH,DIRECTORY,"/pland-elev_location-year.csv",sep=""))

# prediction surface covariates
pland_elev_pred <- inner_join(pland_coords, elev_pred, by = "id")
write_csv(pland_elev_pred, paste(PATH,DIRECTORY,"/pland-elev_prediction-surface.csv",sep=""))

glimpse(pland_elev_pred)