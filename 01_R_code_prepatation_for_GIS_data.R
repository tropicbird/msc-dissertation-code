##########################
# Most of the code in this file is the modification of the code from
# "Best Practices for Using eBirdData. Version 1.0." (Strimas-Mackey et al., 2020)
# (https://cornelllabofornithology.github.io/ebird-best-practices/)
##########################


rm(list=ls(all=T))
gc(reset = TRUE)
library(sf)
library(rnaturalearth)
library(dplyr)

# file to save spatial data
PATH='C:/Users/DAISUKE/Dropbox/Leeds/Modules/dissertation/eBirdR/'

# select one directory
DIRECTORY='data_recwoo'
DIRECTORY='data_nutwoo'
DIRECTORY='data_lewoo'


gpkg_dir <- "/data"
if (!dir.exists(paste(PATH,DIRECTORY,gpkg_dir,sep=''))) {
  dir.create(paste(PATH,DIRECTORY,gpkg_dir,sep=''))
}

f_ne <- file.path(paste(PATH,DIRECTORY,gpkg_dir,sep=''), "gis-data.gpkg")

tmp_dir <- normalizePath(tempdir())
tmp_bcr <- file.path(tmp_dir, "bcr.zip")
paste0("https://www.birdscanada.org/download/gislab/", 
       "bcr_terrestrial_shape.zip") %>% 
  download.file(destfile = tmp_bcr)
unzip(tmp_bcr, exdir = tmp_dir)


bcr <- file.path(tmp_dir, "BCR_Terrestrial_master_International.shp") %>%
  read_sf() %>%
  select(bcr_code = BCR, bcr_name = LABEL) %>%
  filter(bcr_code == 9|bcr_code == 10)
# bcr_code == 32 for 'data_nutwoo'
# bcr_code == 27 for 'data_recwoo', 
# bcr_code == 9|bcr_code == 10 for 'data_lewwoo', 


# clean up
list.files(tmp_dir, "bcr", ignore.case = TRUE, full.names = TRUE) %>%
  unlink()

ne_land<-ne_download(scale=50,category='cultural',
                     type='admin_0_countries_lakes',
                     returnclass='sf')%>%
  filter(CONTINENT=="North America")%>%
  st_set_precision(1e6)%>%
  st_union()

# country lines
# downloaded globally then filtered to north america with st_intersect()
ne_country_lines<-ne_download(scale = 50, category = "cultural",
                              type = "admin_0_boundary_lines_land",
                              returnclass = "sf") %>%
  st_geometry()

ne_country_lines <- st_intersects(ne_country_lines, ne_land, sparse = FALSE) %>%
  as.logical() %>%
  {ne_country_lines[.]}

ne_state_lines <- ne_download(scale = 50, category = "cultural",
                              type = "admin_1_states_provinces_lines",
                              returnclass = "sf") %>%
  filter(adm0_a3 %in% c("USA", "CAN")) %>%
  mutate(iso_a2 = recode(adm0_a3, USA = "US", CAN = "CAN")) %>%
  select(country = adm0_name, country_code = iso_a2)

unlink(f_ne)
write_sf(ne_land, f_ne, "ne_land")
write_sf(ne_country_lines, f_ne, "ne_country_lines")
write_sf(ne_state_lines, f_ne, "ne_state_lines")
write_sf(bcr, f_ne, "bcr")