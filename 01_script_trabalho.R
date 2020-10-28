#' ---
#' title: Trabalho final - Análises Espaciais - Prof. Maurício Vancine
#' author: Gedimar Barbosa e Camila Vieira
#' date: 2020-10-27
#' ---

#Loading required packages

library(tidyverse)
library(igraph)
library(raster)
library(geobr)


# download dataset
dir.create("02_data")

download.file(url = "https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.1818&file=ecy1818-sup-0002-DataS1.zip",
              destfile = here::here("02_data", "atlantic_frugivory.zip"), mode = "wb")

# unzip
unzip(zipfile = here::here("02_data", "atlantic_frugivory.zip"),
      exdir = here::here("02_data"))

setwd("D:/Ge/github/disciplina-mauricio_trabalho-final/02_data")

dataset <- read.csv("ATLANTIC_frugivory.csv")

#Cleaning data set

##Arranging dataset by Latitude and Longitude (lat_lon) and droping NAS
da <- dataset %>% 
  tidyr::drop_na(Latitude:Longitude) %>% 
  tidyr::unite("lon_lat", Longitude:Latitude, sep = "_", remove = F) %>% 
  dplyr::arrange(lon_lat)

##Filtering interactions with precise location (lat_lon)
da <- da %>%
  dplyr::filter(Precision == "preciso")


#Creating a dataset for part 1 - networks of interactions

##da_int (interaction colum)
da_int <- da %>% 
  tidyr::unite("interaction", Frugivore_Species:Plant_Species, sep = "_", remove = F)


#Creating a dataset for part 2 #40 communities

##da_net is a tibble list where each community is defined by lat_lon
da_net <- da %>% 
  group_by(lon_lat) %>% 
  group_split(lon_lat) 


##creating a df of edges for the multilayer network by lat_lon
net_edges <- da %>% 
  dplyr::select(Frugivore_Species, Plant_Species, lon_lat, Longitude:Latitude) %>% 
  dplyr::rename(from = Frugivore_Species, to = Plant_Species)

net_edges

net_edges <- net_edges %>% 
  group_by(lon_lat) %>% 
  group_split(lon_lat) 


##creating a df of nodes for the networks by lat_lon
net_nodes <- da %>%
  dplyr::select(Frugivore_Species, Plant_Species, lon_lat, Longitude:Latitude) %>% 
  tidyr::pivot_longer(Frugivore_Species:Plant_Species,
                      names_to = "sp_type", 
                      values_to = "species") %>%
  dplyr::relocate(species, .before = lon_lat) %>% 
  dplyr::relocate(sp_type, .before = lon_lat)


net_nodes <- net_nodes %>%
  dplyr::group_by(lon_lat) %>% 
  dplyr::group_split(lon_lat) 

net_nodes <- lapply(net_nodes, dplyr::distinct)


#Ploting networks as an example

##creating an igraph object

net <- igraph::graph_from_data_frame(d=net_edges[[1]], vertices = net_nodes[[1]], directed = F)


##checking igraph attrbutes
E(net)
V(net)
V(net)$sp_type

##Ploting network as an example

net <- igraph::simplify(net, remove.multiple = F, remove.loops = T)

V(net)$color = V(net)$sp_type
V(net)$color = gsub("Frugivore_Species","gold",V(net)$color)
V(net)$color = gsub("Plant_Species","green",V(net)$color)

plot(net, vertex.label=NA, edge.width=1.5, layout = layout_with_kk)

#measuring network metrics


### Spatial analysis
# Atlantic forest raster data

# download Atlantic Forest raster
raster::getData(name = "SRTM", lon = -47, lat = -23,#o getData baixa já para long/lat escolhida 
                path = here::here("02_data", "raster"))

# import raster
ra <- raster::raster(here::here("02_data", "raster", "srtm_27_17.tif"))
ra

# Atlantic forest - filter
af <- geobr::read_biomes(year = 2019) %>% 
  sf::st_transform(crs = 4326) %>% 
  dplyr::filter(name_biome == "Mata Atlântica")

af

# plot
raster::plot(af$geom, col = "tomato")


# environmental data from worldclim

dir.create(here::here("02_data", "raster"))

# download
download.file(url = "https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_10m_elev.zip",
              destfile = here::here("02_data", "raster", "wc2.1_10m_elev.zip"), mode = "wb")

# unzip
unzip(zipfile = here::here("02_data", "raster", "wc2.1_10m_elev.zip"),
      exdir = here::here("02_data", "raster"))

# list files
fi <- dir(path = here::here("02_data", "raster"), pattern = "wc") %>% 
  grep(".tif", ., value = TRUE)
fi

# import raster
elev <- raster::raster(here::here("02_data", "raster", "wc2.1_10m_elev.tif"))


#raster::extract(elev, )

d_lon_lat <- da %>% 
  tidyr::drop_na(Longitude:Latitude) %>% 
  dplyr::select(Longitude:Latitude)
  
d_lon_lat
  
data_ele <- d_lon_lat %>% 
  dplyr::mutate(elev = raster::extract(elev, .))
data_ele

#placing elevation data in da_int

da_int <- da_int %>% 
  dplyr::mutate(elev = data_ele$elev)
da_int

da_int <- da_int %>% 
  tidyr::drop_na(elev)

#where to get different altitude data for atlantic forest

# http://srtm.csi.cgiar.org/srtmdata/
# https://www.earthenv.org/DEM