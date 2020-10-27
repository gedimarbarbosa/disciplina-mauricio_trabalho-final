#' ---
#' title: Trabalho final - Análises Espaciais - Prof. Maurício Vancine
#' author: Gedimar Barbosa e Camila Vieira
#' date: 2020-10-27
#' ---

#Loading required packages

library(tidyverse)
library(igraph)


# download dataset
dir.create("02_data")

download.file(url = "https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.1818&file=ecy1818-sup-0002-DataS1.zip",
              destfile = here::here("02_data", "atlantic_frugivory.zip"), mode = "wb")

# unzip
unzip(zipfile = here::here("02_data", "atlantic_frugivory.zip"),
      exdir = here::here("02_data"))

setwd("D:/Ge/github/disciplina-mauricio_trabalho-final/02_data")

dataset <- read.csv("ATLANTIC_frugivory.csv")

#Cleaning data set for our purpose

#Arranging dataset by Latitude and Longitude (lat_lon) and droping NAS
da <- dataset %>% 
  tidyr::drop_na(Latitude:Longitude) %>% 
  tidyr::unite("lat_lon", Latitude:Longitude, sep = ",") %>% 
  dplyr::arrange(lat_lon)

#Filtering interactions with precise location (lat_lon)
da <- da %>%
  dplyr::filter(Precision == "preciso")


#Creating a dataset for part 1 - networks of interactions

#da_int (interaction colum)
da_int <- da %>% 
  tidyr::unite("interaction", Frugivore_Species:Plant_Species, sep = "_", remove = F)


#Creating a dataset for part 2 #49 communities

#da_net is a tibble list where each community is defined by lat_lon
da_net <- da %>% 
  group_by(lat_lon) %>% 
  group_split(lat_long) 


#creating a df of edges for the multilayer network by lat_lon
net_edges <- da %>% 
  dplyr::select(Frugivore_Species, Plant_Species, lat_lon) %>% 
  dplyr::rename(from = Frugivore_Species, to = Plant_Species)

net_edges

net_edges <- net_edges %>% 
  group_by(lat_lon) %>% 
  group_split(lat_long) 


#creating a df of nodes for the networks by lat_lon
net_nodes <- da %>%
  dplyr::select(Frugivore_Species, Plant_Species, lat_lon) %>% 
  tidyr::pivot_longer(Frugivore_Species:Plant_Species,
                      names_to = "sp_type", 
                      values_to = "species") %>%
  dplyr::relocate(species, .before = lat_lon) %>% 
  dplyr::relocate(sp_type, .before = lat_lon)


net_nodes <- net_nodes %>%
  group_by(lat_lon) %>% 
  group_split(lat_long) 

net_nodes <- lapply(net_nodes, dplyr::distinct)


#Ploting networks

#creating an igraph object (example)

net <- igraph::graph_from_data_frame(d=net_edges[[17]], vertices = net_nodes[[17]], directed = F)


#checking igraph attrbutes for net_1
E(net)
V(net)
V(net)$sp_type

#Ploting network as an example

net <- igraph::simplify(net, remove.multiple = F, remove.loops = T)

V(net)$color = V(net)$sp_type
V(net)$color = gsub("Frugivore_Species","gold",V(net)$color)
V(net)$color = gsub("Plant_Species","green",V(net)$color)

plot(net, vertex.label=NA, edge.width=1.5, layout = layout_with_kk)


