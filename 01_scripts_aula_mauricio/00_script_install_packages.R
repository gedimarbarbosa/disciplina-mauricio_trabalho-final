#' ---
#' title: install packages
#' author: mauricio vancine
#' date: 2020-10-19
#' ---

# github packages ----------------------------------------------------------
if(!require(devtools)) install.packages("devtools", dep = TRUE)

# data manipulation and visualization -------------------------------------
# manipulation and visualization
if(!require(tidyverse)) install.packages("tidyverse", dep = TRUE)
if(!require(lubridate)) install.packages("lubridate", dep = TRUE)

# directory
if(!require(here)) install.packages("here", dep = TRUE)

# xlsx
if(!require(openxlsx)) install.packages("openxlsx", dep = TRUE)
if(!require(readxl)) install.packages("readxl", dep = TRUE)
if(!require(writexl)) install.packages("writexl", dep = TRUE)

# geospatial data ---------------------------------------------------------
# vectors
if(!require(geobr)) install.packages("geobr")
if(!require(rnaturalearth)) install.packages("rnaturalearth", dep = TRUE)
if(!require(sf)) install.packages("sf", dep = TRUE)

# raster
if(!require(raster)) install.packages("raster", dep = TRUE)

# projection
if(!require(lwgeom)) install.packages("lwgeom", dep = TRUE)

# maps
if(!require(ggspatial)) install.packages("ggspatial", dep = TRUE)
if(!require(tmap)) install.packages("tmap", dep = TRUE)
if(!require(viridis)) install.packages("viridis", dep = TRUE)
if(!require(RColorBrewer)) install.packages("RColorBrewer", dep = TRUE)
if(!require(cptcity)) install.packages("cptcity", dep = TRUE)
if(!require(wesanderson)) install.packages("wesanderson")

# end ---------------------------------------------------------------------