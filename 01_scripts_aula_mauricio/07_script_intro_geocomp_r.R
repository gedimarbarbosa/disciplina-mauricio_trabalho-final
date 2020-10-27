#' ---
#' title: aula 07 - estrutura e manipulação de dados matriciais
#' author: mauricio vancine
#' date: 2020-10-22
#' ---

# topics ------------------------------------------------------------------
# 7.1 pacotes
# 7.2 sados raster
# 7.3 classes raster
# 7.4 importar dados matriciais
# 7.5 descricao de objetos raster
# 7.6 converter crs
# 7.7 manipulação de dados raster
# 7.8 opera_rccao espaciais
# 7.9 opera_rccao geometricas
# 7.10 intera_rccoes raster-vetor
# 7.11 conversoes raster-vetor
# 7.12 exportar dados matriciais

# packages ----------------------------------------------------------------
library(raster)
library(sf)
library(tidyverse)
library(geobr)
library(rnaturalearth)
library(viridis)

# 7.1 pacotes -------------------------------------------------------------
# raster
# install.packages("raster")
# library(raster)

# terra_rc
# install.packages("terra_rc")
# library(terra_rc)

# stars
# install.packages("stars")
# library(stars)

# 7.3 classes raster ------------------------------------------------------
# volcano
volcano

# rasterlayer
ra_lay <- raster::raster(volcano)
ra_lay

# plot
raster::plot(ra_lay)

# plot
raster::plot(ra_lay, col = viridis::viridis(10))
raster::plot(ra_lay, col = viridis::viridis(100))

# stack (multicamadas de rasters/ bandas de arquivos fonte diferentes)
set.seed(42)
ra_sta <- raster::stack(raster::raster(volcano), 
                        raster::raster(matrix(rnorm(prod(dim(volcano))), nrow = 87)),
                        raster::raster(matrix(rpois(prod(dim(volcano)), 10), nrow = 87)),
                        raster::raster(matrix(rbinom(prod(dim(volcano)), 1, .5), nrow = 87)))
ra_sta

# plot
raster::plot(ra_sta, col = viridis::viridis(10))

# brick (parecido com stack, mas utilizado para importar unico arquivo com muitas bandas, multicamadas)
#quando se tem muitas camadas, o interesse é transformar de stack pra bricks, antes de processar
set.seed(42)
ra_bri <- raster::brick(raster::raster(volcano), 
                        raster::raster(matrix(rnorm(prod(dim(volcano))), nrow = 87)),
                        raster::raster(matrix(rpois(prod(dim(volcano)), 10), nrow = 87)),
                        raster::raster(matrix(rbinom(prod(dim(volcano)), 1, .5), nrow = 87)))
ra_bri

# plot
raster::plot(ra_bri, col = viridis::viridis(10))

# 7.4 importar dados matriciais -------------------------------------------
# create directory
dir.create(here::here("03_dados", "raster"))

# elevation data
# increase time to download
options(timeout = 600)

# download
raster::getData(name = "SRTM", lon = -47, lat = -23,#o getData baixa já para long/lat escolhida
                path = here::here("03_dados", "raster"))

# import raster
ra <- raster::raster(here::here("03_dados", "raster", "srtm_27_17.tif"))
ra

# rio claro
rc_2019 <- geobr::read_municipality(code_muni = 3543907, year = 2019) %>% 
  sf::st_transform(crs = 4326)
rc_2019

# plot
raster::plot(ra, col = viridis::viridis(10))
plot(rc_2019$geom, col = adjustcolor("red", .5), add = TRUE)

# mask - cortando so para o municipio de Rio Claro
ra_rc <- raster::crop(ra, rc_2019)
ra_rc

# plot
raster::plot(ra_rc, col = viridis::viridis(10))
plot(rc_2019$geom, col = adjustcolor("red", .3), lwd = 2, add = TRUE)

# climate data
# download
download.file(url = "https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_10m_bio.zip",
              destfile = here::here("03_dados", "raster", "wc2.0_10m_bio.zip"), mode = "wb")

# unzip
unzip(zipfile = here::here("03_dados", "raster", "wc2.0_10m_bio.zip"),
      exdir = here::here("03_dados", "raster"))

# list files
fi <- dir(path = here::here("03_dados", "raster"), pattern = "wc") %>% 
  grep(".tif", ., value = TRUE)
fi

# import stack
st <- raster::stack(here::here("03_dados", "raster", fi))
st

# map
raster::plot(st[[1:2]], col = viridis::viridis(10))

# 7.5 descricao de objetos raster -----------------------------------------
# raster
ra_rc

# class
class(ra_rc)

# dimension
dim(ra_rc)

# number of layers
nlayers(ra_rc)

# number of rows
nrow(ra_rc)

# number of columns
ncol(ra_rc)

# number of cells
ncell(ra_rc)

# raster resolution
res(ra_rc)

# stack resolution
res(st)

# raster extent
extent(ra_rc)

# stack extent
extent(st)

# projection
projection(ra_rc)

# projection
projection(st)

# raster names
names(ra_rc)

# stack names
names(st)

# raster values
getValues(ra_rc) %>% head
values(ra_rc) %>% head
ra_rc[] %>% head

# raster values histogram
ra_rc %>% 
  raster::values() %>% 
  hist(col = "steelblue", border = "white", main = NA, xlab = "Elevação (m)", ylab = "Frequência")

# stack values
values(st[[1:3]]) %>% head

# stack values - pairs
st[[1:3]] %>% 
  raster::values() %>% 
  tibble::as_tibble() %>% 
  dplyr::sample_n(1e3) %>% 
  pairs(cex = 1.4, pch = 20, col = adjustcolor("steelblue", .7))

# 7.6 conversao do crs ----------------------------------------------------
# projection
ra_rc

# proj4string utm 23 s
utm23 <- "+proj=utm +zone=23 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
utm23

# reprojection raster
ra_rc_utm23 <- raster::projectRaster(ra_rc, crs = utm23)
ra_rc_utm23

# reprojection raster
ra_rc_utm23 <- raster::projectRaster(ra_rc, crs = utm23, res = 90, method = "bilinear")
ra_rc_utm23

# reprojection vector
rc_2019_utm23 <- sf::st_transform(rc_2019, crs = utm23)
rc_2019_utm23

# plot
raster::plot(ra_rc_utm23, col = viridis::viridis(10))
plot(rc_2019_utm23$geom, col = adjustcolor("red", .3), lwd = 2, add = TRUE)

# crs global
# WGS84/GCS
st$wc2.1_10m_bio_1

# plot
raster::plot(st$wc2.1_10m_bio_1, col = viridis::viridis(10))

# proj4string mollweide
moll <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
moll

# reprojection
bio01_moll <- raster::projectRaster(st$wc2.1_10m_bio_1, crs = moll, res = 25e3, method = "bilinear")
bio01_moll

# plot
raster::plot(bio01_moll, col = viridis::viridis(10))

# proj4string winkel tripel
wintri <- "+proj=wintri +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
wintri

# reprojection
bio01_wintri <- raster::projectRaster(st$wc2.1_10m_bio_1, crs = wintri, res = 25e3, method = "bilinear")
bio01_wintri

# plot
raster::plot(bio01_wintri, col = viridis::viridis(10))

# proj4string eckert iv
eck4 <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
eck4

# reprojection
bio01_eck4 <- raster::projectRaster(st$wc2.1_10m_bio_1, crs = eck4, res = 25e3, method = "bilinear")
bio01_eck4

# plot
raster::plot(bio01_eck4, col = viridis::viridis(10))

# proj4string lambert 
laea <- "+proj=laea +x_0=0 +y_0=0 +lon_0=0 +lat_0=0"
laea

# reprojection
bio01_laea <- raster::projectRaster(st$wc2.1_10m_bio_1, crs = laea, res = 25e3, method = "bilinear")
bio01_laea

# plot
raster::plot(bio01_laea, col = viridis::viridis(10))

# 7.7 manipulacao de dados raster -----------------------------------------
# indexacao de linha-coluna ou id da celula
# raster - row 1, column 1
ra[1, 1]

# cell is 1
ra[1]

# selection layer in stack #é possível fazer um subset com a camada escolhida
st_bio01 <- raster::subset(st, "wc2.1_10m_bio_1")
st_bio01

# selection layer in stack
st_bio01 <- raster::raster(st, layer = 1)
st_bio01

# selection layer in stack
st_bio01 <- st[["wc2.1_10m_bio_1"]]
st_bio01

# selection layer in stack
st_bio01 <- st[[1]]
st_bio01

# selection layer in stack
st_bio01 <- st$wc2.1_10m_bio_1
st_bio01

# map
raster::plot(st_bio01, col = viridis::viridis(10))

# rename
# names
names(ra_rc)

# rename
names(ra_rc) <- "elevation"

# names
names(ra_rc)

# name
names(st)

# rename
names(st) <- c("bio01", paste0("bio", 10:19), paste0("bio", 2:9))

# names
names(st)

# summarize information
# mean from all cells
raster::cellStats(ra_rc, mean)

# mean from all cells
raster::cellStats(st[[1:3]], mean)

# count cells
raster::freq(ra_rc)

# count cells
raster::freq(st[[1:3]])

# 7.8 operacoes espaciais -------------------------------------------------
# sum
ra_rc2 <- ra_rc + ra_rc
ra_rc2

# map
raster::plot(ra_rc2, col = viridis::viridis(10))
plot(rc_2019$geom, col = adjustcolor("red", .2), add = TRUE)

# division
ra_rc_km <- ra_rc / 1e3
ra_rc_km

# map
raster::plot(ra_rc_km, col = viridis::viridis(10))
plot(rc_2019$geom, col = adjustcolor("red", .2), add = TRUE)

# log10
ra_rc_log10 <- log10(ra_rc)
ra_rc_log10

# map
raster::plot(ra_rc_log10, col = viridis::viridis(10))
plot(rc_2019$geom, col = adjustcolor("red", .2), add = TRUE)

# upper 600
ra_rc_up_600 <- ra_rc > 600
ra_rc_up_600

# map
raster::plot(ra_rc_up_600)
plot(rc_2019$geom, col = adjustcolor("red", .2), add = TRUE)

# bottom
ra_rc_bot_600 <- ra_rc <= 600
ra_rc_bot_600

# map
raster::plot(ra_rc_bot_600)
plot(rc_2019$geom, col = adjustcolor("red", .2), add = TRUE)

# values
ra_rc[] %>% head

# copy
ra_rc_up_700 <- ra_rc
ra_rc_up_700

# values bottom 700 == NA
ra_rc_up_700[ra_rc_up_700 < 700] <- NA
ra_rc_up_700

# map
raster::plot(ra_rc_up_700, col = viridis::viridis(10))
plot(rc_2019$geom, col = adjustcolor("red", .2), add = TRUE)

# reclassify
# matrix for reclassification
rcl  <- matrix(c(400,500,4, 
                 500,600,5, 
                 600,700,6,
                 700,800,7,
                 800,900,8,
                 900,1000,9), 
               ncol = 3, byrow = TRUE)
rcl

# reclassify - reclassificando os valores (é possível tornar os pixels semelhantes ou diferentes aqui)
ra_rc_rcl <- raster::reclassify(ra_rc, rcl = rcl)
ra_rc_rcl

# plot
raster::plot(ra_rc_rcl, col = viridis::viridis(10))
plot(rc_2019$geom, col = adjustcolor("red", .2), add = TRUE)

# moving window
ra_rc_mw_sd <- raster::focal(ra_rc, w = matrix(1, nrow = 3, ncol = 3), fun = sd)
ra_rc_mw_sd

# plot
raster::plot(ra_rc_mw_sd, col = viridis::viridis(10))
plot(rc_2019$geom, col = adjustcolor("red", .1), add = TRUE)

# zonal statistics
raster::zonal(ra_rc, ra_rc_rcl, fun = "mean")

# zonal statistics
ra_rc_stats <- data.frame(raster::zonal(ra_rc, ra_rc_rcl, fun = "summary"))
colnames(ra_rc_stats) <- c("zone", "min", "1qt", "median", "mean", "3qt", "max")
ra_rc_stats

# 7.9 operacoes geometricas -----------------------------------------------
# aggregation
# resolution
res(ra_rc)[1]

# aggregation - increases the pixel size of the raster
ra_rc_agg_mean <- raster::aggregate(ra_rc, fact = 10, fun = "mean")
ra_rc_agg_mean

# plot
raster::plot(ra_rc_agg_mean, col = viridis::viridis(10))
plot(rc_2019$geom, col = adjustcolor("red", .2), add = TRUE)

# disaggregate
# resolution
res(ra_rc)[1]

# disaggregate - decreases the pixel size of the raster
ra_rc_dis_bil <- raster::disaggregate(ra_rc, fact = 2, fun = "bilinear")
ra_rc_dis_bil

# plot
raster::plot(ra_rc_dis_bil, col = viridis::viridis(10))
plot(rc_2019$geom, col = adjustcolor("red", .2), add = TRUE)

# 7.10 interacoes raster vetor --------------------------------------------
# raster cropping
# crop - adjust extension
ra_rc_crop <- raster::crop(ra, rc_2019)
ra_rc_crop

# plot
raster::plot(ra_rc_crop, col = viridis::viridis(10))
plot(rc_2019$geom, col = adjustcolor("red", .2), add = TRUE)

# mask - adjust limit
ra_rc_mask <- raster::mask(ra, rc_2019)
ra_rc_mask

# plot
raster::plot(ra_rc_mask, col = viridis::viridis(10))
plot(rc_2019$geom, col = adjustcolor("red", .2), add = TRUE)

# crop and mask - adjust extension and limit
ra_rc_crop_mask <- ra %>% 
  raster::crop(rc_2019) %>% 
  raster::mask(rc_2019)
ra_rc_crop_mask

# plot
raster::plot(ra_rc_crop_mask, col = viridis::viridis(10))
plot(rc_2019$geom, col = adjustcolor("red", .2), add = TRUE)

# crop and mask inverse - adjust extension and limit
ra_rc_crop_mask_inv <- ra %>% 
  raster::crop(rc_2019) %>% 
  raster::mask(rc_2019, inverse = TRUE)
ra_rc_crop_mask_inv

# plot
raster::plot(ra_rc_crop_mask_inv, col = viridis::viridis(10))
plot(rc_2019$geom, col = adjustcolor("red", .2), add = TRUE)

# download atlantic forest
af <- geobr::read_biomes(year = 2004) %>% 
  dplyr::filter(name_biome == "Mata Atlântica")
af

# plot
plot(af$geom, col = "gray", main = NA, axes = TRUE, gratidule = TRUE)

# crop and mask - adjust extension and limit
st_af_crop_mask <- st %>% 
  raster::crop(af) %>% 
  raster::mask(af)
st_af_crop_mask

# plot
raster::plot(st_af_crop_mask[[1:4]], col = viridis::viridis(10))

# import points
rc_spr <- sf::st_read(here::here("03_dados", "vetor", "SP_3543907_NASCENTES.shp"), quiet = TRUE) %>% 
  sf::st_transform(crs = 4326)
rc_spr

# plot
plot(rc_spr$geometry, pch = 20, col = "blue", main = NA, axes = TRUE, graticule = TRUE)
plot(rc_2019$geom, add = TRUE)

# extraction
rc_spr_ele <- raster::extract(ra_rc, rc_spr)
rc_spr_ele

# extraction
rc_spr_ele <- rc_spr %>% 
  dplyr::mutate(elev = raster::extract(ra_rc, .))
rc_spr_ele

# plot
plot(rc_spr_ele["elev"], pch = 20, main = NA, axes = TRUE, graticule = TRUE)

# histogram
rc_spr_ele %>% 
  dplyr::pull(elev) %>% 
  hist(col = "steelblue", border = "white", main = NA, xlab = "Elevação (m)", ylab = "Frequência")

# zonal statistics
# buffers
set.seed(42)
rc_spr_buf <- rc_spr %>% 
  dplyr::sample_n(10) %>% 
  sf::as_Spatial() %>% 
  raster::buffer(width = 1000, dissolve = FALSE) %>% 
  sf::st_as_sf()
rc_spr_buf

# plot
plot(rc_spr_buf$geometry, col = adjustcolor("steelblue", .7), pch = 20, main = NA, axes = TRUE, graticule = TRUE)
plot(rc_2019$geom, add = TRUE)

# zonal statistics
raster::extract(x = ra_rc, y = rc_spr_buf, fun = mean, na.rm = TRUE, df = TRUE)

# zonal statistics
rc_spr_buf <- rc_spr_buf %>% 
  dplyr::mutate(elev_mean = raster::extract(x = ra_rc, y = rc_spr_buf, fun = mean, na.rm = TRUE))
rc_spr_buf

# plot
plot(rc_spr_buf["elev_mean"], pch = 20, main = NA, axes = TRUE, graticule = TRUE)

# 7.11 Conversoes raster-vetor --------------------------------------------
# rasterize points
rc_spr_rast <- raster::rasterize(x = rc_spr, y = ra_rc_agg_mean, field = 1, fun = "count")
rc_spr_rast

# plot
raster::plot(rc_spr_rast, col = viridis::viridis(10))
plot(rc_spr$geometry, pch = 20, cex = 1, col = adjustcolor("gray50", .6), add = TRUE)

# statistics
plot(ra_rc_agg_mean[], rc_spr_rast[], col = adjustcolor("gray30", .5), 
     pch = 20, cex = 1.5, bty = "l", xlab = "Elevação (m)", ylab = "Número de nascentes")

# import polygons
rc_use <- sf::st_read(here::here("03_dados", "vetor", "SP_3543907_USO.shp"), quiet = TRUE) %>% 
  dplyr::mutate(CLASSE_USO = as.factor(CLASSE_USO)) %>% 
  sf::st_transform(crs = 4326)

# rasterize polygons
rc_use_rast <- raster::rasterize(x = rc_use, y = ra_rc_agg_mean, field = "CLASSE_USO")
rc_use_rast

# plot
raster::plot(rc_use_rast, col = viridis::viridis(10))
plot(rc_use$geom, add = TRUE)

# package fasterize
# install.packages("fasterize")
library(fasterize) #melhor porque processa mais rápido

# rasterize with fasterize
rc_use_fast <- fasterize::fasterize(sf = rc_use, raster = ra_rc_agg_mean, field = "CLASSE_USO")
rc_use_fast

# plot
raster::plot(rc_use_fast, col = viridis::viridis(10))
plot(rc_use$geom, add = TRUE)

# vectorization points
ra_rc_agg_mean_points = raster::rasterToPoints(ra_rc_agg_mean, spatial = TRUE) %>% 
  sf::st_as_sf()
ra_rc_agg_mean_points

# plot
raster::plot(ra_rc_agg_mean, col = viridis::viridis(10))
plot(ra_rc_agg_mean_points, pch = 20, main = FALSE, add = TRUE)
plot(rc_2019$geom, col = adjustcolor("red", .2), add = TRUE)

# vectorization points
ra_rc_agg_mean_lines = raster::rasterToContour(ra_rc_agg_mean) %>% 
  sf::st_as_sf()
ra_rc_agg_mean_lines

# plot
raster::plot(ra_rc_agg_mean, col = viridis::viridis(10))
plot(ra_rc_agg_mean_lines$geometry, col = "white", pch = 20, lwd = 2, main = FALSE, add = TRUE)
plot(rc_2019$geom, col = adjustcolor("red", .2), add = TRUE)

# vectorization points
rc_use_fast_polygon = raster::rasterToPolygons(rc_use_fast) %>% 
  sf::st_as_sf() %>% 
  dplyr::group_by(layer) %>% 
  summarise()
rc_use_fast_polygon

# plot
raster::plot(rc_use_fast, col = viridis::viridis(10))
plot(rc_use_fast_polygon, col = viridis::viridis(5), pch = 20, lwd = 2, main = FALSE, add = TRUE)

# export raster
raster::writeRaster(x = ra_rc, 
                    filename = here::here("03_dados", "raster", "srtm_27_17_rc"),
                    format = "GTiff",
                    progress = "text",
                    overwrite = TRUE)#apaga o que tem com mesmo nome na pasta

# export stack
raster::writeRaster(x = st_af_crop_mask, 
                    filename = here::here("03_dados", "raster", names(st_af_crop_mask)), 
                    bylayer = TRUE, 
                    format = "GTiff",
                    progress = "text",
                    overwrite = TRUE)

# end ---------------------------------------------------------------------
