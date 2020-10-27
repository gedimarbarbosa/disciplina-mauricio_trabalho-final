#' ---
#' title: aula 06 - estrutura e manejo de dados vetoriais 
#' author: mauricio vancine
#' date: 2020-10-22
#' ---

# topics ------------------------------------------------------------------
# 6.1 pacotes
# 6.2 geometrias sf
# 6.3 classes sf
# 6.4 importar dados vetoriais
# 6.5 descricao de objetos sf
# 6.6 converter dados para sf
# 6.7 converter crs
# 6.8 operacoes de atributos
# 6.9 operacoes espaciais
# 6.10 operacoes geometricas
# 6.11 exportar dados vetoriais

# packages ----------------------------------------------------------------
library(sp)
library(sf)
library(geobr)
library(rnaturalearth) #pacote com dados globais
library(tidyverse)

# 6.1 pacotes -------------------------------------------------------------
# sp
# install.packages("sp")
# library(sp)

# sf
# install.packages("sf")
# library(sf)

# 6.3 classes sf ----------------------------------------------------------
# simple feature geometries (sfg)

# simple
# sf::st_point()
# sf::st_linestring()
# sf::st_polygon()

# multi
# sf::st_multipoint()
# sf::st_multilinestring()
# sf::st_multipolygon()

# collections
# sf::st_geometrycollection()

# vector - point
vec <- c(5, 2)
vec

po <- sf::st_point(vec)
po

# plot
plot(po, pch = 20, cex = 4, axes = TRUE, graticule = TRUE)

# matrix - multipoint
multipoint_matrix <-  rbind(c(5, 2), c(1, 3), c(3, 4), c(3, 2))
multipoint_matrix

po_mul <- sf::st_multipoint(multipoint_matrix)
po_mul

# plot
plot(po_mul, pch = 20, cex = 4, axes = TRUE, graticule = TRUE)

# matrix - multipoint
multipoint_matrix <- rbind(c(5, 2), c(1, 3), c(3, 4), c(3, 2))
multipoint_matrix

lin <- sf::st_linestring(multipoint_matrix)
lin

# plot
plot(lin, lwd = 2, axes = TRUE, graticule = TRUE)

# list - polygon
polygon_list <- list(rbind(c(1, 5), c(2, 2), c(4, 1), c(4, 4), c(1, 5)))
polygon_list

pol <- sf::st_polygon(polygon_list)
pol

# plot
plot(pol, col = "gray", axes = TRUE, graticule = TRUE)

# simple feature columns (sfc)
# sf::st_sfc()
# sf::st_geometry_type()
# sf::st_crs()

# sfc point
point1 <- sf::st_point(c(5, 2))
point1

point2 <- sf::st_point(c(1, 3))
point2

points_sfc <- sf::st_sfc(point1, point2)
points_sfc

sf::st_geometry_type(points_sfc)

# plot
plot(points_sfc, pch = 20, cex = 4, axes = TRUE, graticule = TRUE)

# sfc geometry
po_pol_sfc <- sf::st_sfc(po, pol)
po_pol_sfc

sf::st_geometry_type(po_pol_sfc)

# plot
plot(po_pol_sfc, pch = 20, cex = 4, lwd = 2, col = "gray", axes = TRUE, graticule = TRUE)

# epgs definition (codigo para o sistema de datum)
points_sfc_wgs <- sf::st_sfc(point1, point2, crs = 4326)
points_sfc_wgs

sf::st_crs(points_sfc_wgs)

# proj4string definition
points_sfc_wgs <- sf::st_sfc(point1, point2, crs = "+proj=longlat +datum=WGS84 +no_defs")
points_sfc_wgs

sf::st_crs(points_sfc_wgs)

# plot
plot(points_sfc_wgs, pch = 20, cex = 4, axes = TRUE, graticule = TRUE)

# class sf
rc_point <- sf::st_point(c(-47.57,-22.39))         # sfg object
rc_point

rc_geom <- sf::st_sfc(rc_point, crs = 4326)        # sfc object
rc_geom

rc_attrib <- data.frame(                       # data.frame object
  name = "Rio Claro",
  temperature = 19,
  date = as.Date("2020-10-13")
)
rc_attrib

rc_sf <- sf::st_sf(rc_attrib, geometry = rc_geom)  # sf object
rc_sf

class(rc_sf)

# plot
plot(rc_sf[1], pch = 20, cex = 4, axes = TRUE, graticule = TRUE) # rc_sf[1] - plotar a primeira coluna

# plot
plot(rc_sf[2], pch = 20, cex = 4, axes = TRUE, graticule = TRUE) # rc_sf[2] - plotar a segunda coluna

# plot
plot(rc_sf$geometry, pch = 20, cex = 4, axes = TRUE, graticule = TRUE) # rc_sf$geometry - plotar apenas a geometria

# 6.4 importar dados ------------------------------------------------------
# create directory
dir.create(here::here("03_dados", "vetor"))

# increase time to download #para o R nao cancelar
options(timeout = 600)

# download points
for(i in c(".dbf", ".prj", ".shp", ".shx")){
  download.file(url = paste0("http://geo.fbds.org.br/SP/RIO_CLARO/HIDROGRAFIA/SP_3543907_NASCENTES", i),
                destfile = here::here("03_dados", "vetor", paste0("SP_3543907_NASCENTES", i)), mode = "wb")
}

# download lines
for(i in c(".dbf", ".prj", ".shp", ".shx")){
  download.file(url = paste0("http://geo.fbds.org.br/SP/RIO_CLARO/HIDROGRAFIA/SP_3543907_RIOS_SIMPLES", i),
                destfile = here::here("03_dados", "vetor", paste0("SP_3543907_RIOS_SIMPLES", i)), mode = "wb")
}

# download polygons
for(i in c(".dbf", ".prj", ".shp", ".shx")){
  download.file(url = paste0("http://geo.fbds.org.br/SP/RIO_CLARO/USO/SP_3543907_USO", i),
                destfile = here::here("03_dados", "vetor", paste0("SP_3543907_USO", i)), mode = "wb")
}

# import points
rc_spr <- sf::st_read(here::here("03_dados", "vetor", "SP_3543907_NASCENTES.shp"), quiet = TRUE)
rc_spr

# plot
plot(rc_spr$geometry, pch = 20, col = "blue", main = NA, axes = TRUE, graticule = TRUE)

# import lines
rc_riv <- sf::st_read(here::here("03_dados", "vetor", "SP_3543907_RIOS_SIMPLES.shp"), quiet = TRUE)
rc_riv

# plot
plot(rc_riv$geometry, col = "steelblue", main = NA, axes = TRUE, graticule = TRUE)

# import polygons
rc_use <- sf::st_read(here::here("03_dados", "vetor", "SP_3543907_USO.shp"), quiet = TRUE)
rc_use

# plot
plot(rc_use[5], col = c("blue", "orange", "gray30", "forestgreen", "green"), main = NA, axes = TRUE, graticule = TRUE) #definido pela coluna 5 - [5]

# import gps data
gps_gpx <- sf::read_sf(here::here("03_dados", "vetor", "waypoints.gpx"), layer = "waypoints")
gps_gpx

# plot
plot(gps_gpx$geometry, cex = 4, pch = 20, col = "red", main = NA, axes = TRUE, graticule = TRUE)

# import gps data
gps_kml <- sf::read_sf(here::here("03_dados", "vetor", "waypoints.kml"))
gps_kml

# plot
plot(gps_kml$geometry, cex = 4, pch = 20, col = "red", main = NA, axes = TRUE, graticule = TRUE)

# import data from packages
# brazil 2019
br_2019 <- geobr::read_country(year = 2019, showProgress = FALSE)
br_2019

# plot
plot(br_2019$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)

# brazil 1872
br_1872 <- geobr::read_country(year = 1872, showProgress = FALSE)
br_1872

# plot
plot(br_1872$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)

# sao paulo state
sp_2019 <- geobr::read_state(code_state = "SP", year = 2019, showProgress = FALSE)
sp_2019

# plot
plot(sp_2019$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)

# sao paulo municipalities
sp_mun_2019 <- geobr::read_municipality(code_muni = "SP", year = 2019, showProgress = FALSE)
sp_mun_2019

# plot
plot(sp_mun_2019$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)

# rio claro
rc_2019 <- geobr::read_municipality(code_muni = 3543907, year = 2019, showProgress = FALSE)#puxa pelo codigo
rc_2019

# plot
plot(rc_2019$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)

# biomes
bi_2019 <- geobr::read_biomes(year = 2019, showProgress = FALSE)
bi_2019

# plot
plot(bi_2019$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)

# list all datasets available in the geobr package
geobr::list_geobr()

# south america
sa <- rnaturalearth::ne_countries(scale = "small", continent = "South America", returnclass = "sf")
sa

# plot
plot(sa$geometry, col = "gray", main = NA, axes = TRUE, graticule = TRUE)

# coast lines
coastline <- rnaturalearth::ne_coastline(scale = "small", returnclass = "sf")
coastline

# plot
plot(coastline$geometry, col = "blue", main = NA)

# 6.5 descricao de objetos sf ---------------------------------------------
# south america
sa

# geometry
sf::st_geometry_type(sa)

# extention
sf::st_bbox(sa)

# coordinate reference system
sf::st_crs(sa)

# acessar a tabela de atributos
sa_tab <- sf::st_drop_geometry(sa)
sa_tab

# classe
class(sa_tab)

# 6.6 converter dados sf --------------------------------------------------
# import table
si <- readr::read_csv(here::here("03_dados", "tabelas", "ATLANTIC_AMPHIBIANS_sites.csv"))
si

# add column
si <- si %>% 
  dplyr::mutate(lon = longitude, lat = latitude, .before = 1)
si

# convert to sf
si_ve <- si %>% 
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326)
si_ve

# plot
plot(si_ve$geometry, pch = 20, main = NA, axes = TRUE, graticule = TRUE)

# countries sp
co110_sp <- rnaturalearth::countries110
co110_sp

# countries sf
co110_sf <- sf::st_as_sf(co110_sp)
co110_sf

# plot
plot(co110_sf$geometry, col = "gray", main = NA)

# countries sp
co110_sp <- sf::as_Spatial(co110_sf)
co110_sp

# 6.7 conversao do crs --------------------------------------------------
# crs local
# convert coordinate system
rc_2019_sirgas2000_utm23s <- sf::st_transform(rc_2019, crs = 31983)
rc_2019_sirgas2000_utm23s

# convert datum and coordinate system
rc_2019_wgs84_utm23s <- sf::st_transform(rc_2019, crs = 32723)
rc_2019_wgs84_utm23s

# convert datum
rc_2019_wgs84_gcs <- sf::st_transform(rc_2019, crs = 4326)
rc_2019_wgs84_gcs

# crs global
# countries - WGS84/GCS
co110_sf

# plot
plot(co110_sf$geometry, col = "gray", graticule = TRUE)

# mollweide projection
co110_sf_moll <- sf::st_transform(co110_sf, crs = "+proj=moll")
co110_sf_moll

# plot
plot(co110_sf_moll$geometry, col = "gray", graticule = TRUE)

# winkel tripel projection
co110_sf_wintri <- lwgeom::st_transform_proj(co110_sf, crs = "+proj=wintri")
co110_sf_wintri

# plot
plot(co110_sf_wintri$geometry, col = "gray")

# eckert iv projection
co110_sf_eck4 <- sf::st_transform(co110_sf, crs = "+proj=eck4")
co110_sf_eck4

# plot
plot(co110_sf_eck4$geometry, col = "gray", graticule = TRUE)

# lambert projection
co110_sf_laea1 <- sf::st_transform(co110_sf, crs = "+proj=laea +x_0=0 +y_0=0 +lon_0=0 +lat_0=0")
co110_sf_laea1

# plot
plot(co110_sf_laea1$geometry, col = "gray", graticule = TRUE)

# lambert projection (america)
co110_sf_laea2 <- sf::st_transform(co110_sf, crs = "+proj=laea +x_0=0 +y_0=0 +lon_0=-50 +lat_0=0")
co110_sf_laea2

# plot
plot(co110_sf_laea2$geometry, col = "gray", graticule = TRUE)

# 6.8 operacoes de atributos ----------------------------------------------
# 1. attribute subsetting - nesse exemplo estamos selecionando formação florestal
rc_use_forest <- rc_use %>% 
  dplyr::filter(CLASSE_USO == "formação florestal")
rc_use_forest

# plot
plot(rc_2019_sirgas2000_utm23s$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)
plot(rc_use_forest$geometry, col = "forestgreen", add = TRUE)

# 2. attribute joining
# create data
da_class <- tibble::tibble(CLASSE_USO = rc_use$CLASSE_USO, 
                           classe = c("agua", "antropico", "edificado", "floresta", "silvicultura"))
da_class

# attribute joining
rc_use_class <- dplyr::left_join(rc_use, da_class, by = "CLASSE_USO") %>% 
  sf::st_drop_geometry()
rc_use_class

# 3. attribute aggregation
rc_spr_n <- rc_spr %>% 
  dplyr::group_by(MUNICIPIO, HIDRO) %>% 
  dplyr::summarise(n = n())
rc_spr_n

# plot
plot(rc_2019_sirgas2000_utm23s$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)
plot(rc_spr_n$geometry, pch = 20, col = "blue", add = TRUE)

# 4. attribute create - columns
rc_use_use_area <- rc_use %>% 
  dplyr::mutate(classe_area = paste0(CLASSE_USO, " (", AREA_HA, " ha)")) %>% 
  sf::st_drop_geometry()
rc_use_use_area

# 4. attribute create - area
rc_use_area <- rc_use %>% 
  dplyr::mutate(area_m2 = sf::st_area(rc_use),
                area_ha = sf::st_area(rc_use)/1e4) %>% #transformando m2 em ha
  sf::st_drop_geometry()
rc_use_area

#selecionando nascentes aleatoriamente
plot(dplyr::sample_frac(rc_spr, .3)[1], pch = 20)

# feature manipulation
# dplyr::filter()
# dplyr::distinc()
# dplyr::slice()
# dplyr::n_sample()
# dplyr::group_by()
# dplyr::summarise()

# attribute table manipulation
# dplyr::select()
# dplyr::pull()
# dplyr::rename()
# dplyr::mutate()
# dplyr::arrange()
# dplyr::*_join()

# 6.9 operacoes espaciais -------------------------------------------------
# 1. spatial subsetting
sf::st_intersects(x = rc_spr, y = rc_use_forest)

# 1. spatial subsetting
sf::st_intersects(x = rc_spr, y = rc_use_forest, sparse = FALSE)

# 1. spatial subsetting - inside
rc_spr_forest <- rc_spr %>% 
  dplyr::filter(sf::st_intersects(x = ., y = rc_use_forest, sparse = FALSE)) # . dentro do pipe repete objeto rc_spr
rc_spr_forest

# 1. spatial subsetting - inside
rc_spr_forest <- rc_spr[rc_use_forest, ]
rc_spr_forest

# plot
plot(rc_2019_sirgas2000_utm23s$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)
plot(rc_use_forest$geometry, col = "forestgreen", add = TRUE)
plot(rc_spr_forest$geometry, col = "blue", pch = 20, cex = 1, add = TRUE)

# 1. spatial subsetting - outside
rc_spr_forest_out <- rc_spr %>% 
  dplyr::filter(!sf::st_intersects(x = ., y = rc_use_forest, sparse = FALSE)) # indica não
rc_spr_forest_out

# plot
plot(rc_2019_sirgas2000_utm23s$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)
plot(rc_use_forest$geometry, col = "forestgreen", add = TRUE)
plot(rc_spr_forest_out$geometry, col = "steelblue", pch = 20, cex = 1, add = TRUE)
plot(rc_spr_forest$geometry, col = "blue", pch = 20, cex = 1, add = TRUE)

# 2. spatial join
rc_spr_use <- rc_spr %>% 
  sf::st_join(x = ., y = rc_use)
rc_spr_use

# 2. spatial join
rc_spr_use %>% 
  sf::st_drop_geometry() %>% 
  dplyr::count(CLASSE_USO)

# 2. attribute join
col <- tibble::tibble(CLASSE_USO = c("água", "área antropizada", "área edificada", "formação florestal", "silvicultura"),
                      color = c("blue", "orange", "gray30", "forestgreen", "green"))
col

# 2. attribute join
rc_spr_use <- dplyr::left_join(rc_spr_use, col, by = "CLASSE_USO")
rc_spr_use

# plot
plot(rc_2019_sirgas2000_utm23s$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)
plot(rc_spr_use[10], col = rc_spr_use$color, pch = 20, add = TRUE)

# 3. random points
rc_2019_sirgas2000_utm23s_rp <- sf::st_sample(rc_2019_sirgas2000_utm23s, 1e3)
rc_2019_sirgas2000_utm23s_rp

# plot
plot(rc_2019_sirgas2000_utm23s$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)
plot(rc_2019_sirgas2000_utm23s_rp, pch = 20, cex = .5, col = "red", add = TRUE)

# 3. random points
rc_use_forest_rp <- sf::st_sample(rc_use_forest, 1e3)
rc_use_forest_rp

# plot
plot(rc_2019_sirgas2000_utm23s$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)
plot(rc_use_forest$geom, col = "forestgreen", add = TRUE)
plot(rc_use_forest_rp, pch = 20, cex = .5, col = "red", add = TRUE)

# 4. grid
rc_2019_sirgas2000_utm23s_grid <- rc_2019_sirgas2000_utm23s %>% #funciona melhor se estiver em UTM
  sf::st_make_grid(cellsize = 2000) %>% 
  sf::st_as_sf()
rc_2019_sirgas2000_utm23s_grid

# plot
plot(rc_2019_sirgas2000_utm23s$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)
plot(rc_2019_sirgas2000_utm23s_grid, col = adjustcolor("red", .1), add = TRUE)

# 4. grid - subset
rc_2019_sirgas2000_utm23s_grid_in <- rc_2019_sirgas2000_utm23s_grid[rc_2019_sirgas2000_utm23s, ]
rc_2019_sirgas2000_utm23s_grid_in

# plot
plot(rc_2019_sirgas2000_utm23s$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)
plot(rc_2019_sirgas2000_utm23s_grid_in, col = adjustcolor("red", .1), add = TRUE)

# 4. grid - count spring by cell
rc_2019_sirgas2000_utm23s_grid_in_spr_count <- rc_2019_sirgas2000_utm23s_grid_in %>% 
  dplyr::mutate(n = sf::st_intersects(x = ., rc_spr) %>% lengths())
rc_2019_sirgas2000_utm23s_grid_in_spr_count

# map
plot(rc_2019_sirgas2000_utm23s_grid_in_spr_count["n"], axes = TRUE, graticule = TRUE)

# 5. hexagon
rc_2019_sirgas2000_utm23s_hex <- rc_2019_sirgas2000_utm23s %>% 
  sf::st_make_grid(cellsize = 2000, square = FALSE) %>% 
  sf::st_as_sf()
rc_2019_sirgas2000_utm23s_hex

# plot
plot(rc_2019_sirgas2000_utm23s$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)
plot(rc_2019_sirgas2000_utm23s_hex, col = adjustcolor("blue", .1), add = TRUE)

# 5. hexagon - count spring by hexagon
rc_2019_sirgas2000_utm23s_hex_spr_count <- rc_2019_sirgas2000_utm23s_hex %>% 
  dplyr::mutate(n = sf::st_intersects(x = ., rc_spr) %>% lengths())
rc_2019_sirgas2000_utm23s_hex_spr_count

# map
plot(rc_2019_sirgas2000_utm23s_hex_spr_count["n"], axes = TRUE, graticule = TRUE)

# 6.10 operacoes geometricas ----------------------------------------------
# 1. simplification
rc_riv_simp <- sf::st_simplify(rc_riv, dTolerance = 1e4)
rc_riv_simp

# plot
plot(rc_2019_sirgas2000_utm23s$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)
plot(rc_riv$geometry, col = "steelblue", add = TRUE)
plot(rc_riv_simp$geometry, col = "red", add = TRUE)

# 2. centroids
rc_2019_cent <- sf::st_centroid(rc_2019_sirgas2000_utm23s)
rc_2019_cent

# plot
plot(rc_2019_sirgas2000_utm23s$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)
plot(rc_2019_cent$geom, cex = 3, col = "red", pch = 20, add = TRUE)

# 2. centroids
rc_2019_sirgas2000_utm23s_hex_cent <- sf::st_centroid(rc_2019_sirgas2000_utm23s_hex)
rc_2019_sirgas2000_utm23s_hex_cent

# plot
plot(rc_2019_sirgas2000_utm23s$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)
plot(rc_2019_sirgas2000_utm23s_hex, col = adjustcolor("blue", .1), add = TRUE)
plot(rc_2019_sirgas2000_utm23s_hex_cent, col = "blue", pch = 20, add = TRUE)

# 3. buffer
rc_spr_forest_buf <- rc_spr_forest %>% 
  sf::st_buffer(dist = 500)
rc_spr_forest_buf

# plot
plot(rc_2019_sirgas2000_utm23s$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)
plot(rc_spr_forest$geom, pch = 20, cex = .5, col = "blue", add = TRUE)
plot(rc_spr_forest_buf$geom, col = adjustcolor("white", 0), add = TRUE)

# 4. union
rc_spr_forest_buf_union <- sf::st_union(rc_spr_forest_buf)
rc_spr_forest_buf_union

# plot
plot(rc_2019_sirgas2000_utm23s$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)
plot(rc_spr_forest$geom, pch = 20, cex = .5, col = "blue", add = TRUE)
plot(rc_spr_forest_buf_union, col = adjustcolor("white", 0), add = TRUE)

# 5. clipping - intersection
rc_riv_spr_forest_buf_int <- sf::st_intersection(x = rc_riv, y = rc_spr_forest_buf_union)
rc_riv_spr_forest_buf_int

# plot
plot(rc_2019_sirgas2000_utm23s$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)
plot(rc_spr_forest_buf_union, col = adjustcolor("blue", .1), add = TRUE)
plot(rc_riv_spr_forest_buf_int$geometry, col = "blue", add = TRUE)

# 5. clipping - difference
rc_riv_spr_forest_buf_dif <- sf::st_difference(x = rc_riv, y = rc_spr_forest_buf_union)
rc_riv_spr_forest_buf_dif

# plot
plot(rc_2019_sirgas2000_utm23s$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)
plot(rc_spr_forest_buf_union, col = adjustcolor("blue", .1), add = TRUE)
plot(rc_riv_spr_forest_buf_dif$geometry, col = "blue", add = TRUE)

# 6.11 exportar dados -----------------------------------------------------
# export shapefile
sf::write_sf(obj = rc_spr_forest_buf, 
             dsn = here::here("03_dados", "vetor", "rc_spr_forest_buf.shp"),
             layer = "rc_spr_forest_buf.shp",
             driver = "ESRI Shapefile")

#geopackage integra todos os shapes em um único arquivo - mais fácil e leve pra trabalhar
# export geopackage
sf::write_sf(obj = rc_spr_forest_buf, 
             dsn = here::here("03_dados", "vetor", "rc_spr_forest_buf.gpkg"),
             layer = "rc_spr_forest_buf")

# export geopackage
sf::write_sf(obj = rc_spr_forest_buf_union, 
             dsn = here::here("03_dados", "vetor", "rc_spr_forest_buf.gpkg"),
             layer = "rc_spr_forest_buf_union")

# end ---------------------------------------------------------------------
