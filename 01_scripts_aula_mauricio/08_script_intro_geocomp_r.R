#' ---
#' title: aula 08 - Produção de mapas
#' author: mauricio vancine
#' date: 2020-10-23
#' ---

# topics ------------------------------------------------------------------
# 8.1 Elementos de um mapa
# 8.2 Pacotes para produção de mapas
# 8.3 Pacote ggplot2
# 8.4 Pacote tmap
# 8.5 Mapas vetoriais
# 8.6 Mapas matriciais
# 8.7 Mapas estáticos
# 8.8 Mapas animados
# 8.9 Mapas interativos
# 8.10 Exportar mapas

# packages ----------------------------------------------------------------
library(raster)
library(sf)
library(tidyverse)
library(geobr)
library(rnaturalearth)
library(tmap)
library(viridis)
library(wesanderson)
library(cptcity)

# 8.1 elementos de um mapa
# elementos 
# - titutlo
# - mapa principal
# - mapa secundario
# - legenda
# - coordenadas
# - orientacao (Norte)
# - barra de escala
# - src
# - fontes

# 8.2 pacotes para producao de mapas --------------------------------------
# principais pacotes para producao de mapas no R
# - tmap
# - ggplot2
# - ggspatial
# - cartography
# - googleway
# - leaflet
# - mapview
# - plotly
# - rasterVis

# 8.3 pacote ggplot2 ------------------------------------------------------
# ggplot(...) 
# aes(...) 
# geom_(...) 
# stats_(...) 
# coord_(...) 
# facet_(...) 
# theme_(...)

# 8.4 pacote tmap ---------------------------------------------------------
# tm_shape(...)
# tm_*(...)
# tm_facets(...)
# tm_layout(...)

# tm_style(...)
# tmap_arrange(...)
# tm_mode(...)

# 8.5 mapas vetoriais -----------------------------------------------------
# rio claro
rc_2019 <- geobr::read_municipality(code_muni = 3543907, year = 2019)
rc_2019
plot(rc_2019$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)

# rio claro utm
rc_2019_utm <- geobr::read_municipality(code_muni = 3543907, year = 2019) %>% sf::st_transform(crs = 31983)
rc_2019_utm
plot(rc_2019_utm$geom, col = "gray", main = NA, axes = TRUE, graticule = TRUE)

# import points
rc_spr <- sf::st_read(here::here("03_dados", "vetor", "SP_3543907_NASCENTES.shp"), quiet = TRUE)
rc_spr
plot(rc_spr$geometry, col = "blue", main = NA, axes = TRUE, graticule = TRUE)

# import lines
rc_riv <- sf::st_read(here::here("03_dados", "vetor", "SP_3543907_RIOS_SIMPLES.shp"), quiet = TRUE)
rc_riv
plot(rc_riv$geometry, col = "steelblue", main = NA, axes = TRUE, graticule = TRUE)

# import polygons
rc_use <- sf::st_read(here::here("03_dados", "vetor", "SP_3543907_USO.shp"), quiet = TRUE) %>% sf::st_transform(crs = 31983)
rc_use
plot(rc_use$geometry, col = c("blue", "orange", "gray30", "forestgreen", "green"), main = NA, axes = TRUE, graticule = TRUE)


## ggplot2
# rio claro limit
ggplot() +
  geom_sf(data = rc_2019_utm)

# rio claro limit fill + land use
ggplot() +
  geom_sf(data = rc_2019_utm, color = "black", fill = NA) +
  geom_sf(data = rc_use)

# rio claro limit fill + land use with colors
ggplot() +
  geom_sf(data = rc_2019_utm, color = "black", fill = NA) +
  geom_sf(data = rc_use, aes(fill = CLASSE_USO), color = NA)

# land use with colors + rio claro limit fill
ggplot() +
  geom_sf(data = rc_use, aes(fill = CLASSE_USO), color = NA) +
  geom_sf(data = rc_2019_utm, color = "black", fill = NA) 

# land use and choose colors + rio claro limit fill 
ggplot() +
  geom_sf(data = rc_use, aes(fill = CLASSE_USO), color = NA) +
  geom_sf(data = rc_2019_utm, color = "black", fill = NA) +
  scale_fill_manual(values = c("blue", "orange", "gray30", "forestgreen", "green")) #usado para mudar a cor dos objeivos dentro da aesthetic

# land use and choose colors + rio claro limit fill + coords
ggplot() +
  geom_sf(data = rc_use, aes(fill = CLASSE_USO), color = NA) +
  geom_sf(data = rc_2019_utm, color = "black", fill = NA) +
  scale_fill_manual(values = c("blue", "orange", "gray30", "forestgreen", "green")) +
  coord_sf(datum = 31983)

# land use and choose colors + rio claro limit fill + coords + themes
ggplot() +
  geom_sf(data = rc_use, aes(fill = CLASSE_USO), color = NA) +
  geom_sf(data = rc_2019_utm, color = "black", fill = NA) +
  scale_fill_manual(values = c("blue", "orange", "gray30", "forestgreen", "green")) +
  coord_sf(datum = 31983) +
  theme_bw()

# land use and choose colors + rio claro limit fill + coords + themes + scalebar + north
ggplot() +
  geom_sf(data = rc_use, aes(fill = CLASSE_USO), color = NA) +
  geom_sf(data = rc_2019_utm, color = "black", fill = NA) +
  scale_fill_manual(values = c("blue", "orange", "gray30", "forestgreen", "green")) +
  coord_sf(datum = 31983) +
  theme_bw() +
  annotation_scale(location = "br", width_hint = .3) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0, "cm"), pad_y = unit(.7, "cm"),
                         style = north_arrow_fancy_orienteering)

# land use and choose colors + rio claro limit fill + coords + themes + scalebar + north
ggplot() +
  geom_sf(data = rc_use, aes(fill = CLASSE_USO), color = NA) +
  geom_sf(data = rc_2019_utm, color = "black", fill = NA) +
  scale_fill_manual(values = c("blue", "orange", "gray30", "forestgreen", "green")) +
  coord_sf(datum = 31983) +
  theme_bw() +
  annotation_scale(location = "br", width_hint = .3) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0, "cm"), pad_y = unit(.7, "cm"),
                         style = north_arrow_fancy_orienteering)

# land use and choose colors + rio claro limit fill + coords + themes + scalebar + north + names
map_use_gg <- ggplot() +
  geom_sf(data = rc_use, aes(fill = CLASSE_USO), color = NA) +
  geom_sf(data = rc_2019_utm, color = "black", fill = NA) +
  scale_fill_manual(values = c("blue", "orange", "gray30", "forestgreen", "green")) +
  coord_sf(datum = 31983) +
  theme_bw() +
  annotation_scale(location = "br", width_hint = .3) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0, "cm"), pad_y = unit(.7, "cm"),
                         style = north_arrow_fancy_orienteering) +
  labs(x = "Longitude", y = "Latitude", title = "Cobertura da terra Rio Claro/SP (2015)", fill = "Legenda") +
  theme(legend.position = c(.18,.18),
        legend.box.background = element_rect(colour = "black"),
        axis.text.y = element_text(angle = 90, hjust = .5))
map_use_gg


# springs 
ggplot() +
  geom_sf(data = rc_2019_utm, color = "black", fill = "gray") +
  geom_sf(data = rc_spr, shape = 20, color = "blue") +
  coord_sf(datum = 31983) +
  theme_bw() +
  annotation_scale(location = "br", width_hint = .3) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0, "cm"), pad_y = unit(.7, "cm"),
                         style = north_arrow_fancy_orienteering) +
  labs(x = "Longitude", y = "Latitude", title = "Nascentes Rio Claro/SP (2015)") +
  theme(legend.position = c(.18,.18),
        legend.box.background = element_rect(colour = "black"),
        axis.text.y = element_text(angle = 90, hjust = .5))

# rivers
ggplot() +
  geom_sf(data = rc_2019_utm, color = "black", fill = "gray") +
  geom_sf(data = rc_riv, shape = 20, color = "steelblue") +
  coord_sf(datum = 31983) +
  theme_bw() +
  annotation_scale(location = "br", width_hint = .3) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0, "cm"), pad_y = unit(.7, "cm"),
                         style = north_arrow_fancy_orienteering) +
  labs(x = "Longitude", y = "Latitude", title = "Rios Rio Claro/SP (2015)") +
  theme(legend.position = c(.18,.18),
        legend.box.background = element_rect(colour = "black"),
        axis.text.y = element_text(angle = 90, hjust = .5))



  ## tmap
# rio claro limit
tm_shape(rc_2019_utm) +
  tm_polygons() #se mudar o tm, tipo tm_fill, tm_borders, é possível mudar como o grafico aparece (cores tb)

# rio claro limit fill + land use
tm_shape(rc_2019_utm) +
  tm_polygons() +
  tm_shape(rc_use) +
  tm_polygons()

# rio claro limit fill + land use with colors
tm_shape(rc_2019_utm) +
  tm_polygons() +
  tm_shape(rc_use) +
  tm_polygons(col = "CLASSE_USO")

# land use and choose colors + rio claro limit fill 
tm_shape(rc_use) +
  tm_fill(col = "CLASSE_USO", 
          pal = c("blue", "orange", "gray30", "forestgreen", "green"), title = "Legenda") +
  tm_shape(rc_2019_utm) +
  tm_borders() 

#borda mais grossa no mapa
tm_shape(rc_use) +
  tm_fill(col = "CLASSE_USO", 
          pal = c("blue", "orange", "gray30", "forestgreen", "green"), title = "Legenda") +
  tm_shape(rc_2019_utm) +
  tm_borders(lwd = 2, col = "black") 

# land use and choose colors + rio claro limit fill + coords
tm_shape(rc_2019_utm) +
  tm_polygons() +
  tm_shape(rc_use) +
  tm_polygons(col = "CLASSE_USO", pal = c("blue", "orange", "gray30", "forestgreen", "green"), title = "Legenda") +
  tm_grid(lines = FALSE, labels.format = list(big.mark = ""), labels.rot = c(0, 90))

tm_shape(rc_2019_utm) +
  tm_polygons() +
  tm_shape(rc_use) +
  tm_polygons(col = "CLASSE_USO", pal = c("blue", "orange", "gray30", "forestgreen", "green"), title = "Legenda") +
  tm_grid(lines = FALSE, labels.format = list(big.mark = ""), labels.rot = c(0, 90)) +
  tm_compass() +
  tm_scale_bar()

# land use and choose colors + rio claro limit fill + coords + themes + scalebar + north + names
map_use_tmap <- tm_shape(rc_2019_utm) +
  tm_polygons() +
  tm_shape(rc_use) +
  tm_polygons(col = "CLASSE_USO", pal = c("blue", "orange", "gray30", "forestgreen", "green"), title = "Legenda") +
  tm_grid(lines = FALSE, labels.format = list(big.mark = ""), labels.rot = c(0, 90)) +
  tm_compass() +
  tm_scale_bar() +
  tm_xlab("Longitude") +
  tm_ylab("Latitude") +
  tm_layout(main.title = "Cobertura da terra Rio Claro/SP (2015)")
map_use_tmap

# springs 
tm_shape(rc_2019_utm) +
  tm_polygons() +
  tm_shape(rc_spr) +
  tm_bubbles(col = "HIDRO", pal = "blue", size = .5, alpha = .5) + #col=" " adiciona legend, pal muda a cor
  tm_grid(lines = FALSE, labels.format = list(big.mark = ""), labels.rot = c(0, 90)) +
  tm_compass() +
  tm_scale_bar() +
  tm_xlab("Longitude") +
  tm_ylab("Latitude") +
  tm_layout(main.title = "Nascentes Rio Claro/SP (2015)")

# rivers
tm_shape(rc_2019_utm) +
  tm_polygons() +
  tm_shape(rc_riv) +
  tm_lines(col = "HIDRO", pal = "steelblue") +
  tm_grid(lines = FALSE, labels.format = list(big.mark = ""), labels.rot = c(0, 90)) +
  tm_compass() +
  tm_scale_bar() +
  tm_xlab("Longitude") +
  tm_ylab("Latitude") +
  tm_layout(main.title = "Rios Rio Claro/SP (2015)")

# 8.6 mapas matriciais --------------------------------------------------
# import raster
elev <- raster::raster(here::here("03_dados", "raster", "srtm_27_17_rc.tif")) %>% 
  raster::mask(rc_2019)
elev
plot(elev, col = viridis::viridis(10))

## ggplot2
# raster to tibble
da_elev <- raster::rasterToPoints(elev) %>% 
  tibble::as_tibble() %>% 
  dplyr::rename(elev = srtm_27_17_rc)
head(da_elev)

# elevation map
map_elev_gg <- ggplot() +
  geom_raster(data = da_elev, aes(x = x, y = y, fill = elev)) +
  geom_sf(data = rc_2019, color = "black", fill = NA) +
  scale_fill_gradientn(colors = viridis::viridis(10)) +
  theme_bw() +
  annotation_scale(location = "br", width_hint = .3) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0, "cm"), pad_y = unit(.7, "cm"),
                         style = north_arrow_fancy_orienteering) +
  labs(x = "Longitude", y = "Latitude", title = "Elevação Rio Claro/SP (2007)", fill = "Legenda") +
  theme(legend.position = c(.18,.18),
        legend.box.background = element_rect(colour = "black"),
        axis.text.y = element_text(angle = 90,hjust = 0.5))
map_elev_gg


## tmap
# elevation map
map_elev_tmap <- tm_shape(elev) +
  tm_raster(title = "Legenda") +
  tm_shape(rc_2019) +
  tm_borders() +
  tm_grid(lines = FALSE, labels.format = list(big.mark = ""), labels.rot = c(0, 90)) +
  tm_compass() +
  tm_scale_bar() +
  tm_xlab("Longitude") +
  tm_ylab("Latitude") +
  tm_layout(legend.position = c("left", "bottom"), 
            main.title = "Elevação Rio Claro/SP (2015)")
map_elev_tmap

# elevation map #mudando as cores e tals
map_elev_tmap <- tm_shape(elev) +
  tm_raster(pal = wesanderson::wes_palette("Zissou1"), title = "Legenda") +
  tm_shape(rc_2019) +
  tm_borders() +
  tm_grid(lines = FALSE, labels.format = list(big.mark = ""), labels.rot = c(0, 90)) +
  tm_compass() +
  tm_scale_bar() +
  tm_xlab("Longitude") +
  tm_ylab("Latitude") +
  tm_layout(legend.position = c("left", "bottom"), 
            main.title = "Elevação Rio Claro/SP (2015)")
map_elev_tmap

# elevation map
map_elev_tmap <- tm_shape(elev) +
  tm_raster(pal = cptcity::cpt(pal = "gmt_GMT_dem4"), 
            # breaks = c(400, 500, 600, 700, 800, 900), 
            n = 20, title = "Legenda") + #função break especifica as camadas que está no n 
  tm_shape(rc_2019) +
  tm_borders() +
  tm_grid(lines = FALSE, labels.format = list(big.mark = ""), labels.rot = c(0, 90)) +
  tm_compass() +
  tm_scale_bar() +
  tm_xlab("Longitude") +
  tm_ylab("Latitude") +
  tm_layout(legend.position = c("left", "bottom"), legend.outside = TRUE,
            main.title = "Elevação Rio Claro/SP (2015)")
map_elev_tmap

# 8.8 Mapas animados ------------------------------------------------------
# download
update.packages(ask = F)

year <- geobr::list_geobr()[3, 3] %>% 
  dplyr::pull() %>% 
  stringr::str_split(", ", simplify = TRUE) %>% 
  as.character()
year

br <- NULL

for(i in year){
  
  br <- geobr::read_state("all", i) %>% 
    dplyr::mutate(year = i) %>% 
    dplyr::bind_rows(br, .)
  
}

br$year %>% table

# create facet
br_years <- tm_shape(br) + 
  tm_polygons() + 
  tm_facets(along = "year", free.coords = FALSE)

# export
tmap_animation(br_years, filename = "geo_br_years.gif", delay = 25)

# 8.9 mapas interativos ---------------------------------------------------
# change plot tmap
tmap_mode("view")
map_use_tmap

# change plot tmap
tmap_mode("view")
map_elev_tmap

# 8.10 exportar mapas ---------------------------------------------------
## ggplot2
# export 
ggsave(map_use_gg, 
       filename = "map_rio_claro_land_use_gg.png",
       path = here::here("03_dados"),
       width = 20, 
       height = 20, 
       units = "cm", 
       dpi = 300)

# export
ggsave(map_elev_gg, 
       filename = "map_rio_claro_elevation_gg.png",
       path = here::here("03_dados"),
       width = 20, 
       height = 20, 
       units = "cm", 
       dpi = 300)

## tmap
# export
tmap_save(map_use_tmap, 
          filename = here::here("03_dados", "map_rio_claro_land_use_tmap.png"),
          width = 20, 
          height = 20, 
          units = "cm", 
          dpi = 300)

# export
tmap_save(map_elev_tmap, 
          filename = here::here("03_dados", "map_rio_claro_elevation_tmap.png"),
          width = 20, 
          height = 20, 
          units = "cm", 
          dpi = 300)


# muito obrigado pela paciencia e atencao =]

# end ---------------------------------------------------------------------