# K4u <- extract(ugrid[[DAY]], newpts3, method = "bilinear")
# K4v <- extract(vgrid[[DAY]], newpts3, method = "bilinear")
# Ku <-(K1u+2*K2u+2*K3u+K4u)/6
# Kv <-(K1v+2*K2v+2*K3v+K4v)/6
#
# newpts <- cbind(newpts[,1] + (metres2lon(Ku,  newpts[,2]) * time_step) ,
#                 newpts[,2] + (metres2lat(Kv) * time_step))
#

library(readxl)
waypoint <- read_excel("waypoint_navigation/waypoint_navigation_matchup.xlsx")
waypoint <- waypoint[-1, ]
library(tidyverse)
library(raster)

pproj <- "+proj=laea +lon_0=80 +lat_0=-67 +ellps=WGS84 +no_defs"



if (FALSE) {
#' @rdname read-maps
read_i_u <- function(i, files) {
  ex <- extent(50, 104, -80, -40)
  dat <- crop(raster(files$fullname[i], varname = "u"), ex, snap = "out")
  raadtools:::.rotate(dat)
}
#' @rdname read-maps
read_i_v <- function(i, files) {
  ex <- extent(50, 104, -80, -40)
  dat <- crop(raster(files$fullname[i], varname = "v"), ex, snap = "out")
  raadtools:::.rotate(dat)
}

ff <- raadtools::currentsfiles() %>% as_tibble() %>%  filter(between(date, min(waypoint$utc), max(waypoint$utc)))

if (!file.exists("ubulk.grd")) {
  ubulk <- aceecostats::build_bulk_file(ff, "ubulk.grd", read_i_raster = read_i_u, layer_prefix = "u_vel")
}
if (!file.exists("vbulk.grd")) {
  vbulk <- aceecostats::build_bulk_file(ff, "vbulk.grd", read_i_raster = read_i_v, layer_prefix = "v_vel")
}

ubulk <- brick("ubulk.grd")
vbulk <- brick("vbulk.grd")

## read everything in as a big table
d <- coordinates(raster(ubulk)) %>% as_tibble() %>% setNames(c("lon", "lat")) %>%
  mutate(cell_ = row_number())

## take out the map projection here
## it just complicates things

d <- local({
  XY <- rgdal::project(cbind(d$lon, d$lat), pproj)
  d$X <- XY[, 1]
  d$Y <- XY[, 2]
  d
})
dim(d)
d <- bind_rows(replicate(nlayers(ubulk), d, simplify = FALSE)) %>%
  mutate(date = rep(as.numeric(as.POSIXct(getZ(ubulk), tz = "GMT")), each = ncell(ubulk)),
         u = values(readAll(ubulk)), v = values(readAll(vbulk)))

saveRDS(d, "bulk_uv.rds")

}
d <- readRDS("bulk_uv.rds")

library(nabor)
nn_fun <- nabor::WKNNF(d %>% dplyr::select(X, Y, date) %>%
                         mutate(date = as.numeric(date)) %>%
                         as.matrix())


nn_fun$query(cbind(147, -42, 16812), k = 1, eps = 0)

find_nn <- function(x, ...) {
  UseMethod("find_nn")
}
find_nn.matrix <- function(x, ...) {
  nn_fun$query(x, k = 1, eps = 0)
}

find_nn.data.frame <- function(x, ...) {
  tibble::as_tibble(lapply(find_nn(as.matrix(x %>% dplyr::select(X, Y, TIME))), c))
}

tr <- waypoint %>% transmute(LONGITUDE, LATITUDE,
                             TIME = as.numeric(utc))

curr_track <- local({
  XY <- rgdal::project(cbind(tr$LONGITUDE, tr$LATITUDE), pproj)
  tr %>% mutate(X = XY[, 1], Y = XY[, 2])
})

n_iter <- 5000
trim <- 100
time_step <- 60
l <- vector("list", n_iter/trim)
count <- 0
for (i in seq_len(n_iter)) {

  idx <- find_nn(curr_track)
  uv_incr <- d %>% slice(idx[["nn.idx"]]) %>% dplyr::select(uu = u, vv = v) %>% mutate(n = row_number())
  tab <- inner_join(curr_track %>% mutate(n = row_number()), uv_incr, "n")
  tab$X = tab$X + tab$uu * time_step
  tab$Y = tab$Y + tab$vv * time_step
  tab$TIME = tab$TIME + time_step

  curr_track <- tab %>% dplyr::select(X, Y, TIME)
  if (i %% trim == 0) {
    count <- count + 1
    print(count)
    l[[count]] <- curr_track
  }
}

dd <- bind_rows(l)
dd <- local({
  XY <- rgdal::project(cbind(dd$X, dd$Y), pproj, inv = TRUE)
  dd$lon <- XY[,1]
  dd$lat <- XY[, 2]
  dd
})
library(ggplot2)
ggplot(dd, aes(lon, lat, colour = TIME)) + geom_point(pch = ".")




###############################################################################################
## utility functions
metres2lon <- function(m, lat, radius = 6378137.0) {
  (m * cos
   (lat * pi/180)) / (2 * pi * radius / 360)
}
metres2lat <- function(m, radius = 6378137.0) {
  m / (2 * pi * radius / 360)
}
raad_ufield <- function(date, ...) {

}
get_u_field <- function(date, ...) {
  files <- list.files(path, pattern = sprintf("%s_.*grd$", option), full.names = TRUE)
  filename <- grep(sprintf("%s_Dec", as.character(year)), files, value = TRUE)
  brick(filename)
}
