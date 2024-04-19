# function to divide a polygon into subregions of roughly equal size
# steps: random sampling of points, clustering via kmeans, then voronoi tesselation
# approach taken from here: https://gis.stackexchange.com/questions/375345/dividing-polygon-into-parts-which-have-equal-area-using-r
# type_sampling is random by default, but all the options in st_sample can be chosen

library(dismo)
library(sf)
library(terra)
library(rnaturalearth)
library(deldir)

divide.into.regions = function(boundary_regions, nbregions, nbsamples = NULL, type_sampling = "random", itermax_kmeans = 50){
  if(is.null(nbsamples)) nbsamples = nbregions * 1000
  boundary_regions = st_as_sf(boundary_regions[,1])
  points_rd = st_sample(boundary_regions, size = nbsamples, type = type_sampling)
  points_rd = st_coordinates(points_rd)
  k_means = kmeans(points_rd, centers = nbregions, iter.max = itermax_kmeans,nstart = 10)
  if(k_means$ifault == 4){
    k_means = kmeans(points_rd, centers = nbregions, iter.max = itermax_kmeans,nstart = 10,algorithm = "MacQueen")
  }
  polygons_voronoi = dismo::voronoi(k_means$centers, ext = boundary_regions)
  polygons_voronoi = st_as_sf(polygons_voronoi)
  st_crs(polygons_voronoi) = st_crs(boundary_regions)
  polygons_regions = vect(st_intersection(polygons_voronoi, boundary_regions))
  polygons_regions = polygons_regions[,1]
  return(polygons_regions)
}

# example with UK
borders_UK = vect(ne_countries(scale = 110, type = "countries",country = "United Kingdom",returnclass = "sf"))
plot(borders_UK)

# split into 100 roughly equal areas
# standard deviation is typically 10% of mean area
# random sampling
nbregions = 100
UK_clusters = divide.into.regions(borders_UK, nbregions = nbregions,itermax_kmeans = 10) # can take a while
UK_clusters$area = expanse(UK_clusters)/1000000 # in sqkm
plot(UK_clusters,"area", type = "continuous", col = viridis(100))
mean(UK_clusters$area)
sd(UK_clusters$area)

# # regular sampling
# UK_clusters = divide.into.regions(borders_UK, nbregions = nbregions,itermax_kmeans = 10, type_sampling = "regular") # can take a while
# UK_clusters$area = expanse(UK_clusters)/1000000 # in sqkm
# plot(UK_clusters,"area", type = "continuous", col = viridis(100))
# mean(UK_clusters$area)
# sd(UK_clusters$area)




