rm(list = ls())
gc()

options(verbose = F, warn = F)

require(nhdplusTools)
require(raster)
require(sf)
require(elevatr)
require(concaveman)
require(plyr)
require(stplanr)
require(igraph)
require(smoothr)
#require(exactextractr)

#====================================================================================================
# download critical habitat polygons

# # Download critical habitat polygons
# temp = tempfile(tmpdir = getwd())
# download.file("https://ecos.fws.gov/docs/crithab/crithab_all/crithab_all_shapefiles.zip", temp)
# unzip(temp,  overwrite = TRUE, exdir = file.path('crithab_all_shapefiles'))
# file.remove(temp)

#====================================================================================================
# UNIT TESTING

# polygon example with broad (aquatic) extent:

# load critical habitat polys
file = 'FCH_Notropis_topeka_20040727.zip'

# unzip a shapefile to 'temp_unzip' subdirectory
unzip(file.path('crithab_all_shapefiles', file), exdir = file.path('crithab_all_shapefiles', 'unzipped'))

# read unzipped shapefile
shp = try(st_read(file.path('crithab_all_shapefiles', 'unzipped', gsub('.zip', '.shp', file))))

# dissolve adjoining polygons
shp = st_as_sf(st_cast(st_union(shp), "POLYGON"))


wbd = "WBD_National_GDB.gdb"

if(!file.exists(wbd)){
  download_wbd(".")
}

# load hu8 basins
huc10 = st_read(wbd, layer = 'WBDHU10')



# http://village.anth.wsu.edu
vepPolygon <- polygon_from_extent(raster::extent(672800,740000,4102000,4170000), 
                                  proj4string='+proj=utm +datum=NAD83 +zone=12')

# Get the NHD (USA ONLY)
NHD <- get_nhd(template=vepPolygon, label='VEPIIN', raw.dir = getwd(), extraction.dir = getwd(), force.redo = T)

# Plot the VEP polygon
plot(vepPolygon)

# Plot the NHD data
plot(NHD$NHDFlowline, add=T)
plot(NHD$NHDLine, add=T)
plot(NHD$NHDArea, col='black', add=T)
plot(NHD$NHDWaterbody, col='black', add=T)
# }

poly = shp[1,]


require(FedData)
get_nhd()


# retain hucs that intersect with the habitat polygons
poly_hucs = huc10[apply(st_intersects(huc10, st_transform(poly, st_crs(huc10)), sparse = F), 1, any),]

gsub("(.{2})", "\\1 ", poly_hucs$huc10)

require(ggplot2)
ggplot()+
  geom_sf(data = poly, fill = 'red')+
  geom_sf(data = poly_hucs, fill = 'transparent', colour = 'blue')

# need to get flowlines that intersect with shp bounds

#minAnnualFlowDists(r, shp, search_dist)


# wbd = "WBD_National_GDB.gdb"
# 
# if(!file.exists(wbd)){
#   download_wbd(".")
# }


#=======================================================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# WORKING HERE --- Next try to subset the NHD by the bounding box of each polygon. Then clip to the polygon.
#                  Then identify points at which flowlines intersect the polygon. Identify unique drainage
#                  Apply routing algorthim to each drainage.

download_nhdplusv2(
  './data',
  url = paste0("https://s3.amazonaws.com/edap-nhdplus/NHDPlusV21/",
               "Data/NationalData/NHDPlusV21_NationalData_Seamless", "_Geodatabase_Lower48_07.7z")
)

output_file <- tempfile(fileext = ".gpkg")

# TEST = subset_nhdplus(
#   #output_file = output_file,
#   nhdplus_data = './data/NHDPlusNationalData/NHDPlusV21_National_Seamless_Flattened_Lower48.gdb',
#   bbox = st_bbox(shp[1,]),
#   simplified = FALSE,
#   overwrite = FALSE,
#   return_data = TRUE,
#   status = TRUE,
#   flowline_only = TRUE,
#   streamorder = FALSE,
#   out_prj = 4269
# )

nhd_layers = st_layers('./data/NHDPlusNationalData/NHDPlusV21_National_Seamless_Flattened_Lower48.gdb')

conus_flowlines = st_read('./data/NHDPlusNationalData/NHDPlusV21_National_Seamless_Flattened_Lower48.gdb', layer = "NHDFlowline_Network")




# load hu10 basins
huc10 = st_read(wbd, layer = 'WBDHU10')

TEST = st_intersects(st_transform(shp, st_crs(huc10)), huc10, sparse = TRUE)

# retain hucs that intersect with the habitat polygons
shp_hucs = huc10[apply(st_intersects(huc10, st_transform(shp, st_crs(huc10)), sparse = F), 1, any),]

get_huc8(shp)


# retain hucs that intersect with the habitat polygons
shp_hucs = huc[apply(st_intersects(huc, shp, sparse = F), 1, any),]

ggplot()+
  geom_sf(data = shp, fill = 'red')+
  geom_sf(data = shp_hucs, fill = 'transparent', colour = 'blue')


# disolve adjacent huc polygons
shp_hucs = st_sf(st_union(shp_hucs))