rm(list = ls())

#====================================================================================================
# load required packages, installing any that have not yet been installed

packages = c(
  "nhdplusTools",
  "raster",
  "sf"
)

install.packages(setdiff(packages, rownames(installed.packages()))) 
lapply(packages, require, character.only = TRUE)


options(verbose = F, warn = F)

#====================================================================================================
# load functions

source('HUCconversionArea.R')
source('bufferConversionArea.R')

#====================================================================================================
# download s35 raster (first crop)

# Authenticate BOX connection using Client ID and Secret in renviron
source('box_auth.R')

# download the usxp s35 mtr raster from BOX
box_dl(file_id = '778783353386', overwrite = T, file_name = 's35_fc.tif')

# load the usxp s35 mtr raster (5 land use change classes)
r = raster("s35_fc.tif")

#====================================================================================================
# download critical habitat polygons

# Download critical habitat polygons
temp = tempfile(tmpdir = getwd())
download.file("https://ecos.fws.gov/docs/crithab/crithab_all/crithab_all_shapefiles.zip", temp)
unzip(temp,  overwrite = TRUE, exdir = file.path('crithab_all_shapefiles'))
file.remove(temp)

files = list.files(path = 'crithab_all_shapefiles')

#====================================================================================================
# download watershed boundaries

wbd = "WBD_National_GDB.gdb"

if(!file.exists(wbd)){
  download_wbd(".")
}

# get HUC 10 and HUC 12 Boundaries
huc10 = st_read(wbd, layer = 'WBDHU10')
huc12 = st_read(wbd, layer = 'WBDHU12')

hucs = list('huc10' = huc10, 'huc12' = huc12)

#====================================================================================================
#====================================================================================================

# select crops for specific tabulations
crops = c(1, 5)
buffers = c(0, 1.6, 16) # km

# create list in which to accumulate results
outList = list()

# iterate through each critical habitat file
for(file in files){
  
  cat('file ', match(file, files), ' of ', length(files), '\n')
  
  # create name for outlist element name
  name = gsub('.zip', '', file)
  
  # unzip a shapefile to 'temp_unzip' subdirectory
  unzip(file.path('crithab_all_shapefiles', file), exdir = file.path('crithab_all_shapefiles', 'unzipped'))
  
  # read unzipped shapefile
  shp = try(st_read(file.path('crithab_all_shapefiles', 'unzipped', gsub('.zip', '.shp', file))))
  
  # if error occurs opening shp, create NULL output
  if(is(shp, "try-error")){
    
    outList[[name]] = NULL
    
  }else{
  
    # calculate conversion area within hucs
    huc_areas = try(HUCconversionArea(r, shp, hucs, crops))
    if(is(huc_areas, "try-error")) huc_areas = NULL else huc_areas = huc_areas
    
    # calculate conversion area within buffers
    buffer_areas = try(bufferConversionArea(r, shp, buffers, crops))
    if(is(buffer_areas, "try-error")) buffer_areas = NULL else buffer_areas = buffer_areas
    
    outList[[name]] = c(huc_areas, buffer_areas)
      
  }
  
}

# merge into data frame, matching column names
allCols = unique(unlist(lapply(outList, function(v) names(v))))
outDF = as.data.frame(do.call(rbind, lapply(outList, function(x) x[match(allCols, names(x))])))
outDF = cbind(CH_file = row.names(outDF), outDF)
row.names(outDF) = NULL

write.csv(outDF, 'CH_huc&buffer_analysis.csv', row.names = F)

#====================================================================================================
# merge with species attributs from 'ECOS USFWS Threatened  Endangered Species Active Critical Habitat Report.csv'
# downloaded manually from https://ecos.fws.gov/ecp/report/table/critical-habitat.html

ch_attribs = read.csv('ECOS USFWS Threatened  Endangered Species Active Critical Habitat Report.csv')
ch_attribs['key'] = gsub(" ", "_", ch_attribs$Scientific.Name)


outDF['key'] = gsub('NCH_', '', gsub('PCH_', "", gsub('FCH_', "", gsub('_$', '', gsub('[[:digit:]]', '', outDF$CH_file)))))

outDF = merge(outDF, ch_attribs, by = 'key', all = T)

write.csv(outDF, 'criticalHabitat_huc&buffer_output.csv', row.names = T)
  
# # # load critical habitat polys
# # file = 'FCH_Notropis_topeka_20040727.zip'
# 
# 
# # load critical habitat polys
# file = 'FCH_Notropis_topeka_20040727.zip'
# 
# # unzip a shapefile to 'temp_unzip' subdirectory
# unzip(file.path('crithab_all_shapefiles', file), exdir = file.path('crithab_all_shapefiles', 'unzipped'))
# 
# # read unzipped shapefile
# shp = try(st_read(file.path('crithab_all_shapefiles', 'unzipped', gsub('.zip', '.shp', file))))
# 
# bufferConversionArea(r, shp, buffers = c(0, 1.6, 16), crops)
# require(units)
# 
# # reproject to match raster (puts units in m for more transparent buffering)
# shp = st_transform(shp, st_crs(r))
# 
# # buffer and disolve
# shp_buff = st_buffer(shp, set_units(1.6, km))
# shp_buff = st_as_sf(st_cast(st_union(shp_buff), "POLYGON"))
# 
# # count all conversion pixels (value != 0) in each huc; sum 
# all_conv = sum(exact_extract(
#   x = r,
#   y = shp_buff,
#   fun = function(value, cov_frac) sum(ifelse(value != 0,1,0)*cov_frac, na.rm = T)
# ), na.rm = T)
# 
# # iterate through each crop 
# crop_conv = lapply(crops, function(crop, r, huc_shp){
#   
#   # count all conversion pixels of given crop in huc; sum
#   conv = sum(exact_extract(
#     x = r,
#     y = shp_buff,
#     fun = function(value, cov_frac) sum(ifelse(value != 0,1,0)*cov_frac, na.rm = T)
#   ), na.rm = T)
#   
#   return(conv)
#   
# }, r = r, huc_shp = huc_shp)
# 
# # name elements of huc_crop_conv output
# names(crop_conv) = crops
# crop_conv = unlist(crop_conv)
# 
# # calculate area of huc polygons (units taken from raster crs)
# buff_area = sum(st_area(shp_buff), na.rm = T) # units are those of r: m
# 
# # prepare output: convert conversion area to meters; then everthing to hectares
# out = c('all' = all_conv, crop_conv)*prod(res(r)) # convert to m2
# out = c(out, 'area' = huc_area)*0.0001 # convert to ha
# 
# 
# 
# 
# 
# ggplot()+
#   geom_sf(data = shp, fill = 'red')+
#   geom_sf(data = shp_buff, fill = 'transparent', colour = 'blue')
# 
# 
# 
# hucs_conv = lapply(hucs, function(huc, r, shp, crops){
#   
#   # retain hucs that overlap with shp
#   huc_shp = huc[apply(st_intersects(huc, st_transform(shp, st_crs(huc)), sparse = F), 1, any),]
#   
#   # count all conversion pixels (value != 0) in each huc; sum 
#   huc_all_conv = sum(exact_extract(
#     x = r,
#     y = st_transform(huc10_shp, st_crs(r)),
#     fun = function(value, cov_frac) sum(ifelse(value != 0,1,0)*cov_frac, na.rm = T)
#   ), na.rm = T)
#   
#   # iterate through each crop 
#   huc_crop_conv = lapply(crops, function(crop, r, huc_shp){
#     
#     # count all conversion pixels of given crop in huc; sum
#     huc_crop = sum(exact_extract(
#       x = r,
#       y = st_transform(huc10_shp, st_crs(r)),
#       fun = function(value, cov_frac) sum(ifelse(value == crop,1,0)*cov_frac, na.rm = T)
#     ), na.rm = T)
# 
#     return(huc_crop)
#     
#   }, r = r, huc_shp = huc_shp)
#   
#   # name elements of huc_crop_conv output
#   names(huc_crop_conv) = crops
#   huc_crop_conv = unlist(huc_crop_conv)
#   
#   # calculate area of huc polygons (units taken from raster crs)
#   huc_area = sum(st_area(st_transform(huc10_shp, st_crs(r))), na.rm = T) # units are those of r: m
#   
#   # prepare output: convert conversion area to meters; then everthing to hectares
#   out = c('all' = huc_all_conv, huc_crop_conv)*prod(res(r)) # convert to m2
#   out = c(out, 'area' = huc_area)*0.0001 # convert to ha
#   
#   return(out)
#   
# }, r = r, shp = shp, crops = crops)
# 
# return(unlist(hucs_conv))
# 
# 
# 
# 
# 
# huc_conversion_area = function(shp, huc10, huc12, crops){
#   
#   # retain hucs that intersect with the habitat polygons
#   huc10_shp = huc10[apply(st_intersects(huc10, st_transform(shp, st_crs(huc10)), sparse = F), 1, any),]
#   huc12_shp = huc12[apply(st_intersects(huc12, st_transform(shp, st_crs(huc12)), sparse = F), 1, any),]
#   
#   # calculating the number of conversion pixels (value != 0) in each huc 
#   huc10_conv = exact_extract(
#     x = r,
#     y = st_transform(huc10_shp, st_crs(r)),
#     fun = function(value, cov_frac) sum(ifelse(value != 0,1,0)*cov_frac, na.rm = T)
#   )
#   
#   huc12_conv = exact_extract(
#     x = r,
#     y = st_transform(huc12_shp, st_crs(r)),
#     fun = function(value, cov_frac) sum(ifelse(value != 0,1,0)*cov_frac, na.rm = T)
#   )
#   
#   # get sum of conversion in each huc category
#   total_conv_hucs = c('all.huc10' = sum(huc10_conv), 'all.huc12' = sum(huc12_conv))
#   
#   
#   
#   # iterate through each crop for which conversion totals should be tabulated
#   crop_conv_hucs = lapply(crops, function(crop){
#     
#     # calculating the number of pixels of crop pixels in huc10 and huc12 basincs
#     huc10_conv = exact_extract(
#       x = r,
#       y = st_transform(huc10_shp, st_crs(r)),
#       fun = function(value, cov_frac) sum(ifelse(value == crop,1,0)*cov_frac, na.rm = T)
#     )
#     
#     huc12_conv = exact_extract(
#       x = r,
#       y = st_transform(huc12_shp, st_crs(r)),
#       fun = function(value, cov_frac) sum(ifelse(value == crop,1,0)*cov_frac, na.rm = T)
#     )
#     
#     # return the sum of conversion in each huc category
#     return(c('huc10' = sum(huc10_conv), 'huc12' = sum(huc12_conv)))
#     
#   })
#   
#   # name list elements
#   names(crop_conv_hucs) = crops
#   
#   # unlist to vector (fame format becomes: [crop].huc[XX])
#   crop_conv_hucs = unlist(crop_conv_hucs)
#   
#   
#   # calculate area of huc polygons (units taken from raster crs)
#   huc10_area = st_area(st_transform(huc10_shp, st_crs(r))) # units are those of r: m
#   huc12_area = st_area(st_transform(huc12_shp, st_crs(r))) # units are those of r: m
#   
#   huc_areas = c('area.huc10' = sum(huc10_area), 'area.huc12' = sum(huc12_area))
#   
#   
#   # prepare output: convert conversion area to meters; then everthing to hectares
#   out = c(total_conv_hucs, crop_conv_hucs)*prod(res(r)) # convert to m2
#   out = c(out, huc_areas)*0.0001 # convert to ha
#   
#   return(out)
#   
# }
# 
# 
# 
# 
# 
# 
# 
# 
# # calculating the number of pixels  of conversion (mtr = 3) in each county
# test = exact_extract(
#   x = r,
#   y = st_transform(huc10_shp, st_crs(r)),
#   fun = function(value, cov_frac) sum(ifelse(value == cropVal,1,0)*cov_frac, na.rm = T)
# )
# 
# test = exact_extract(
#   x = r,
#   y = st_transform(huc12_shp, st_crs(r)),
#   fun = function(value, cov_frac) sum(ifelse(value == cropVal,1,0)*cov_frac, na.rm = T)
# )
# 
# 
# ggplot()+
#   geom_sf(data = shp, fill = 'red')+
#   geom_sf(data = huc10_shp, fill = 'transparent', colour = 'blue')+
#   geom_sf(data = huc12_shp, fill = 'transparent', colour = 'blue')
# 
# 
# i = 1
# # iterate through each polygon in shp
# #for(i in 1:nrow(shp)){
#   
#   cat('Polygon ', i, ' of ', nrow(shp), "\n")
#   
#   # get individual polygon
#   poly = shp[i,]
#   
#   # fill holes in poly to prevent topology issues
#   poly = fill_holes(poly, threshold = st_area(poly))
#   
#   # retain hucs that intersect with the habitat polygons
#   huc10_poly = huc10[apply(st_intersects(huc10, st_transform(poly, st_crs(huc10)), sparse = F), 1, any),]
#   huc12_poly = huc12[apply(st_intersects(huc12, st_transform(poly, st_crs(huc12)), sparse = F), 1, any),]
#   
#   
#   
#   
#   
# 
# require(ggplot2)
# ggplot()+
#   geom_sf(data = shp, fill = 'red')
# 
# 
# # dissolve adjoining polygons
# shp = st_as_sf(st_cast(st_union(shp), "POLYGON"))
# 
# # function(shp, huc8, huc10)
# 
# packages = c(
#   "sf",
#   "raster",
#   "smoothr",
#   "exactextractr"
# )
# 
# install.packages(setdiff(packages, rownames(installed.packages()))) 
# lapply(packages, require, character.only = TRUE)
# 
# # fill holes in poly to prevent topology issues
# poly = fill_holes(poly, threshold = st_area(poly))
# 
# # retain hucs that intersect with the habitat polygons
# huc10_poly = huc10[apply(st_intersects(huc10, st_transform(poly, st_crs(huc10)), sparse = F), 1, any),]
# huc12_poly = huc12[apply(st_intersects(huc12, st_transform(poly, st_crs(huc12)), sparse = F), 1, any),]
# 
# ggplot()+
#   geom_sf(data = poly, fill = 'red')+
#   geom_sf(data = poly_hucs, fill = 'transparent', colour = 'blue')
# 
