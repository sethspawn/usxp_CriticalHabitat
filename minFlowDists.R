
# need a seperate condition for line string features
# calculate total area within polygon and polygon area
# calculate within search drainage area of upstream flow dist
# calculate area in hucs


rm(list = ls())

#====================================================================================================
# load required packages, installing any that have not yet been installed

packages = c(
  "boxr",
  "raster",
  "sf",
  "units",
  "nhdplusTools",
  # "concaveman",
  # "plyr",
  "stplanr",
  "igraph",
  "smoothr",
  "exactextractr",
  "ggplot2"
)

install.packages(setdiff(packages, rownames(installed.packages())))
lapply(packages, require, character.only = TRUE)

source('./findDrainageOutlets.R')
source('./conversionPatches.R')

#====================================================================================================

# Authenticate BOX connection using Client ID and Secret in renviron

Sys.setenv(BOX_CLIENT_ID = "bsxau952md09y1t6t0s6jyvx1ri0f68n")
Sys.setenv(BOX_CLIENT_SECRET = "OjKUIBQnsCR2iRF9JSVwz37BIUrvRmcH")

box_auth()

#====================================================================================================

defineSearchBasin = function(poly, comids, search_dist_km, crs){
  
  # get the merged extent of each comids drainage basin
  for(i in 1:length(comids)){
    
    outlet_comid = list(featureSource = "comid", featureID = comids[i])
    
    basin = get_nldi_basin(nldi_feature = outlet_comid)
    
    if(i == 1){
      
      search_basins = basin
      
    }else{
      
      search_basins = st_union(search_basins, basin)
      
    }
    
  }
  
  # reproject to match raster crs
  search_basins = st_transform(search_basins, st_crs(crs))
  
  # fill holes in poly to prevent topology issues
  poly = fill_holes(poly, threshold = st_area(poly))
  
  # # reproject poly to match raster crs
  poly = st_transform(poly, st_crs(crs))
  
  # buffer poly to search distance
  poly_buff = st_buffer(poly, set_units(search_dist_km, km))
  
  # clip search basin to no greater than search distance from polygon
  search_basins_clip = st_intersection(search_basins, poly_buff)
  
  return(search_basins_clip)
  
}


#----------------------------------------------------------------------------------------------------

minFlowDist = function(poly, comid, patch_pts, upstream_search_dist_km){
  
  if(is.null(patch_pts)){
    
    annual_min_dists = NA
    
  }else{
    
    # fill holes in poly to prevent topology issues
    poly = fill_holes(poly, threshold = st_area(poly))
    
    # set max network distance such that the min max distance (at top of reach) = search_dist
    distance_km = round(as.numeric((st_length(st_cast(poly, 'LINESTRING'))/2/1000)), digits = 0) + upstream_search_dist_km
    distance_km = set_units(ifelse(distance_km > 9999, 9999, distance_km), km)
    
    # get flowlines from comid
    outlet_comid = list(featureSource = "comid", featureID = comid)
    flowlines = navigate_nldi(outlet_comid, mode = "upstreamTributaries", distance_km = distance_km)
    
    # get drainage basin of comid outlet point
    dbasin = get_nldi_basin(nldi_feature = outlet_comid)
    
    # reproject all data to match flowlines (wgs84)
    patch_pts = st_transform(patch_pts, crs = st_crs(flowlines))
    dbasin = st_transform(dbasin, crs = st_crs(flowlines))
    poly = st_transform(poly, st_crs(flowlines))
    
    # get point sf of drainage outlet (comid)
    outlet_pt = st_intersection(poly, st_cast(flowlines[flowlines$nhdplus_comid == comid,], 'LINESTRING'))
    
    # only retain patch_pts within the drainage basin of comid
    patch_pts = patch_pts[st_within(patch_pts, dbasin, sparse = F),]
    
    # if no patch points within drainage basin, return NA
    if(nrow(patch_pts) == 0){
      
      annual_min_dists = NA
      
    }else{
      
      # get points at which flowlines enter critical habitat:
      ch_inlet_pts = st_cast(st_intersection(poly, flowlines), 'POINT')
      
      # add segment length attribute to flowlines (this will be network weight)
      flowlines$length = st_length(flowlines) # units: meters
      
      # convert flowlines to a network
      sln = SpatialLinesNetwork(flowlines)
      
      # identify network nodes closest to the end (habitat inlets) points
      end_nodes = na.omit(unique(find_network_nodes(sln, x = st_coordinates(ch_inlet_pts))))
      
      # create empty comtainers in which to accumulate minimum dists for each year
      #comid_annual_min_dists = c()
      all_annual_min_dists = list()
      
      # iterate through each year of conversion
      for(year in na.omit(unique(patch_pts$year))){
        
        # subset patch_pts by year
        sub_patch_pts = patch_pts[patch_pts$year == year,]
        
        # identify network nodes closest to the start (conversion patches) points
        #start_nodes = na.omit(unique(find_network_nodes(sln = sln, x = st_coordinates(sub_patch_pts), maxdist = 10000)))
        start_nodes = find_network_nodes(sln = sln, x = st_coordinates(sub_patch_pts), maxdist = 10000)
        
        # calculate exhaustive distance matrix from start to end nodes
        dist = distances(sln@g, v = start_nodes, to = end_nodes, weight = NULL) # NULL defers to the networks weight attr = 'length'
        
        # calculate the min distance to CH for each node; add as attribute to sub_patch_pts
        sub_patch_pts[,'min_dist'] = do.call(pmin, as.data.frame(dist))
        
        
        all_annual_min_dists[[as.character(year)]] = sub_patch_pts
      
        
        # # identify minimum dist. If it's greater than the search distance (in m), set to NA
        # min_dist = ifelse(set_units(min(dist)/1000, km) > set_units(upstream_search_dist_km, km), NA, set_units(min(dist)/1000, km))
        # 
        # # save min_dist to annual min_dists_container
        # comid_annual_min_dists = c(comid_annual_min_dists, min_dist)
        
      }
      
      # # name the numeric
      # names(comid_annual_min_dists) = paste0('minDist_year', na.omit(unique(patch_pts$year)))
      
      # get all possible unique column names in outList
      all_names = unique(unlist(lapply(all_annual_min_dists, function(e) names(e))))
      
      # combine all list elements by matching column names
      all_annual_min_dists = do.call(rbind, lapply(all_annual_min_dists, function(x) x[match(all_names, names(x))]))
      
      all_annual_min_dists[all_annual_min_dists$min_dist <= upstream_search_dist_km*1000,]
      
    }
    
  }
  
  # return(annual_min_dists)
  return(all_annual_min_dists)
  
}

#====================================================================================================
# remove summarys just return all min dists as a sf data frame  ***
# calculate the summary at the end min by year for each polygon and total area by per year
# do you use the river distance for the threshold of the area calculation or the buffer distance? Some how need to determine
# watershed area in which conversion area is tabulated.

#----------------------------------------------------------------------------------------------------

r = r_ytc


annualMinFlowDists = function(poly, comids, patch_pts, search_dist_km){
  
  # fill holes in poly to prevent topology issues
  poly = fill_holes(poly, threshold = st_area(poly))
  
  # set max network distance such that the min max distance (at top of reach) = ~ search_dist
  distance_km = round(as.numeric((st_length(st_cast(poly, 'LINESTRING'))/2/1000)), digits = 0) + search_dist_km
  distance_km = set_units(ifelse(distance_km > 9999, 9999, distance_km), km)
  
  # create empty list in which to accumulate each comid's annual min dists
  outList = list()
  
  # iterate through each comid
  for(comid in comids){
    
    # get flowlines from comid
    outlet_comid = list(featureSource = "comid", featureID = comid)
    flowlines = navigate_nldi(outlet_comid, mode = "upstreamTributaries", distance_km = distance_km)
    
    # get drainage basin of comid outlet point
    dbasin = get_nldi_basin(nldi_feature = outlet_comid)
    
    # reproject all data to match flowlines (wgs84)
    patch_pts = st_transform(patch_pts, crs = st_crs(flowlines))
    dbasin = st_transform(dbasin, crs = st_crs(flowlines))
    poly = st_transform(poly, st_crs(flowlines))
    
    # get point sf of drainage outlet (comid)
    outlet_pt = st_intersection(poly, st_cast(flowlines[flowlines$nhdplus_comid == comid,], 'LINESTRING'))
    
    # only retain patch_pts within the drainage basin of comid
    patch_pts_in_db = patch_pts[st_within(patch_pts, dbasin, sparse = F),]
    
    # if no patch points within drainage basin, return NA
    if(nrow(patch_pts_in_db ) == 0){
      
      next
      
    }else{
      
      # get points at which flowlines enter critical habitat:
      ch_inlet_pts = st_cast(st_intersection(poly, flowlines), 'POINT')
      
      # add segment length attribute to flowlines (this will be network weight)
      flowlines$length = st_length(flowlines) # units: meters
      
      # convert flowlines to a network
      sln = SpatialLinesNetwork(flowlines)
      
      # identify network nodes closest to the end (habitat inlets) points
      end_nodes = na.omit(unique(find_network_nodes(sln, x = st_coordinates(ch_inlet_pts))))
      
      # identify network nodes closest to the start (conversion patches) points   <<========================= CAN WE CALCULATE HOW FAR EACH PT WAS FROM NODE?
      start_nodes = find_network_nodes(sln = sln, x = st_coordinates(patch_pts_in_db), maxdist = 10000)
      
      # create key to address NAs
      start_node_key = data.table(clump_id2 = patch_pts_in_db$clump_id2, start_nodes = start_nodes)
      start_node_key = start_node_key[complete.cases(start_node_key),]
      
      # calculate exhaustive distance matrix from start to end nodes
      dist = distances(sln@g, v = start_node_key$start_nodes, to = end_nodes, weight = NULL) # NULL defers to the networks weight attr = 'length'
      
      # calculate the min distance to CH for each node; add as attribute to sub_patch_pts
      start_node_key[,'min_dist'] = do.call(pmin, as.data.frame(dist))
      start_node_key[,start_nodes:=NULL] # remove start_node column
      
      # add comid's pts to outList
      outList[[as.character(comid)]] = merge(patch_pts_in_db, start_node_key, by = "clump_id2")
      
    }
    
  }
  
  # get all possible unique column names in outList
  all_names = unique(unlist(lapply(outList, function(e) names(e))))
  
  # combine all list elements by matching column names
  annual_min_dists = do.call(rbind, lapply(outList, function(x) x[match(all_names, names(x))]))
  
  return(annual_min_dists)
  
}
  
  
#   
#   
#     
#     if(is.null(patch_pts)){
#       
#       annual_min_dists = NA
#       
#     }else{
#       
#       # fill holes in poly to prevent topology issues
#       poly = fill_holes(poly, threshold = st_area(poly))
#       
#       # set max network distance such that the min max distance (at top of reach) = search_dist
#       distance_km = round(as.numeric((st_length(st_cast(poly, 'LINESTRING'))/2/1000)), digits = 0) + upstream_search_dist_km
#       distance_km = set_units(ifelse(distance_km > 9999, 9999, distance_km), km)
#       
#       # get flowlines from comid
#       outlet_comid = list(featureSource = "comid", featureID = comid)
#       flowlines = navigate_nldi(outlet_comid, mode = "upstreamTributaries", distance_km = distance_km)
#       
#       # get drainage basin of comid outlet point
#       dbasin = get_nldi_basin(nldi_feature = outlet_comid)
#       
#       # reproject all data to match flowlines (wgs84)
#       patch_pts = st_transform(patch_pts, crs = st_crs(flowlines))
#       dbasin = st_transform(dbasin, crs = st_crs(flowlines))
#       poly = st_transform(poly, st_crs(flowlines))
#       
#       # get point sf of drainage outlet (comid)
#       outlet_pt = st_intersection(poly, st_cast(flowlines[flowlines$nhdplus_comid == comid,], 'LINESTRING'))
#       
#       # only retain patch_pts within the drainage basin of comid
#       patch_pts = patch_pts[st_within(patch_pts, dbasin, sparse = F),]
#       
#       # if no patch points within drainage basin, return NA
#       if(nrow(patch_pts) == 0){
#         
#         annual_min_dists = NA
#         
#       }else{
#         
#         # get points at which flowlines enter critical habitat:
#         ch_inlet_pts = st_cast(st_intersection(poly, flowlines), 'POINT')
#         
#         # add segment length attribute to flowlines (this will be network weight)
#         flowlines$length = st_length(flowlines) # units: meters
#         
#         # convert flowlines to a network
#         sln = SpatialLinesNetwork(flowlines)
#         
#         # identify network nodes closest to the end (habitat inlets) points
#         end_nodes = na.omit(unique(find_network_nodes(sln, x = st_coordinates(ch_inlet_pts))))
#         
#         # create empty comtainers in which to accumulate minimum dists for each year
#         #comid_annual_min_dists = c()
#         all_annual_min_dists = list()
#         
#         # iterate through each year of conversion
#         for(year in na.omit(unique(patch_pts$year))){
#           
#           # subset patch_pts by year
#           sub_patch_pts = patch_pts[patch_pts$year == year,]
#           
#           # identify network nodes closest to the start (conversion patches) points
#           #start_nodes = na.omit(unique(find_network_nodes(sln = sln, x = st_coordinates(sub_patch_pts), maxdist = 10000)))
#           start_nodes = find_network_nodes(sln = sln, x = st_coordinates(sub_patch_pts), maxdist = 10000)
#           
#           # calculate exhaustive distance matrix from start to end nodes
#           dist = distances(sln@g, v = start_nodes, to = end_nodes, weight = NULL) # NULL defers to the networks weight attr = 'length'
#           
#           # calculate the min distance to CH for each node; add as attribute to sub_patch_pts
#           sub_patch_pts[,'min_dist'] = do.call(pmin, as.data.frame(dist))
#           
#           
#           all_annual_min_dists[[as.character(year)]] = sub_patch_pts
#           
#           
#           # # identify minimum dist. If it's greater than the search distance (in m), set to NA
#           # min_dist = ifelse(set_units(min(dist)/1000, km) > set_units(upstream_search_dist_km, km), NA, set_units(min(dist)/1000, km))
#           # 
#           # # save min_dist to annual min_dists_container
#           # comid_annual_min_dists = c(comid_annual_min_dists, min_dist)
#           
#         }
#         
#         # # name the numeric
#         # names(comid_annual_min_dists) = paste0('minDist_year', na.omit(unique(patch_pts$year)))
#         
#         # get all possible unique column names in outList
#         all_names = unique(unlist(lapply(all_annual_min_dists, function(e) names(e))))
#         
#         # combine all list elements by matching column names
#         all_annual_min_dists = do.call(rbind, lapply(all_annual_min_dists, function(x) x[match(all_names, names(x))]))
#         
#         all_annual_min_dists[all_annual_min_dists$min_dist <= upstream_search_dist_km*1000,]
#         
#       }
#       
#     }
#     
#     # return(annual_min_dists)
#     return(all_annual_min_dists)
#     
#     # # if search distance is set to 0, skip this and the remaining iterations
#     # if(search_dist_km == 0){   # <<================ THIS ISN'T RELEVANT IF THE DOWNSIZING CODE IS REMOVED.  ***
#     #   
#     #   next
#     #   
#     # }else{
#     #   
#     #   # calculate the annual min dists for comid
#     #   min_dists = minFlowDist(poly, comid = comid, patch_pts = patch_pts, upstream_search_dist_km = search_dist_km)
#     #   
#     #   # add to outList
#     #   outList[[as.character(comid)]] = min_dists
#       
#       # # if there are no NAs in min dists, all years are conderdered, and the maximum min dist is less than search dist...
#       # if(!any(is.na(min_dists)) & length(min_dists) == length(unique(patch_pts$year)) & max(min_dists) < search_dist_km){
#       #   
#       #   # ...reset search dist to max min dist
#       #   search_dist_km = round(max(min_dists) + 1, digits = 0)
#       #   
#       #   cat('COMID: ', comid, '. New search dist = ', search_dist_km, '\n')
#       #   
#       # }
#       
#     }
#     
#   }
#   
#   # get all possible unique column names in outList
#   all_names = unique(unlist(lapply(outList, function(e) names(e))))
#   
#   # combine all list elements by matching column names
#   min_dists = do.call(rbind, lapply(outList, function(x) x[match(all_names, names(x))]))
#   colnames(min_dists) = all_names
#   
#   # calculate the mimimum distance in each dataframe column (year)
#   min_dists = apply(min_dists, 2, function(x) min(x, na.rm = T))
#   
#   # if min_dist is infinate, change to NA
#   min_dists = replace(min_dists, is.infinite(min_dists), NA)
#   
#   return(min_dists)
#   
# }


#----------------------------------------------------------------------------------------------------


flowDist2CH = function(shp, r, search_dist_km, allYears){
  
#### IF EXACT EXTRACTER RETURNS THAT ALL YEARS ARE REPRESENTED BY CONVERSION WITHIN POLYGON, CALCULATE AREAS

#### ELSE...  
   
  # create list in which to accumulate results for conversion within and outside of each polygon
  out_poly_list = list()
  in_poly_list  = list()
  
  # iterate through each polygon in shp
  for(i in 1:nrow(shp)){
    
    cat('Polygon ', i, ' of ', nrow(shp), "\n")
    
    # get individual polygon
    poly = shp[i,]
    
    # determinE outlet comids
    comids = findDrainageOutlets(poly = poly)
    
    # get maximum search_basin
    search_basin = defineSearchBasin(poly = poly, comids = comids, search_dist_km = search_dist_km, crs = crs(r))
    
    # get patch points within search_basin
    patch_pts = conversionPatches(r = r, region = search_basin, years = allYears)
    
    # if no conversion patches within the search basin, skip to next polygon
    if(nrow(patch_pts) == 0){
      
      next
      
    } else{
      
      #### INSERT: FILTER OUT PATCHES WITHIN POLYGON 
      patch_pts_out_poly = patch_pts[!st_within(patch_pts, st_transform(poly, st_crs(patch_pts)), sparse = F),]
      
      # calculate minimum distances between all conversion patches outside of polygon and nearest CH
      out_poly_list[[as.character(i)]] = annualMinFlowDists(poly = poly, comids = comids, patch_pts = patch_pts_out_poly, search_dist_km = search_dist_km)
      
      in_poly_list[[as.character(i)]]  = patch_pts[st_within(patch_pts, st_transform(poly, st_crs(patch_pts)), sparse = F),]
      
    }
    
  }
  
  # get all possible unique column names in out_poly_list
  all_names_out_poly = unique(unlist(lapply(out_poly_list, function(e) names(e))))
  
  # combine all list elements by matching column names
  out_poly_annual_min_dists = do.call(rbind, lapply(out_poly_list, function(x) x[match(all_names_out_poly, names(x))]))
  
  # remove those patches whose flow distance is greater than the search distance
  out_poly_annual_min_dists = out_poly_annual_min_dists[out_poly_annual_min_dists $min_dist <= search_dist_km*1000,]
  
  # get all possible unique column names in in_poly_list
  all_names_in_poly = unique(unlist(lapply(in_poly_list, function(e) names(e))))
  
  # combine all list elements by matching column names
  in_poly_annual_min_dists = do.call(rbind, lapply(in_poly_list, function(x) x[match(all_names_in_poly, names(x))]))
  
  return(list('in' = in_poly_annual_min_dists, 'out' = out_poly_annual_min_dists))
  
}


#====================================================================================================


if(!file.exists('s35_ytc.tif')){
  
  # download s35 ytc raster
  box_dl(file_id = '754823568329', overwrite = T, file_name = 's35_ytc.tif')
  
}

# load the usxp s35 ytc raster (year of each patches conversion to cropland)
r_ytc = raster("s35_ytc.tif")

# load critical habitat polys
file = 'FCH_Notropis_topeka_20040727.zip'
#file = 'FCH_Gila_cypha_19940321.zip'

# unzip a shapefile to 'temp_unzip' subdirectory
unzip(file.path('crithab_all_shapefiles', file), exdir = file.path('crithab_all_shapefiles', 'unzipped'))

# read unzipped shapefile
shp = try(st_read(file.path('crithab_all_shapefiles', 'unzipped', gsub('.zip', '.shp', file))))

# dissolve adjoining polygons
shp = st_as_sf(st_cast(st_union(shp), "POLYGON"))

test = flowDist2CH(shp = shp, r = r_ytc, search_dist_km = 100, allYears = 2009:2016)  ## would be faster if we started by searching the smallest search basins





# 
# 
# 
# poly = shp[1,]
# comid = comids[1]
# 
# ggplot()+geom_sf(data = poly)
# 
# getDrainageBasin = function(poly, comids, upstream_dist_km, crs){
#   
#   
#   
#   # set max network distance such that the min max distance (at top of reach) = search_dist
#   upstream_dist_km = round(sum(as.numeric(st_length(st_cast(poly, 'LINESTRING'))))/1000, digits = 0)
#   upstream_dist_km = set_units(ifelse(distance_km > 9999, 9999, distance_km), km)
#   
#   # get flowlines from comid
#   outlet_comid = list(featureSource = "comid", featureID = comid)
#   flowlines = navigate_nldi(outlet_comid, mode = "upstreamTributaries", distance_km = distance_km)
#   
#   # get the merged extent of each comids drainage basin
#   for(i in 1:length(comids)){
#     
#     outlet_comid = list(featureSource = "comid", featureID = comids[i])
#     
#     basin = get_nldi_basin(nldi_feature = outlet_comid)
#     
#     if(i == 1){
#       
#       search_basins = basin
#       
#     }else{
#       
#       search_basins = st_union(search_basins, basin)
#       
#     }
#     
#   }
#   
#   # reproject to match raster crs
#   search_basins = st_transform(search_basins, st_crs(crs))
#   
#   # fill holes in poly to prevent topology issues
#   poly = fill_holes(poly, threshold = st_area(poly))
#   
#   # # reproject poly to match raster crs
#   poly = st_transform(poly, st_crs(crs))
#   
#   # buffer poly to search distance
#   poly_buff = st_buffer(poly, set_units(search_dist_km, km))
#   
#   # clip search basin to no greater than search distance from polygon
#   search_basins_clip = st_intersection(search_basins, poly_buff)
#   
#   return(search_basins_clip)
#   
# }