#rm(list = ls())

#====================================================================================================
# load required packages, installing any that have not yet been installed

packages = c(
  "sf",
  "nhdR",
  "nhdplusTools",
  "smoothr"
)

install.packages(setdiff(packages, rownames(installed.packages()))) 
lapply(packages, require, character.only = TRUE)

#================================================================================================================


findDrainageOutlets = function(poly){
  
  # explicitly cast poly as a polygon
  poly = st_cast(poly, "POLYGON")
  
  # fill any internal holes in poly (throw errors when calculating intersection)
  poly = fill_holes(poly, threshold = st_area(poly))
  
  # load vpus for querying nhd
  vpus = nhdR::vpu_shp
  
  # get vpus that overlap with poly
  poly_vpus = vpus[apply(st_intersects(vpus, st_transform(poly, st_crs(vpus)), sparse = F), 1, any),]
  poly_vpus = as.character(poly_vpus$UnitID[poly_vpus$UnitType == 'VPU'])
  
  # iterate through each vpu
  outlet_comids = lapply(poly_vpus, function(vpu){
    
    # download the nhd for the given vpu
    nhd_plus_get(vpu = vpu, "NHDSnapshot")
    
    # load vpu's flowlines
    vpu_flowlines = nhd_plus_load(vpu = vpu, "NHDSnapshot", "NHDFlowline")
    
    #vpu_flowlines_clipped = st_intersection(vpu_flowlines, poly)
    
    # get points at which habitat polygon intersects the flowlines.
    pts = st_cast(st_intersection(st_cast(poly, 'LINESTRING'), vpu_flowlines), 'POINT')
    
    # capitalization is inconsistent across VPUs - standardize native colnames as captitolized
    colnames(pts) = toupper(colnames(pts))
    
    # create columns in which to report the outlet comid and the number of upstream comids
    pts['outlet_comid'] = rep(NA, nrow(pts))
    pts['upstrm_comid_count'] = rep(0, nrow(pts))
    
    # iterate through each comid id
    for(comid in pts$COMID){
      
      # calculate a search distance (here, based on the perimete of the polygon)
      distance_km = round(as.numeric((st_length(st_cast(poly, 'LINESTRING'))/1000)), digits = 0) # polygon's perimeter length
      
      # get upstream flowlines of comid
      flowlines = navigate_nldi(list(featureSource = "comid", 
                                     featureID = comid), 
                                mode = "upstreamTributaries", 
                                distance_km =  distance_km )
      
      # identify upstream comids
      upstream_comids = pts$COMID[which(pts$COMID %in% flowlines$nhdplus_comid)]
      
      # if points are upstream and their prior upstream comid count is less than current
      pts$outlet_comid[which(pts$COMID %in% upstream_comids & pts$upstrm_comid_count < length(upstream_comids))] = comid
      
      # set upstream comid count for all selected comids (current and upstream) to current comid count
      pts$upstrm_comid_count[which(pts$COMID %in% upstream_comids & pts$upstrm_comid_count < length(upstream_comids))] = length(upstream_comids)
      
    }
    
    # if outlet_comid not found, assign outlet_comid as self.
    pts$outlet_comid[which(is.na(pts$outlet_comid))] = pts$COMID[which(is.na(pts$outlet_comid))]
    
    return(unique(pts$outlet_comid))
    
  })
  
  return(unlist(outlet_comids))
  
}

#================================================================================================================

# # load critical habitat polys
# file = 'FCH_Notropis_topeka_20040727.zip'
# 
# # unzip a shapefile to 'temp_unzip' subdirectory
# unzip(file.path('crithab_all_shapefiles', file), exdir = file.path('crithab_all_shapefiles', 'unzipped'))
# 
# # read unzipped shapefile
# shp = try(st_read(file.path('crithab_all_shapefiles', 'unzipped', gsub('.zip', '.shp', file))))
# 
# # dissolve adjoining polygons
# shp = st_as_sf(st_cast(st_union(shp), "POLYGON"))
# 
# outlets = lapply(split(shp, 1:nrow(shp)), findDrainageOutlets)

#================================================================================================================

# require(ggplot2)
# ggplot()+
#   #geom_sf(data = poly_vpus)+
#   geom_sf(data = poly, fill = 'red')+
#   geom_sf(data = vpu_flowlines_clipped)+
#   geom_sf(data = pts, aes(colour = outlet_comid))+
#   geom_sf(data = flowlines, color = 'green')+
#   geom_sf(data = pts[which(is.na(pts$outlet_comid)),], colour = 'green', size = 3)

