
#rm(list = ls())

#====================================================================================================
# load required packages, installing any that have not yet been installed

packages = c(
  "raster",
  "sf",
  "plyr",
)

install.packages(setdiff(packages, rownames(installed.packages())))
lapply(packages, require, character.only = TRUE)

source('./findDrainageOutlets.R')

#----------------------------------------------------------------------------------------------------


conversionPatches = function(r, region, years = NULL){
  
  # confirm projections of r and region match
  if(crs(r) != st_crs(region)){
    
    region = st_transform(region, st_crs(r))
    
  }
  
  # clip r to the clipped search baisn extent
  r_clip = crop(r, region)
  
  # calculate patch ids for conversion patches within watershed
  # from: https://stackoverflow.com/questions/15632630
  r_ids = clump(r_clip,  directions = 4) # using 4 directions to seperate fields that might be on a grid
  clump_id = getValues(r_ids)
  xy = xyFromCell(r_ids,1:ncell(r_ids))
  df = data.frame(xy, clump_id, is_clump = r_ids[] %in% freq(r_ids, useNA = 'no')[,1])
  df = df[df$is_clump == T, ]
  
  if(nrow(df) > 0){
    
    # add corresponding year as attribute
    df['year'] = extract(r_clip, df[,c('x', 'y')])
    df = df[complete.cases(df),]
    
    # some clumps include conversion in more than one year --> make truly unique id:
    df['clump_id2'] = paste0(df$clump_id, '_', df$year)
    
    # get centroid coordinates and year of each patch       <<<======= could be faster
    dfm = ddply(df, .(clump_id2), summarise, xm = mean(x), ym = mean(y), year = mean(year))
    
    # convert to sf points
    patch_pts = st_as_sf(dfm, coords = c('xm', 'ym'), crs = crs(r))
    
    # filter points so only those within regions_clip are retained
    patch_pts = patch_pts[st_within(patch_pts, region, sparse = F),]
    
    # if years are provided, filter them
    if(!is.null(years)){
      
      # filter points to only include those in "years"
      patch_pts = patch_pts[which(patch_pts$year %in% years),]
      
    }
    
  }else{
    
    patch_pts = NULL
    
  }
  
  return(patch_pts)
  
}

