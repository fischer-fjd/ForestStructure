# classify forest example
library(terra)
library(sf)
library(data.table)

dir_current = dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
dir_current
setwd(dir_current)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### 1. Function definition  ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

classify.forest = function(file_chm, path_output = "", path_output_cover = "", minheight_canopy = 2, mincover_canopy = 0.1, minarea_forest = 5000, factor_refine = 0.0, buffer_scan = 50, minarea_patch = 0, maxarea_hole = 1000000, ratio_patch_to_hole = 5, factor_preagg = 1, plot_layer = F){
  
  # define a minimum area for the analysis: this is the minimum size of patches we want to analysis and may be relevant for metrics that only make sense at large spatial scales (gap size frequency distributions) or when sampling large continuous areas
  # this is independent of the area for forest definition, but defaults to it when not provided
  if(minarea_patch < minarea_forest) minarea_patch = minarea_forest
  
  # read in chm and delineate boundaries
  chm = rast(file_chm)
  chm_background = clamp(chm, upper = 1, lower = 1, values = T)
  
  if(factor_preagg > 1){
    chm_background = aggregate(chm_background, fact = factor_preagg, fun = "mean", na.rm = T)
  }
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  cat("Delineate study area\n")
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  
  polygons_background = as.polygons(chm_background, values = F)
  polygons_background = disagg(polygons_background)
  polygons_background = fillHoles(polygons_background)
  polygons_background = st_as_sf(polygons_background)
  polygons_background = st_simplify(polygons_background, preserveTopology = T, dTolerance = buffer_scan * 0.25)
  polygons_background = st_buffer(polygons_background, -buffer_scan)
  polygons_background = polygons_background[!st_is_empty(polygons_background),]
  polygons_background = vect(polygons_background)
  polygons_background$isforest = 0
  
  # the trim function in terra does not seem to work at the moment
  # so do trimming manually
  # in future: chm = trim(mask(chm, polygons_background))
  extent_background = ext(polygons_background)
  chm = crop(chm, extent_background)
  chm = mask(chm, polygons_background)
  
  # get total chm area and 
  chm_area = length(cells(chm))
  if(minarea_forest  > 0.5 * chm_area) minarea_forest  = as.integer(0.5 * chm_area)
  
  #%%%%%%%%%%%%%%%%%%%%%%%#
  cat("Get canopy cover\n")
  #%%%%%%%%%%%%%%%%%%%%%%%#
  
  iscanopy = ifel(chm < minheight_canopy, 0, 1) # canopy / non-canopy
  
  if(factor_preagg > 1){
    iscanopy = aggregate(iscanopy, fact = factor_preagg, fun = "mean", na.rm = T)
  }
  
  # convert canopy area into window length
  minradius_forest = sqrt(minarea_forest/pi)
  window_forest = focalMat(iscanopy, d = minradius_forest, type = "circle", fillNA = TRUE)
  window_forest[!is.na(window_forest)] = 1
  
  # compute mean canopy cover
  # issue: focal is very slow, potential speed up by aggregating first
  cover_canopy = focal(iscanopy, w = window_forest, fun = "mean", na.rm = T, na.policy = "omit")
  
  # issue (again) : focal is very slow, potential speed up by aggregating first
  isforest = ifel(cover_canopy >= mincover_canopy,1,0)
  isforest = focal(isforest, w = window_forest, fun = "min", na.rm = T, na.policy = "omit") # we only accept as forest that does not lie within reach of non-forest area
  names(isforest) = "isforest"
  
  if(factor_refine > 0.0){
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    cat("Refining canopy cover threshold adaptively\n")
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    iscanopy_refine = ifel(chm < minheight_canopy, 0, chm/minheight_canopy) 
    if(factor_preagg > 1){
      iscanopy_refine = aggregate(iscanopy_refine, fact = factor_preagg, fun = "mean", na.rm = T)
    } 
    
    cover_canopy_refine = focal(iscanopy_refine, w = window_forest, fun = "mean", na.rm = T, na.policy = "omit")
    cover_canopy_refine = ifel(isforest == 1, cover_canopy_refine,NA)
    
    median_cover = global(cover_canopy_refine, median,na.rm = T)
    mincover_canopy_refine = max(as.numeric(factor_refine * median_cover), mincover_canopy)
    # cat("Overall canopy threshold updated to:",round(mincover_canopy_refine,2),"\n") # not real canopy cover
    
    isforest = ifel(cover_canopy_refine >= mincover_canopy_refine,1,NA)
    isforest = focal(isforest, w = window_forest, fun = "min", na.rm = T, na.policy = "omit") # we only accept as forest that does not lie within reach of non-forest area
    names(isforest) = "isforest"
    isforest[is.na(isforest)] = 0
  }
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  cat("Delineate forest patches\n")
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  
  patches_analysis = as.polygons(isforest)
  patches_analysis = disagg(patches_analysis)
  
  cat("Filter out patches below minimum area and non-forest patches\n")
  patches_analysis$area = expanse(patches_analysis, transform = F)
  patches_analysis = patches_analysis[(patches_analysis$isforest == 1 & patches_analysis$area >= minarea_patch)]
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  cat("Fill holes in patches\n")
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  
  if(nrow(patches_analysis) == 0){
    # case when there are no analysis patches
    patches_analysis = polygons_background
  } else {
    # final fixing of polygons (removing holes, recalculating area)
    # we use smoothr package function, as it allows to set a threshold (thus preventing overlap of forest/noforest polygons)
    patches_analysis = vect(lapply(1:nrow(patches_analysis),function(x){threshold_hole = min(patches_analysis[x]$area/ratio_patch_to_hole,maxarea_hole); y = smoothr::fill_holes(patches_analysis[x], threshold_hole); return(y)}))
    patches_analysis$area = NULL
    
    # remove analysis patches from background and combine again (if they are not coextensive anyways)
    diff_background = erase(polygons_background, patches_analysis)
    patches_analysis_background = disagg(diff_background)
    if(nrow(patches_analysis_background) > 0){
      patches_analysis = rbind(patches_analysis, patches_analysis_background)
    } 
  }
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  cat("Finalize patches and write to file\n")
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  
  patches_analysis = simplifyGeom(patches_analysis, tolerance = 1)
  patches_analysis$area = expanse(patches_analysis, transform = F) 
  
  # adding an ID (starting with forest fragments, then all else)
  patches_analysis$ID = 1:nrow(patches_analysis)
  
  if(path_output == ""){
    path_output = dirname(filename_chm)
  }
  
  writeRaster(cover_canopy, filename = file.path(path_output_cover, "cover_canopy.tif"), overwrite = T)
  writeVector(patches_analysis, filename = file.path(path_output, "patches_analysis.shp"), overwrite = T)
  
  if(plot_layer != F){
    # determine how to do the plotting 
    chm = rast(file_chm)
    ratio_dim_xy = dim(chm)[1]/dim(chm)[2]
    # check if plot should be written to file
    if(plot_layer != T){
      width_pdf = NULL
      height_pdf = NULL
      if(ratio_dim_xy < 1.0){
        height_pdf = 20
        width_pdf = 0.5 * height_pdf / ratio_dim_xy # divide by 2, as we combine two plots
      }else {
        width_pdf = 20
        height_pdf = 0.5 * width_pdf * ratio_dim_xy # divide by 2, as we combine two plots
      }
      pdf(file = plot_layer, width = width_pdf, height = height_pdf)
    } else {
      par_mfrow_previous = par()$mfrow
    }
    if(ratio_dim_xy < 1.0){par(mfrow = c(2,1))} else { par(mfrow = c(1,2))}
    
    plot(ext(chm), main = "Canopy height (m)", mar = c(3,3,3,8))
    plot(chm, col = viridis(100), legend = T, add = T) 
    plot(ext(chm), mar = c(3,3,3,8), main = paste0("Forest/non-forest (",round(minheight_canopy,1),"m,",round(100*mincover_canopy),"%)\n"))
    plot(patches_analysis, "isforest", legend = "bottomleft", plg = list(bg = "white", bty = "o", box.col = "black"), add = T)
    
    # reset
    if(plot_layer != T){
      dev.off()
    } else {
      par(mfrow = par_mfrow_previous)
    }
  } 
  return(patches_analysis)
}

#%%%%%%%%%%%%%%%%%%%%%#
#### 2. Execution  ####
#%%%%%%%%%%%%%%%%%%%%%#

# output directory
dir_output = "03_classification"
file_chm = "01_chm/chm_lspikefree_multi3.1_slope1.75_offset2.1.tif"

# parameters
minheight_canopy = 5 # in m, height that defines canopy cover (FAO default 2m)
mincover_canopy = 0.05 # fraction, minimum fraction of canopy cover that defines forest (FAO default 0.1)
minarea_forest = 5000 # in m, minimum scale at which something becomes a forest (FAO default 0.5 ha)
buffer_scan = 50 # buffer around scan area (avoid edge effects)
minarea_patch = 250000 # minimum size of analysis patches (25ha)
maxarea_hole = 1000000 # overall maximum area for holes within forest patches (100ha, ~maximum gap size, should also be valid for boreal forests, cf. Goodbody et al. 2020)
ratio_patch_to_hole = 5 # holes within patches are filled as long as they are 5 times smaller than the patch itself
factor_preagg = 5 # aggregation factor for rasters

#%%%%%%%%%%%%%%%%%%%%%%%#
#### 3. Non-adaptive ####
#%%%%%%%%%%%%%%%%%%%%%%%#

# some of the above parameters should be chosen according to study site, alternatively use the adaptive version

# create analysis patches
path_output = file.path(dir_output,"patches_analysis")
if(!dir.exists(path_output)) dir.create(path_output)
path_output_cover = path_output
dir_figures = path_output

factor_refine = 0.0 # if > 0.0, the mincover_canopy parameter is first used to delineate potential forest area; then the median canopy cover of the potential forest area is calculated and multiplied with the factor_adaptive to update the mincover_canopy 

patches_analysis = classify.forest(file_chm = file_chm, path_output = path_output, path_output_cover = path_output_cover, minheight_canopy = minheight_canopy, mincover_canopy = mincover_canopy, minarea_forest = minarea_forest, factor_refine = factor_refine, buffer_scan = buffer_scan, minarea_patch = minarea_patch, maxarea_hole = maxarea_hole, ratio_patch_to_hole = ratio_patch_to_hole, factor_preagg = factor_preagg, plot_layer = file.path(dir_figures,"patches_analysis.pdf"))

#%%%%%%%%%%%%%%%%%%%#
#### 4. Adaptive ####
#%%%%%%%%%%%%%%%%%%%#

path_output = file.path(dir_output,"patches_analysis_adaptive")
if(!dir.exists(path_output)) dir.create(path_output)
path_output_cover = path_output
dir_figures = path_output

factor_refine = 0.5 # if > 0.0, the mincover_canopy parameter is first used to delineate potential forest area; then the median canopy cover of the potential forest area is calculated and multiplied with the factor_adaptive to update the mincover_canopy 

patches_analysis_adaptive = classify.forest(file_chm = file_chm, path_output = path_output, path_output_cover = path_output_cover, minheight_canopy = minheight_canopy, mincover_canopy = mincover_canopy, minarea_forest = minarea_forest, factor_refine = factor_refine, buffer_scan = buffer_scan, minarea_patch = minarea_patch, maxarea_hole = maxarea_hole, ratio_patch_to_hole = ratio_patch_to_hole, factor_preagg = factor_preagg, plot_layer = file.path(dir_figures,"patches_analysis_adaptive.pdf"))

