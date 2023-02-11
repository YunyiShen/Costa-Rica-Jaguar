./res/js_lock_stan_fit.rda: ./clean_data/js_stan_data.rda
	Rscript ./R/run_JS_lock.R

./res/js_stan_fit.rda: ./clean_data/js_stan_data.rda
	Rscript ./R/run_JS.R

./clean_data/js_stan_data.rda: ./clean_data/jaguar_trap_mats_js.rda \
								./clean_data/grid_objs_data.rda \
								./clean_data/rasters/2018_wlprai.tif \
								./clean_data/rasters/dist_prkbndry.tif
	Rscript ./R/prepare_js.R

./clean_data/rasters/2018_wlprai.tif: ./data/SpatialCovariates/WLP_RAI_Rasters/2018_wlprai/hdr.adf
	Rscript ./R/pre_processing_rasters.R

./clean_data/rasters/dist_prkbndry.tif: ./data/SpatialCovariates/dist_prkbndry/hdr.adf
	Rscript ./R/pre_processing_rasters.R


./res/scr_stan_fit.rda: ./clean_data/scr_stan_data.rda
	Rscript ./R/run_SCR.R

./clean_data/scr_stan_data.rda: ./clean_data/jaguar_trap_mats_scr.rda \
								./clean_data/grid_objs_data.rda
	Rscript ./R/prepare_scr.R

./clean_data/grid_objs_data.rda: ./clean_data/jaguar_trap_mats_scr.rda ./clean_data/lulc_2017.tif
	Rscript ./R/form_grids.R

./clean_data/jaguar_trap_mats_scr.rda: ./clean_data/CT_loc.rda	./clean_data/jaguar.rda
	Rscript ./R/get_detection_hist.R

./clean_data/jaguar_trap_mats_js.rda: ./clean_data/CT_loc.rda	./clean_data/jaguar.rda
	Rscript ./R/get_detection_hist_js.R
	
./clean_data/lulc_2017.tif: ./data/Camera_loc/Costa*	\
							./data/detections/JaguarEvents2015-2021_EDIT.csv \
							./data/mask/NASA_OSA/lulc_2017.tif 
	Rscript ./R/pre_processing_loc.R

./clean_data/CT_loc.rda: ./data/Camera_loc/Costa*	\
						 ./data/detections/JaguarEvents2015-2021_EDIT.csv \
						 ./data/mask/NASA_OSA/lulc_2017.tif 
	Rscript ./R/pre_processing_loc.R

