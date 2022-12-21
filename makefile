./res/scr_stan_fit.rda: ./clean_data/scr_stan_data.rda
	Rscript ./R/run_SCR.R

./clean_data/scr_stan_data.rda: ./clean_data/jaguar_trap_mats_scr.rda \
																./clean_data/lulc_2017.tif
	Rscript ./R/prepare_scr.R

./clean_data/jaguar_trap_mats_scr.rda: ./clean_data/CT_loc.rda	./clean_data/jaguar.rda
	Rscript ./R/get_detection_hist.R
	
./clean_data/lulc_2017.tif: ./data/Camera_loc/Costa*	\
											./data/detections/JaguarEvents2015-2021_EDIT.csv \
											./data/mask/NASA_OSA/lulc_2017.tif 
	Rscript ./R/pre_processing_loc.R

./clean_data/CT_loc.rda: ./data/Camera_loc/Costa*	\
											./data/detections/JaguarEvents2015-2021_EDIT.csv \
											./data/mask/NASA_OSA/lulc_2017.tif 
	Rscript ./R/pre_processing_loc.R

