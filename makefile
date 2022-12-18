./clean_data/scr_stan_data.rda: ./clean_data/jaguar_trap_mats_scr.rda
	Rscript ./R/prepare_scr.R

./clean_data/jaguar_trap_mats_scr.rda: ./clean_data/CT_loc.rda	./clean_data/jaguar.rda
	Rscript ./R/get_detection_hist.R
	
./clean_data/CT_loc.rda: ./data/Costa*	./data/JaguarEvents2015-2021_EDIT.csv
	Rscript ./R/pre_processing_loc.R

