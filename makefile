### Guanecaste SCR

./res/Guanecaste_scr_stan_fit.rda: ./clean_data/Guanecaste/scr_stan_data.rda
	Rscript ./R/run_SCR.R
	rm clean_data/Guanecaste/*.rda

./clean_data/Guanecaste/scr_stan_data.rda: ./clean_data/Guanecaste/jaguar_trap_mats_scr.rda \
								./clean_data/Guanecaste/grid_objs_data.rda
	Rscript ./R/prepare_scr.R

./clean_data/Guanecaste/grid_objs_data.rda: ./clean_data/Guanecaste/jaguar_trap_mats_scr.rda \
																 ./clean_data/shoreline.tif \
																 ./config.json
	Rscript ./R/form_grids.R

./clean_data/Guanecaste/jaguar_trap_mats_scr.rda: ./clean_data/Guanecaste/CT_loc.rda	./clean_data/Guanecaste/jaguar.rda
	Rscript ./R/get_detection_hist.R

./clean_data/Guanecaste/jaguar_trap_mats_js.rda: ./clean_data/Guanecaste/CT_loc.rda	./clean_data/Guanecaste/jaguar.rda
	Rscript ./R/get_detection_hist_js.R
	
./clean_data/shoreline.tif: ./data/Guanecaste/shore_line/CRI_adm0*	\
							./data/Guanecaste/jaguar_capture.csv \
							./data/Guanecaste/Guanecaste_trap_loc.csv 
	Rscript ./R/pre_processing_loc.R

./clean_data/Guanecaste/CT_loc.rda: ./data/Guanecaste/shore_line/CRI_adm0*	\
							./data/Guanecaste/jaguar_capture.csv \
							./data/Guanecaste/Guanecaste_trap_loc.csv 
	Rscript ./R/pre_processing_loc.R

./res/two_beach_model.rda: ./clean_data/Guanecaste/jaguar_trap_mats_scr.rda \
							./clean_data/Guanecaste/jaguar.rda\
							./clean_data/Guanecaste/jaguar_trap_mats_scr.rda
	Rscript ./R/two_beach_model.R
	rm clean_data/Guanecaste/*.rda