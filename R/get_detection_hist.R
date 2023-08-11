## this for getting a detection history of individuals
library(sp)
load("./clean_data/Guanecaste/CT_loc.rda")
load("./clean_data/Guanecaste/jaguar.rda")

config <- jsonlite::fromJSON("config.json")


# deal with dates
date_range <- config$date |> as.Date() # setup a date range
occ_days <- config$occasion # number of days in a secondary occasion
M <- config$max_individual # number of augmented individuals


# deal with dates
combined_CT <- lapply(camera_sites, as.data.frame) |>
  Reduce(f = rbind)

combined_CT$Setup_date <- as.Date(combined_CT$Start_Date.Time, 
                                  format = "%d-%m-%y")
combined_CT$Retrieval_date <- as.Date(combined_CT$End_Date.Time, 
                                      format = "%d-%m-%y")
combined_CT <- combined_CT[(combined_CT$Setup_date >= date_range[1])&
                            (combined_CT$Retrieval_date <= date_range[2]),]
jaguar_det <- as.data.frame(jaguar)
jaguar_det$Date.Time <- as.Date(jaguar_det$Date.Time, 
                            format = "%m/%d/%Y")
jaguar_det <- jaguar_det[(jaguar_det$Date.Time <= date_range[2]) & 
                   (jaguar_det$Date.Time >= date_range[1]),]
jaguar_det <- jaguar_det[jaguar_det$Individual.ID!="Undetermined",]

# get dictionaries for numerical ids of stations and individuals
trap_ids <- data.frame(unique(combined_CT[,c("StationName","x","y")]))
trap_ids$id <- 1:nrow(trap_ids)

ind_ids <- data.frame(jaguar = unique(jaguar_det$Individual.ID))
ind_ids$id <- 1:nrow(ind_ids)

n_traps <- nrow(trap_ids)
n_inds <- nrow(ind_ids)

# get occasions
n_occ <- ceiling(as.numeric(max(combined_CT$Retrieval_date)-
                                min(combined_CT$Setup_date))/
                     occ_days) # how many occasions?

occ_end_dates <- min(combined_CT$Setup_date) + (1:n_occ) * occ_days
occ_start_dates <- occ_end_dates - 7

# get trap deployment matrix (row as trap, col as occasion, 1 for active)
trap_mat <- matrix(0, n_traps, n_occ)

## heck, for loops
for(i in 1:nrow(combined_CT)){
  deployment <- combined_CT[i,]
  station <- trap_ids$id[trap_ids$Station==deployment$Station]
  occ_in <- max(which(occ_start_dates <= deployment$Setup_date)) : 
    min(which(occ_end_dates >= deployment$Retrieval_date)) # hack-y way to determine activity in occasions
  trap_mat[station, occ_in] <- 1
}
#image(trap_mat)

# get capture matrix
#cap_mat <- data.frame(ind_id = rep(0,nrow(jaguar_det)), occa = 0, trap = 0)
cap_mat <- array(0,c(M, n_traps, n_occ))
for(i in 1:nrow(jaguar_det)){
  detection <- jaguar_det[i,c("Date.Time", "x","y","Individual.ID","StationName")]
  ind_id <- ind_ids$id[ind_ids$jaguar==detection$Individual.ID]
  occ_in <- which((detection$Date.Time <= occ_end_dates) & 
                          (detection$Date.Time > occ_start_dates))[1]
  the_trap <- trap_ids$id[trap_ids$StationName==detection$StationName]
  if(length(the_trap)>1){
    the_trap <- trap_ids$id[(trap_ids$x==detection$x) & 
                              (trap_ids$StationName==detection$StationName)]
  }
  cap_mat[ind_id, the_trap, occ_in] <- 1
  
}

jaguar_trap_mats <- list(cap_mat = cap_mat, trap_mat = trap_mat,
                    ids = list(trap_ids = trap_ids, ind_ids = ind_ids),
                    occa_design = list(start_date = occ_start_dates, 
                                       occ_days = occ_days, 
                                       range = date_range
                                       )
                    )

save(jaguar_trap_mats, 
  file = "./clean_data/Guanecaste/jaguar_trap_mats_scr.rda")
