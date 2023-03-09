## this for getting a detection history of individuals
library(sp)
library(jsonlite)
load("./clean_data/CT_loc.rda")
load("./clean_data/jaguar.rda")
config <- jsonlite::fromJSON("config.json")


# deal with dates
date_range <- config$date |> as.Date() # setup a date range
occ_days <- config$occasion # number of days in a secondary occasion
M <- config$max_individual # number of augmented individuals


year_range <- format(date_range,"%Y")
primary_sessions <- 1 + as.numeric(year_range[2]) - as.numeric(year_range[1]) # number of primary sessions
years <- year_range[1]:year_range[2] # years in the range


combined_CT <- lapply(camera_sites, as.data.frame) |>
  Reduce(f = rbind)

combined_CT$Setup_date <- as.Date(combined_CT$Setup_date, 
                                  format = "%m/%d/%Y")
combined_CT$Retrieval_date <- as.Date(combined_CT$Retrieval_date, 
                                      format = "%m/%d/%Y")
combined_CT <- combined_CT[(combined_CT$Setup_date >= date_range[1])&
                            (combined_CT$Retrieval_date <= date_range[2]),]
combined_CT$year <- (format(combined_CT$Setup_date,"%Y"))

jaguar_det <- as.data.frame(jaguar)
jaguar_det$Date.Time <- as.Date(jaguar_det$Date.Time, 
                            format = "%m/%d/%Y %H:%M")
jaguar_det <- jaguar_det[(jaguar_det$Date.Time <= date_range[2]) & 
                   (jaguar_det$Date.Time >= date_range[1]),]
jaguar_det$year <- (format(jaguar_det$Date.Time,"%Y"))
jaguar_det <- jaguar_det[jaguar_det$Individual.ID!="Undetermined",]



# get dictionaries for numerical ids of stations and individuals
trap_ids <- data.frame(unique(combined_CT[,c("Station","x","y")]))
trap_ids$id <- 1:nrow(trap_ids)

ind_ids <- data.frame(jaguar = unique(jaguar_det$Individual.ID))
ind_ids$id <- 1:nrow(ind_ids)

n_traps <- nrow(trap_ids)
n_inds <- nrow(ind_ids)

# get # occasions each year
n_occ <- 0 * years
occ_start_dates <- occ_end_dates <- list()

for(i in 1:length(years)){
  year <- years[i]
  tmp <- combined_CT[combined_CT$year==year,]
  n_occ[i] <- ceiling(as.numeric(max(tmp$Retrieval_date)-
                                min(tmp$Setup_date))/
                     occ_days) # how many occasions?
  occ_end_dates[[i]] <- min(tmp$Setup_date) + (1:n_occ[i]) * occ_days
  occ_start_dates[[i]] <- occ_end_dates[[i]] - 7
}

Kmax <- max(n_occ) # max number of secondary occasions

# get trap deployment matrix (row as trap, col as occasion, 1 for active)
trap_mat <- array(0, c(n_traps, Kmax, length(years)))

## heck, for loops
for(i in 1:nrow(combined_CT)){
  deployment <- combined_CT[i,]
  station <- trap_ids$id[trap_ids$Station==deployment$Station]
  year_id <- which(years==deployment$year)
  occ_in <- max(which(occ_start_dates[[year_id]] <= deployment$Setup_date)) : 
    min(which(occ_end_dates[[year_id]] > deployment$Retrieval_date)) # hack-y way to determine activity in occasions
  trap_mat[station, occ_in, year_id] <- 1
}
#image(trap_mat)

# get capture matrix
cap_mat <- array(0,c(M, n_traps, Kmax, length(years)))
for(i in 1:nrow(jaguar_det)){
  detection <- jaguar_det[i,c("Date.Time", "year","x","y","Individual.ID","Station")]
  inid_id <- ind_ids$id[ind_ids$jaguar==detection$Individual.ID]
  year_id <- which(years==detection$year)
  occ_in <- which((detection$Date.Time <= occ_end_dates[[year_id]]) & 
                          (detection$Date.Time > occ_start_dates[[year_id]]))[1]
  trap_id <- trap_ids$id[trap_ids$Station==detection$Station]
  if(length(trap_id)>1){
    trap_id <- trap_ids$id[(trap_ids$x==detection$x) & 
                              (trap_ids$Station==detection$Station)]
  }
  cap_mat[inid_id, trap_id, occ_in, year_id] <- 1
}

jaguar_trap_mats <- list(cap_mat = cap_mat, trap_mat = trap_mat,
                    ids = list(trap_ids = trap_ids, ind_ids = ind_ids),
                    occa_design = list(start_date = occ_start_dates, 
                                       occ_days = occ_days, 
                                       range = date_range
                                       )
                    )

save(jaguar_trap_mats, file = "./clean_data/jaguar_trap_mats_js.rda")

