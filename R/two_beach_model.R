library(rstan)
options(mc.cores = parallel::detectCores()/2)

load("./clean_data/Guanecaste/CT_loc.rda")
load("./clean_data/Guanecaste/jaguar.rda")
load("./clean_data/Guanecaste/jaguar_trap_mats_scr.rda")
config <- jsonlite::fromJSON("config.json")

get_detection_sum <- function(full_hist, ind_id, trap_id, 
                            max_ind = config$max_individual){
    apply(full_hist[ind_id,trap_id,],1,sum) |>
    c(rep(0, config$max_ind-length(ind_id)))
}

## some setting ups
beach <- c("Naranjo", "Nancite")
stations_in <- as.data.frame(jaguar)[,c("StationName", "Beach")][jaguar$Beach %in% beach,]
stations_in <- stations_in[!duplicated(stations_in$StationName),]

### which station in beach 1
stations_id_beach1 <- jaguar_trap_mats$ids$trap_ids$id[
  jaguar_trap_mats$ids$trap_ids$StationName %in% 
    stations_in$StationName[stations_in$Beach==beach[1]]]

### which station in beach 2
stations_id_beach2 <- jaguar_trap_mats$ids$trap_ids$id[
    jaguar_trap_mats$ids$trap_ids$StationName %in% 
        stations_in$StationName[stations_in$Beach==beach[2]]]

### get_detections
female_id <- jaguar_trap_mats$ids$ind_ids$id[jaguar_trap_mats$ids$ind_ids$sex=="Female"]
male_id <- jaguar_trap_mats$ids$ind_ids$id[jaguar_trap_mats$ids$ind_ids$sex=="Male"]

jaguar_det_male <- jaguar_trap_mats$cap_mat[male_id,,]
jaguar_det_female <- jaguar_trap_mats$cap_mat[female_id,,]

### total number of detections
female_det_sum_beach1 <- get_detection_sum(jaguar_trap_mats$cap_mat, 
                                    female_id, stations_id_beach1)
female_det_sum_beach2 <- get_detection_sum(jaguar_trap_mats$cap_mat, 
                                    female_id, stations_id_beach2)
male_det_sum_beach1 <- get_detection_sum(jaguar_trap_mats$cap_mat, 
                                    male_id, stations_id_beach1)
male_det_sum_beach2 <- get_detection_sum(jaguar_trap_mats$cap_mat, 
                                    male_id, stations_id_beach2)

### total number of trap effort

beach1_deploy_sum <- sum(jaguar_trap_mats$trap_mat[stations_id_beach1,])
beach2_deploy_sum <- sum(jaguar_trap_mats$trap_mat[stations_id_beach2,])

deploy_sum <- c(beach1_deploy_sum, beach2_deploy_sum)

#### run stan model ####
m_init <- stan_model("./stan/two_beach_model.stan")
set.seed(1234)
female_stan_data <- list(deploy_sum = deploy_sum, y1 = female_det_sum_beach1,
                        y2 = female_det_sum_beach2, M = config$max_individual)
female_model <- sampling(m_init,  data = female_stan_data, 
                        iter = 15 * config$stan_iters)

male_stan_data <- list(deploy_sum = deploy_sum, y1 = male_det_sum_beach1,
                        y2 = male_det_sum_beach2, M = config$max_individual)
male_model <- sampling(m_init,  data = male_stan_data, 
                        iter = 15 * config$stan_iters)


save(female_stan_data, female_model,
     male_stan_data, male_model, 
     file = "./res/two_beach_model.rda")


## psi
#stan_dens(male_model, pars = c("psi_reduced"))
#stan_dens(female_model, pars = c("psi_reduced"))

## decoded N
z_male <- extract(male_model, pars = c("state"))$state
NN_male_1 <- rowSums(z_male == 2)
NN_male_2 <- rowSums(z_male == 3)
NN_male_both <- rowSums(z_male == 4)

z_female <- extract(female_model, pars = c("state"))$state
NN_female_1 <- rowSums(z_female == 2)
NN_female_2 <- rowSums(z_female == 3)
NN_female_both <- rowSums(z_female == 4)


png("./res/Figs_Guanecaste/two_beach_model.png", width = 8, height = 3, units = "in", res = 300)
par(mar = c(2.5,2.5,1,.5), mgp = c(1.5, 0.5, 0))
par(mfrow = c(2,4))

hist(NN_female_1, breaks = 0:27-.5, ylim = c(0,0.8),
    freq = FALSE, main = paste("female,", beach[1], "exclusive"), 
    xlab = "# individuals", ylab = "Posterior")
hist(NN_female_2, breaks = 0:27-.5, ylim = c(0,0.8),
    freq = FALSE, main = paste("female,", beach[2], "exclusive"), 
    xlab = "# individuals", ylab = "Posterior")

hist(NN_female_both, breaks = 0:29-.5, ylim = c(0,0.8), 
    freq = FALSE, main = paste("female, both beaches"), 
    xlab = "# individuals", ylab = "Posterior")

hist(NN_female_both + NN_female_1 + NN_female_2, breaks = 5:35-.5, ylim = c(0,0.8), 
    freq = FALSE, main = paste("female, total"), 
    xlab = "# individuals", ylab = "Posterior")


hist(NN_male_1, breaks = 0:27-.5, ylim = c(0,0.8),
    freq = FALSE, main = paste("male,", beach[1], "exclusive"), 
    xlab = "# individuals", ylab = "Posterior")
hist(NN_male_2, breaks = 0:27-.5, ylim = c(0,0.8),
    freq = FALSE, main = paste("male,", beach[2], "exclusive"), 
    xlab = "# individuals", ylab = "Posterior")

hist(NN_male_both, breaks = 0:29-0.5, ylim = c(0,0.8),
    freq = FALSE, main = paste("male, both beaches"), 
    xlab = "# individuals", ylab = "Posterior")

hist(NN_male_both + NN_male_1 + NN_male_2, 
    breaks = 5:35-.5, ylim = c(0,0.8), 
    freq = FALSE, main = paste("male, total"), 
    xlab = "# individuals", ylab = "Posterior")

dev.off()
