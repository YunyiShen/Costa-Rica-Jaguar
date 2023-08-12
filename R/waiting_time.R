# compare two waiting times: 
#   waiting time towards next visit by same individual
#   waiting time towards change of ownership, 
#     e.g. the detection history being 3,1,1,1,2, 
#      we calculate the time difference btween the second and the 5th detection
# TODO: think about how to deal with censoring on the right, 
#       some sites might only 
#       be visited by an individual once 
#       and the waiting time are actually censored but not 0


get_waiting_time_same_inid <- function(dets_in, stations_in){
    waiting_time_next_visit <- c()
    time_visit_by_other_inid <- c()
    for(station in stations_in$StationName){
        dets_in_station <- dets_in[dets_in$StationName==station,]
        dets_in_station <- dets_in_station[order(dets_in_station$Date.Time), ]
        inds_in <- unique(dets_in_station$Individual.ID)
        for(ind in inds_in){
            dets_this_ind <- dets_in_station[dets_in_station$Individual.ID==ind,]
            N_record <- nrow(dets_this_ind)
            if(N_record>1){
                this_waiting_starting_time <- dets_this_ind$Date.Time[1]
                for(i in 2:N_record){
                    waiting_time_next_visit <- c(waiting_time_next_visit, 
                                                 dets_this_ind$Date.Time[i] - this_waiting_starting_time)
                    this_waiting_starting_time <- dets_this_ind$Date.Time[i]

                    # other records within the time period
                    other_records <- dets_in_station$Individual.ID[dets_in_station$Date.Time<=dets_this_ind$Date.Time[i] & dets_in_station$Date.Time>=this_waiting_starting_time]
                    time_visit_by_other_inid <- c(time_visit_by_other_inid, sum(other_records!=ind))
                    #time_visit_by_other_inid <- c(time_visit_by_other_inid, length(unique(other_records))-1)
                }

            }
        }

    }
    return(list(waiting = waiting_time_next_visit, 
                counts = time_visit_by_other_inid))
}

get_waiting_time_change_ownership <- function(dets_in, stations_in){
    waiting_time_next_change_owner <- c()
    time_visit_by_current_owner <- c()
    for(station in stations_in$StationName){
        dets_in_station <- dets_in[dets_in$StationName==station,]
        # to order the records
        dets_in_station <- dets_in_station[order(dets_in_station$Date.Time), ]
        N_record <- nrow(dets_in_station)
        n_uniq_inid <- length(unique(dets_in_station$Individual.ID))
        if(n_uniq_inid>1){
            current_owner <- dets_in_station$Individual.ID[1]
            current_ownership_starts <- dets_in_station$Date.Time[1]
            visit_counts <- 0
            for(i in 2:N_record){
                current_visitor <- dets_in_station$Individual.ID[i]
                if(current_visitor!=current_owner){
                    waiting_time_next_change_owner <- c(waiting_time_next_change_owner, 
                                                        dets_in_station$Date.Time[i] - current_ownership_starts)
                    current_owner <- current_visitor
                    current_ownership_starts <- dets_in_station$Date.Time[i]
                    time_visit_by_current_owner <- c(time_visit_by_current_owner, visit_counts)
                    visit_counts <- 0
                } else {
                    visit_counts <- visit_counts + 1
                }
            }
        }
        return(list(waiting = waiting_time_next_change_owner,
                counts = time_visit_by_current_owner))
    }
}


get_waiting_time_same_sex <- function(dets_in, stations_in, sex = "Female"){
    waiting_time_next_visit <- c()
    time_visit_by_other_sex <- c()
    for(station in stations_in$StationName){
        dets_in_station <- dets_in[dets_in$StationName==station,]
        dets_in_station <- dets_in_station[order(dets_in_station$Date.Time), ]
        
        dets_this_sex <- dets_in_station[dets_in_station$Sex==sex,]
        N_record <- nrow(dets_this_sex)
        if(N_record>1){
            
            this_waiting_starting_time <- dets_this_sex$Date.Time[1]
            for(i in 2:N_record){
                waiting_time_next_visit <- c(waiting_time_next_visit, 
                                             dets_this_sex$Date.Time[i] - this_waiting_starting_time)
                this_waiting_starting_time <- dets_this_sex$Date.Time[i]

                # other records within the time period
                other_records <- dets_in_station$Sex[dets_in_station$Date.Time<=dets_this_sex$Date.Time[i] & dets_in_station$Date.Time>=this_waiting_starting_time]
                time_visit_by_other_sex <- c(time_visit_by_other_sex, sum(other_records!=sex))
                #time_visit_by_other_inid <- c(time_visit_by_other_inid, length(unique(other_records))-1)
            }

        }

    }
    return(list(waiting = waiting_time_next_visit, 
                counts = time_visit_by_other_sex))
}

get_waiting_time_change_sex <- function(dets_in, stations_in, sex = "Female"){
    waiting_time_next_change_owner <- c()
    time_visit_by_current_owner <- c()
    for(station in stations_in$StationName){
        dets_in_station <- dets_in[dets_in$StationName==station,]
        # to order the records
        dets_in_station <- dets_in_station[order(dets_in_station$Date.Time), ]
        
        transition_events <- which(dets_in_station$Sex[2:nrow(dets_in_station)-1]==sex & 
              dets_in_station$Sex[2:nrow(dets_in_station)-1]!= dets_in_station$Sex[2:nrow(dets_in_station)])
        if(length(transition_events)>1){
            waiting_time <- diff(dets_in_station$Date.Time[transition_events])
            num_visit <- diff(transition_events)
            waiting_time_next_change_owner <- c(waiting_time_next_change_owner, waiting_time)
            time_visit_by_current_owner <- c(time_visit_by_current_owner, num_visit)
        }
        

        return(list(waiting = waiting_time_next_change_owner,
                counts = time_visit_by_current_owner))
    }
}







load("./clean_data/Guanecaste/CT_loc.rda")
load("./clean_data/Guanecaste/jaguar.rda")
## some setting ups
beach <- c("Naranjo", "Nancite")
stations_in <- as.data.frame(jaguar)[,c("StationName", "Beach")][jaguar$Beach %in% beach,]
stations_in <- stations_in[!duplicated(stations_in$StationName),]
dets_in <- as.data.frame(jaguar)[jaguar$Beach %in% beach,]
dets_in$Date.Time <- as.Date(dets_in$Date.Time, format = "%m/%d/%Y")



#### waiting time towards next visit ####
waiting_time_next_visit_female <- 
    get_waiting_time_same_inid(dets_in[dets_in$Sex=="Female",], 
                                stations_in)

waiting_time_next_visit_male <- 
    get_waiting_time_same_inid(dets_in[dets_in$Sex=="Male",], 
                                stations_in)


#### waiting time towards change of ownership ####
waiting_time_change_owner_female <- 
    get_waiting_time_change_ownership(dets_in[dets_in$Sex=="Female",], 
                                stations_in)

waiting_time_change_owner_male <- 
    get_waiting_time_change_ownership(dets_in[dets_in$Sex=="Male",], 
                                stations_in)

### make some plots ###
png("./res/Figs_Guanecaste/waiting_time.png", width = 4, height = 3, units = "in", res = 300)
par(mar = c(2.5,2.5,1,.5), mgp = c(1.5, 0.5, 0))
plot(ecdf(waiting_time_next_visit_female$waiting), main = "", 
     xlab = "Waiting time (days)", ylab = "Empirical CDF", xlim = c(0,150))
plot(ecdf(waiting_time_next_visit_male$waiting), add = T, col = "red", lty = 2)
plot(ecdf(waiting_time_change_owner_female$waiting), add = T, col = "blue", lty = 3)
plot(ecdf(waiting_time_change_owner_male$waiting), add = T, col = "#035c03", lty = 4)
legend(x="bottomright", legend = c("female, revisit", 
                "male, revisit", 
                "female, switch individual", 
                "male, switch individual"), 
                col = c("black", "red", "blue", "#035c03"),
                lty = c(1,2,3,4), pch = rep(16,4))
dev.off()

png("./res/Figs_Guanecaste/counts.png", width = 6, height = 4.5, units = "in", res = 300)
par(mar = c(2.5,2.5,1,.5), mgp = c(1.5, 0.5, 0), mfrow = c(2,2))
hist(waiting_time_next_visit_female$counts, freq = F, 
        breaks = 0:14-.5, main = "", xlab = "# visits by others before revisit, female", 
        ylim = c(0,1), ylab = "Proportion")
hist(waiting_time_next_visit_male$counts, breaks = 0:14-.5, freq = F, 
      main = "",  xlab = "# visits by others before revisit, male", 
        ylim = c(0,1), ylab = "Proportion")
hist(waiting_time_change_owner_female$counts, breaks = 0:14-.5, freq = F, 
      main = "",  xlab = "# revisits before switching, female", 
        ylim = c(0,1), ylab = "Proportion")
hist(waiting_time_change_owner_male$counts, breaks = 0:14-.5, freq = F,
      main = "",  xlab = "# revisits before switching, male", 
        ylim = c(0,1), ylab = "Proportion")

dev.off()


### hacky way to do sex encounter
dets_in_hacky <- dets_in
dets_in_hacky$Individual.ID <- dets_in_hacky$Sex

waiting_same_sex <- 
    get_waiting_time_same_inid(dets_in_hacky, 
                                stations_in)

waiting_diff_sex <- 
    get_waiting_time_change_ownership(dets_in_hacky, 
                                stations_in)
png("./res/Figs_Guanecaste/waiting_time_between_sex.png", width = 4, height = 3, units = "in", res = 300)
par(mar = c(2.5,2.5,1,.5), mgp = c(1.5, 0.5, 0))
plot(ecdf(waiting_same_sex$waiting), main = "", 
     xlab = "Waiting time (days)", ylab = "Empirical CDF", xlim = c(0,100))
plot(ecdf(waiting_diff_sex$waiting), add = T, col = "red", lty = 2)
legend(x="bottomright", legend = c("visit by same sex", 
                "visit by different sex"
                ), 
                col = c("black", "red"),
                lty = c(1,2), pch = rep(16,2))
dev.off()

png("./res/Figs_Guanecaste/counts_by_sex.png", width = 6, height = 2.25, units = "in", res = 300)

par(mar = c(2.5,2.5,1,.5), mgp = c(1.5, 0.5, 0), mfrow = c(1,2))
hist(waiting_same_sex$counts, freq = F, 
        breaks = 0:10-.5, main = "", xlab = "# visits by others sex before revisit", 
        ylim = c(0,0.8), ylab = "Proportion")

hist(waiting_diff_sex$counts, breaks = 0:10-.5, freq = F, 
      main = "",  xlab = "# revisits before switching sex", 
        ylim = c(0,0.8), ylab = "Proportion")


dev.off()

### actual by sex comparison
waiting_female_to_female <- 
    get_waiting_time_same_sex(dets_in, 
                                stations_in, 
                                sex = "Female")

waiting_male_to_male <- 
    get_waiting_time_same_sex(dets_in, 
                                stations_in, 
                                sex = "Male")

waiting_female_to_male <- 
    get_waiting_time_change_sex(dets_in, 
                                stations_in, 
                                sex = "Female")

waiting_male_to_female <- 
    get_waiting_time_change_sex(dets_in, 
                                stations_in, 
                                sex = "Male")

png("./res/Figs_Guanecaste/waiting_time_sexes.png", width = 4, height = 3, units = "in", res = 300)
par(mar = c(2.5,2.5,1,.5), mgp = c(1.5, 0.5, 0))
plot(ecdf(waiting_female_to_female$waiting), main = "", 
     xlab = "Waiting time (days)", ylab = "Empirical CDF", xlim = c(0,150))
plot(ecdf(waiting_male_to_male$waiting), add = T, col = "red", lty = 2)
plot(ecdf(waiting_female_to_male$waiting), add = T, col = "blue", lty = 3)
plot(ecdf(waiting_male_to_female$waiting), add = T, col = "#035c03", lty = 4)
legend(x="bottomright", legend = c("female to female", 
                "male to male", 
                "female to male", 
                "male to female"), 
                col = c("black", "red", "blue", "#035c03"),
                lty = c(1,2,3,4), pch = rep(16,4))
dev.off()

png("./res/Figs_Guanecaste/counts_sexes.png", width = 6, height = 4.5, units = "in", res = 300)
par(mar = c(2.5,2.5,1,.5), mgp = c(1.5, 0.5, 0), mfrow = c(2,2))
hist(waiting_female_to_female$counts, freq = F, 
        breaks = 0:14-.5, main = "", xlab = "# of male before female revisit", 
        ylim = c(0,1), ylab = "Proportion")
hist(waiting_male_to_male$counts, breaks = 0:14-.5, freq = F, 
      main = "",  xlab = "# of female before male revisit", 
        ylim = c(0,1), ylab = "Proportion")
hist(waiting_female_to_male$counts, breaks = 0:14-.5, freq = F, 
      main = "",  xlab = "# male revisits before female visit", 
        ylim = c(0,1), ylab = "Proportion")
hist(waiting_male_to_female$counts, breaks = 0:14-.5, freq = F,
      main = "",  xlab = "# female revisits before male visit", 
        ylim = c(0,1), ylab = "Proportion")

dev.off()
