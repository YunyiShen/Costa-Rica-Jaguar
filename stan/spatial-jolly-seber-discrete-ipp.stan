// Jolly-Seber model as a multistate hidden Markov model
// adapted from stan-dev/example-models/blob/master/BPA/Ch.10/js_ms.stan
//--------------------------------------
// States:
// 1 not yet entered
// 2 alive
// 3 dead

functions {
  real log_sum_exp_vec(vector x){
    real max_x = max(x);
    return max_x + log(sum(exp(x - max_x)));
  }

  vector rowSum(matrix x){
    int nrow = rows(x);
    int ncol = cols(x);
    vector[norw] res;
    for (i in 1:nrow){
      res[i] = sum(x[i,]);
    }
  }

  int sample_multinomial_from_log_prob(vector x){
    vector tmp = x - log_sum_exp_vec(x);
    return categorical_rng(exp(tmp),1);
  }
}

data {
  int<lower = 1> M;
  int<lower = 1> T;
  int<lower = 1> Kmax;
  int<lower = 1> n_trap;
  //int<lower = 0> K[n_trap,T];
  array int<lower = 0, upper = 1>[M, n_trap,Kmax, T] y; // detection history, it is ok to ignore time order etc. as we assume secondary occasions are exchangable
  array int<lower = 0, upper = 1>[n_trap,Kmax, T] deploy; //deployment history of all camera traps at all occasions )(primary and secondary)
  matrix[n_trap, 2] X; //trap locations

  //delt with environmental covariates, set up a grid for Poisson intensity approximations
  int<lower = 1> n_grid;
  matrix[n_grid, 2] grid_pts;
  //real pixel_area; // sort of absorbed into psi

  // get enviromental variables at grid points
  int<lower = 1> n_env;
  matrix[n_grid, n_env] envX[T]; // environmental variables at each primary occasion, NO NEED OF intercept!
}

transformed data {
  int Tm1 = T - 1;
  int Tp1 = T + 1;

  // deal with environmental covariates
  //real log_pixel_area = log(pixel_area);
  matrix[n_grid, n_trap] sq_dist;
  for (i in 1:n_grid) {
    for (j in 1:n_trap) {
      sq_dist[i, j] = squared_distance(grid_pts[i, ], X[j, ]);
    }
  }
}

parameters {
  // JC part
  vector<lower = 0,upper = 1>[Tm1] gamma;         // recruitment
  real<lower = 0, upper = 1> psi;                 // initial pr. of being alive
  real<lower = 0, upper = 1> phi;                 // survival
  real p0;
  real<lower = 0> alpha1;

  // deal with environmental covariates
  vector[n_env+1] beta_env; // regression coefficients for activity centers
}

transformed parameters {
  real log_p0 = log(p0);
  real<upper = 0> log_recruit[T]; // log transition from unseen to alive
  real<upper = 0> log_survival[T,2]; // log transition from alive to alive
  // dead can only go dead
  vector[M] log_lik;

  // Define probabilities of state S(t+1) given S(t)
  log_recruit[1] = log(psi);
  log_survival[1] = 0;
  
  for (t in 2:T) {
      log_recruit[t] = log(gamma[t-1]);
      log_survival[t] = log(phi);
  }

  { // begin local scope
    real negINF = -1e20; // local negative inf, to work in log scale avoiding underflow, such a effective negInf is needed to avoid autograd error caused by exp(-Inf) 
    vector[n_grid] log_seen_center_grid;// detection probability of a given trap, given an individual is at a given activity center
    real acc1;
    real acc2[2];
    real acc3[2];
    real log_obs;
    real log_obs_nothere; // more a place holder to check if we have seen an individual at a trap when it should NOT be seen
    vector[3] gam[Tp1];
    real<upper = 0> log_po[n_trap, T];// detection probability, given the state is alive, otherwise is -inf (emission in HMM sense), details later

    vector[n_grid] log_mu; // log Poisson intensity, at pixels, up to a constant (we do not need it since we condition on numer of total points)
    vector[n_grid] log_probs;

    // assume individuals do not move in primary occasions, calculate the probability each individual show up in a certain pixel
    for (t in 1:T) {// loop over primary occasions first
      log_mu = envX[t] * beta_env;// intensity at that year
      log_probs = log_mu - log_sum_exp_vec(log_mu);// a discrete distribution over pixels, probability s at each pixel
    

      // below calculates probability of seeing an individual at given trap, in a given primary occasion, condition on state of that individual being alive, otherwise on just cannot see it. We will assume secondary occasions and individuals are exchangable, thus we will not index by individuals and secondary occasions
      for (i in 1:M) {
          for (j in 1:n_trap) {
            for (l in 1:n_grid) { // marginalize over activity centers
                log_seen_center_grid[l] = log_inv_logit( log_p0 - alpha1 * sq_dist[l,j] );// distance sampling 
                log_seen_center_grid[l] += log_probs[l];// add Poisson intensity, essentially weighted sum so that we marginalize over activity centers
            }
            log_po[j,t] = log_sum_exp_vec(log_seen_center_grid);// margianlize over activity centers by logsumexp
          }
      }
    }
    
    // Forward algorithm
    for (i in 1:M) {
      // All individuals are in state 1 (not recruited) at t=0, we work in log scale
      gam[1, 1] = 0;
      gam[1, 2] = negINF; // not sure if this is a good idea but effectively -inf, exp(-inf) can cause problem in Stan due to autograd
      gam[1, 3] = negINF;

      // we iterate to T + 1, because we inserted a dummy period where 
      // every individual is in the "not recruited" state
      for (t in 2:(Tp1)) {
        // likelihood of seeing/not seeing the individual when it is alive
        log_obs = 0;
        log_obs_nothere = 0;
        // loop over secondary occasions and traps
        for (occasion in 1:Kmax) {
          for(j in 1:n_trap){
            log_obs += deploy[j, occasion, t - 1] * (y[i, j, occasion, t - 1] * log_po[j, t-1] + 
            (1-y[i, j, occasion, t - 1]) * log1m_exp(log_po[j, t-1]));
            log_obs_nothere += y[i, j, occasion, t - 1] * negINF;// should not show up when it is not there
          }
        }
        // change to un-entered
          // get from un-entered
        acc1 = gam[t - 1, 1] + log1m_exp(log_recruit[t - 1]);
        acc1 += log_obs_nothere;
        gam[t, 1] = acc1;
        // change to alive
          // come from un-entered
        acc2[1] = gam[t - 1, 1] + log_recruit[t - 1];// recruited into
          // come from alive
        acc2[2] = gam[t - 1, 2] + log_survival[t - 1];// survived
        acc2 += log_obs;
        gam[t, 2] = log_sum_exp_vec(acc2);
        // change to dead
          // come from alive
        acc3[1] = gam[t - 1, 2] + log1m_exp( log_survival[t - 1] );// dead last year
          // come from dead
        acc3[2] = gam[t - 1, 3]; // dead can only be dead
        acc3 += log_obs_nothere; // should not see if dead
        gam[t, 3] = log_sum_exp_vec(acc3);
      }

      log_lik[i] = log_sum_exp_vec(gam[Tp1]);
    }
  }
}

model {
  p0 ~ normal(0, 3);
  alpha1 ~ lognormal(0, 1);
  beta_env ~ normal(0, 10);
  target += sum(log_lik);
}

// note that we could sample z and N in the generated quantities block

generated quantities {
  // TODO: FFBS/Viterbi algorithm to generate z
  array[M, T] int<lower=0, upper = 3> z; // Latent state of all individuals at all primary occasions
  // TODO: generate s from the marginalization
  array[M, T] int<lower = 1, upper = n_grid> s; // Activity center for all individuals at all primary occasions 
  { // begin local stuff

    matrix[n_grid,n_trap] log_lik_center_grid; // we need to save all traps, given we are not marginalize over grids any more 
    vector[n_grid] log_mu; // log Poisson intensity, at pixels, up to a constant (we do not need it since we condition on numer of total points)
    vector[n_grid] log_probs;
    real<upper = 0> log_po[n_trap, T];// detection probability, 


    // local stuff for FFBS/Viterbi
    real negINF = -1e20; // local negative inf, to work in log scale avoiding underflow, such a effective negInf is needed to avoid autograd error caused by exp(-Inf) 
    vector[n_grid] log_seen_center_grid;
    real acc1;
    real acc2[2];
    real acc3[2];
    real log_obs;
    real log_obs_nothere; // more a place holder to check if we have seen an individual at a trap when it should NOT be seen
    vector[3] gam[Tp1];
    vector int<lower = 1, upper = 3>[3] back_ptr[Tp1];

    real logp_state; 



    // construct emission probablity, as well as activity center for each individual 
    for (i in 1:M) {
      for (t in 1:T) {// loop over primary occasions first
        log_mu = envX[t] * beta_env;// intensity at that year
        log_probs = log_mu - log_sum_exp_vec(log_mu);// a discrete distribution over pixels, probability s at each pixel
    

      // below calculates probability of seeing an individual at given trap, in a given primary occasion, condition on state of that individual being alive, otherwise one just cannot see the individual. We will assume secondary occasions and individuals are exchangable, thus we will not index by individuals and secondary occasions
        for (j in 1:n_trap) {
          for (l in 1:n_grid) { // marginalize over activity centers
              log_seen_center_grid[l] = log_inv_logit( log_p0 - alpha1 * sq_dist[l,j] );// distance sampling 
              log_obs = 0;// reuse this, probability of seening the detection history at that trap given the individual is alive in that primary occasion 
              for (occasion in 1:Kmax) { // product over all secondary occasions
                log_obs += deploy[j, occasion, t] * (y[i, j, occasion, t ] * log_seen_center_grid[l] + (1-y[i, j, occasion, t]) * log1m_exp(log_seen_center_grid[l]));
              } // likelihood of the detection history if the individual is alive in that primary occasion, and the activity center is at that pixel
              log_seen_center_grid[l] += log_probs[l];// add Poisson intensity, essentially weighted sum so that we marginalize over activity centers
              log_lik_center_grid[l,j] = log_seen_center_grid[l] + log_obs;// save for later use   
          }
          log_po[j,t] = log_sum_exp_vec(log_seen_center_grid);// margianlize over activity centers by logsumexp
        }

        // now sum over traps using log_lik_center_grid, we have the probability of seeing the detection history at that primary occasion, given the individual is alive in that primary occasion locating at the certain pixel

        s[i,t] = sample_multinomial_from_log_prob(rowSum(log_lik_center_grid)); // sample activity center for each individual at each primary occasion
      }
    } // end dealing with observation and activaty center

    // Viterbi for latent state
    for (i in 1:M) {
      // All individuals are in state 1 (not recruited) at t=0, we work in log scale
      gam[1, 1] = 0;
      gam[1, 2] = negINF; // not sure if this is a good idea but effectively -inf, exp(-inf) can cause problem in Stan due to autograd
      gam[1, 3] = negINF;
      back_ptr[1,] = negINF;

      // we iterate to T + 1, because we inserted a dummy period where 
      // every individual is in the "not recruited" state
      // forward algorithm, similar to calculating the likelihood
      for (t in 2:(Tp1)) {
        // likelihood of seeing/not seeing the individual when it is alive
        log_obs = 0;
        log_obs_nothere = 0;
        // loop over secondary occasions and traps
        for (occasion in 1:Kmax) {
          for(j in 1:n_trap){
            log_obs += deploy[j, occasion, t - 1] * (y[i, j, occasion, t - 1] * log_po[j, t-1] + 
            (1-y[i, j, occasion, t - 1]) * log1m_exp(log_po[j, t-1]));
            log_obs_nothere += y[i, j, occasion, t - 1] * negINF;// should not show up when it is not there
          }
        }
        // change to un-entered
          // get from un-entered
        acc1 = gam[t - 1, 1] + log1m_exp(log_recruit[t - 1]);
        acc1 += log_obs_nothere;
        gam[t, 1] = acc1;

        back_ptr[t, 1] = 1;// state decoding

        // change to alive
          // come from un-entered
        acc2[1] = gam[t - 1, 1] + log_recruit[t - 1];// recruited into
          // come from alive
        acc2[2] = gam[t - 1, 2] + log_survival[t - 1];// survived
        acc2 += log_obs;

        if(acc2[2]>acc2[1]){
          back_ptr[t, 2] = 2;
          gam[t, 2] = acc2[2];

        } else {
          back_ptr[t, 2] = 1;
          gam[t, 2] = acc2[1];
        }

        // change to dead
          // come from alive
        acc3[1] = gam[t - 1, 2] + log1m_exp( log_survival[t - 1] );// dead last year
          // come from dead
        acc3[2] = gam[t - 1, 3]; // dead can only be dead
        acc3 += log_obs_nothere; // should not see if dead
        if(acc3[2]>acc3[1]){
          back_ptr[t, 3] = 3;
          gam[t, 3] = acc3[2];
        } else {
          back_ptr[t, 3] = 2;
          gam[t, 3] = acc3[1];
        }

      }// end forward algorithm

      // backward pass
      // we start at the last period, and work our way back
      logp_state = max(gamma[Tp1,]);
      for(ss in 1:3){
        if(gam[Tp1,ss] == logp_state){
          state[i,T] = ss;
        }
      }


      for (tt in (T - 1):1) {
        state[i, tt] = back_ptr[tt, state[i, tt+ 1]]; // back_ptr is already indexed +1
      }
    } // end Viterbi
  }// end local stuff
}