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
    vector[n_grid] log_seen_center_grid;
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
    

      // below calculates probability of seeing an individual at given trap, in a given primary occasion, condition on state of that individual. We will assume secondary occasions and individuals are exchangable, thus we will not index by individuals and secondary occasions
      for (i in 1:M) {
          for (j in 1:n_trap) {
            for (l in 1:n_grid) { // marginalize over activity centers
                log_seen_center_grid[l] = log_inv_logit( log_p0 - alpha1 * sq_dist[l,j] );// distance sampling 
                log_seen_center_grid += log_probs[l];// add Poisson intensity, essentially weighted sum so that we marginalize over activity centers
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
        acc2[2] = gam[t - 1, 2] + log_survival[t - 1];// survived
        acc2 += log_obs;
        gam[t, 2] = log_sum_exp_vec(acc2);
        // change to dead
          // come from alive
        acc3[1] = gam[t - 1, 2] + log1m_exp( log_survival[t - 1] );// dead last year
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
  //array[M, T] int<lower=1, upper=3> z; // Latent state
  // TODO: FFBS algorithm to generate z
  // TODO: generate s from the marginalization
}