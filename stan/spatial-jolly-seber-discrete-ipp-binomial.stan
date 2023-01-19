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
    vector[nrow] res;
    for (i in 1:nrow){
      res[i] = sum(x[i,]);
    }
    return res;
  }

}

data {
  int<lower = 1> M;
  int<lower = 1> n_trap;
  int<lower = 1> Kmax;
  int<lower = 1> T;
  
  //int<lower = 0> K[n_trap,T];
  int<lower = 0, upper = 1> y[M, n_trap,Kmax, T]; // detection history, it is ok to ignore time order etc. as we assume secondary occasions are exchangable
  int<lower = 0, upper = 1> deploy[n_trap,Kmax, T]; //deployment history of all camera traps at all occasions )(primary and secondary)
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
  //vector<lower = 0,upper = 1>[Tm1] gamma;         // recruitment
  real<lower = 0, upper = 1> gamma;         // recruitment, not time dependent
  real<lower = 0, upper = 1> psi;                 // initial pr. of being alive
  real<lower = 0, upper = 1> phi;                 // survival
  real<lower = 0, upper = 1> p0;
  real<lower = 0> alpha1;

  // deal with environmental covariates
  vector[n_env] beta_env; // regression coefficients for activity centers
}

transformed parameters {
  real log_p0 = log(p0);
  real<upper = 0> log_init_recruit;
  real<upper = 0> log_recruit; // log transition from unseen to alive
  real<upper = 0> log_survival; // log transition from alive to alive
  // dead can only go dead
  vector[M] log_lik;

  // Define probabilities of state S(t+1) given S(t)
  log_init_recruit = log(psi);
  log_survival = log(phi);
  
  log_recruit = log(gamma);
      

  { // begin local scope
    real negINF = -1e20; // local negative inf, to work in log scale avoiding underflow, such a effective negInf is needed to avoid autograd error caused by exp(-Inf) 
    //real log_po[M, T]; // detection probability, given the state is alive, otherwise is -inf (emission in HMM sense), details later
    real acc1;
    vector[2] acc2;
    vector[2] acc3;
    real log_obs;
    real log_obs_nothere; // more a place holder to check if we have seen an individual at a trap when it should NOT be seen
    vector[3] gam[Tp1];
    int flag = 0;

    
    vector[n_grid] log_seen_center_grid;// detection probability of a given trap, given an individual is at a given activity center
    matrix[n_grid, n_trap] dist_tmp;
    matrix[n_grid, T] log_probs; // log Poisson intensity, at pixels, up to a constant (we do not need it since we condition on numer of total points)
    
    dist_tmp = log_p0 - alpha1 * sq_dist;
    for(t in 1:T){
      log_probs[:, t] = log_softmax( envX[t] * beta_env);// intensity at that year
    }

    

    
    // Forward algorithm
    for (i in 1:M) {
      // All individuals are in state 1 (not recruited) at t=0, we work in log scale
      gam[1, 1] = 0;
      gam[1, 2] = negINF; // not sure if this is a good idea but effectively -inf, exp(-inf) can cause problem in Stan due to autograd
      gam[1, 3] = negINF;
      // manually set the first period's recruitment
      
      log_seen_center_grid = rep_vector(0, n_grid);
      for (j in 1:n_trap) {
        
        log_seen_center_grid += sum(y[i, j, :, 1]) * col( dist_tmp, j )+ 
              (sum(deploy[j, :, 1 ])-sum(y[i, j, :, 1])) * 
              log1m_exp(col( dist_tmp, j )); // prob seeing the detection history at a given activaty center
            
      }
      log_obs = log_sum_exp_vec(log_seen_center_grid + log_probs[:,1]);// margianlize over activity centers by logsumexp
      log_obs_nothere = sum( to_matrix(y[i, :, :, 1])) * negINF;
      // change to un-entered
        // get from un-entered
      acc1 = gam[1, 1] + log1m_exp(log_init_recruit);
      acc1 += log_obs_nothere;
      gam[2, 1] = acc1;
        
      // change to alive
        // come from un-entered
      acc2[1] = gam[1, 1] + log_init_recruit;// recruited into
        // come from alive
      acc2[2] = gam[1, 2] + log_survival;// survived
      acc2 += log_obs;
        
      gam[2, 2] = log_sum_exp_vec(acc2);
      // change to dead
        // come from alive
      acc3[1] = gam[1, 2] + log1m_exp( log_survival );// dead last year
        
          // come from dead
      acc3[2] = gam[1, 3]; // dead can only be dead
      acc3 += log_obs_nothere; // should not see if dead
        
      gam[2, 3] = log_sum_exp_vec(acc3);
      // we iterate to T + 1, because we inserted a dummy period where 
      // every individual is in the "not recruited" state
      for (t in 3:(Tp1)) {
        log_seen_center_grid = rep_vector(0, n_grid);
        for (j in 1:n_trap) {
          
          log_seen_center_grid += sum(y[i, j, :, t - 1]) * col(dist_tmp, j) + 
                (sum(deploy[j, :, t - 1 ])-sum(y[i, j, :, t - 1])) * 
                log1m_exp(col(dist_tmp, j)); // prob seeing the detection history at a given activaty center
            
        }
        log_obs = log_sum_exp_vec(log_seen_center_grid + log_probs[:,t-1]);// margianlize over activity centers by logsumexp
        log_obs_nothere = sum( to_matrix(y[i, :, :, t - 1])) * negINF;
        // change to un-entered
          // get from un-entered
        acc1 = gam[t - 1, 1] + log1m_exp(log_recruit);
        acc1 += log_obs_nothere;
        gam[t, 1] = acc1;
        
        // change to alive
          // come from un-entered
        acc2[1] = gam[t - 1, 1] + log_recruit;// recruited into
          // come from alive
        acc2[2] = gam[t - 1, 2] + log_survival;// survived
        acc2 += log_obs;
        
        gam[t, 2] = log_sum_exp_vec(acc2);
        // change to dead
          // come from alive
        acc3[1] = gam[t - 1, 2] + log1m_exp( log_survival );// dead last year
        
          // come from dead
        acc3[2] = gam[t - 1, 3]; // dead can only be dead
        acc3 += log_obs_nothere; // should not see if dead
        
        gam[t, 3] = log_sum_exp_vec(acc3);
      }

      log_lik[i] = log_sum_exp_vec(gam[Tp1]);
      //print("log_lik[i]",log_lik[i]);
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
  int<lower=0, upper = 3> z[M, T]; // Latent state of all individuals at all primary occasions
  int<lower = 1, upper = n_grid> s[M, T]; // Activity center for all individuals at all primary occasions 
  { // begin local stuff

    matrix[n_grid,n_trap] log_lik_center_grid; // we need to save all traps, given we are not marginalize over grids any more 
    
    real log_po[M, T];// detection probability, 


    // local stuff for FFBS/Viterbi
    real negINF = -1e20; // local negative inf, to work in log scale avoiding underflow, such a effective negInf is needed to avoid autograd error caused by exp(-Inf) 
    vector[n_grid] log_seen_center_grid;
    matrix[n_grid, n_trap] dist_tmp;
    real acc1;
    vector[2] acc2;
    vector[2] acc3;
    real log_obs;
    real log_obs_nothere; // more a place holder to check if we have seen an individual at a trap when it should NOT be seen
    vector[3] gam[Tp1];
    vector[3] bkw; // backward (unnormalized) probability at a state

    real logp_state; 

    matrix[n_grid, T] log_probs; // log Poisson intensity, at pixels, up to a constant (we do not need it since we condition on numer of total points)
    
    dist_tmp = log_p0 - alpha1 * sq_dist;
    for(t in 1:T){
      log_probs[:, t] = log_softmax( envX[t] * beta_env);// intensity at that year
    }

    dist_tmp = log_p0 - alpha1 * sq_dist;


    // construct emission probablity, as well as activity center for each individual 
    for (i in 1:M) {
      for (t in 1:T) {// loop over primary occasions first
        log_seen_center_grid = rep_vector(0, n_grid);; // reset
      // below calculates probability of seeing an individual at given trap, in a given primary occasion, condition on state of that individual being alive, otherwise one just cannot see the individual. We will assume secondary occasions and individuals are exchangable, thus we will not index by individuals and secondary occasions
        for (j in 1:n_trap) {
          
          log_seen_center_grid += sum(y[i, j, :, t]) * col(dist_tmp, j)+ 
                  (sum(deploy[j, :, t ])-sum(y[i, j, :, t])) * 
                    log1m_exp(col(dist_tmp, j)); // prob seeing the detection history at a given activaty center        
        }
        log_po[i, t] = log_sum_exp_vec(log_seen_center_grid + log_probs[:, t]);// margianlize over activity centers by logsumexp
        // now sum over traps using log_lik_center_grid, we have the probability of seeing the detection history at that primary occasion, given the individual is alive in that primary occasion locating at the certain pixel

        s[i,t] = categorical_rng( softmax( log_seen_center_grid + log_probs[:, t])); // sample activity center for each individual at each primary occasion
      }
    } // end dealing with observation and activaty center

    // FFBS for latent state, see https://github.com/probcomp/metaprob/issues/17
    for (i in 1:M) {
      // All individuals are in state 1 (not recruited) at t=0, we work in log scale
      
      // forward pass
      gam[1, 1] = 0;
      gam[1, 2] = negINF; // not sure if this is a good idea but effectively -inf, exp(-inf) can cause problem in Stan due to autograd
      gam[1, 3] = negINF;


      log_obs = log_po[i, 1];
      log_obs_nothere = sum( to_matrix(y[i, :, :, 1])) * negINF;
        
      // change to un-entered
        // get from un-entered
      acc1 = gam[1, 1] + log1m_exp(log_init_recruit);
      acc1 += log_obs_nothere;
      gam[2, 1] = acc1;


      // change to alive
        // come from un-entered
      acc2[1] = gam[1, 1] + log_init_recruit;// recruited into
          // come from alive
      acc2[2] = gam[1, 2] + log_survival;// survived
      acc2 += log_obs;
      gam[2,2] = log_sum_exp_vec(acc2);

        // change to dead
          // come from alive
      acc3[1] = gam[1, 2] + log1m_exp( log_survival );// dead last year
        // come from dead
      acc3[2] = gam[1, 3]; // dead can only be dead
      acc3 += log_obs_nothere; // should not see if dead
      gam[2,3] = log_sum_exp_vec(acc3);

      // we iterate to T + 1, because we inserted a dummy period where 
      // every individual is in the "not recruited" state
      // forward algorithm, calculate P(z_T, o_1,...,o_T)
      for (t in 3:(Tp1)) {
        // likelihood of seeing/not seeing the individual when it is alive
        
        log_obs = log_po[i, t-1];
        log_obs_nothere = sum( to_matrix(y[i, :, :, t - 1])) * negINF;
        
        // change to un-entered
          // get from un-entered
        acc1 = gam[t - 1, 1] + log1m_exp(log_recruit);
        acc1 += log_obs_nothere;
        gam[t, 1] = acc1;


        // change to alive
          // come from un-entered
        acc2[1] = gam[t - 1, 1] + log_recruit;// recruited into
          // come from alive
        acc2[2] = gam[t - 1, 2] + log_survival;// survived
        acc2 += log_obs;

        gam[t,2] = log_sum_exp_vec(acc2);

        // change to dead
          // come from alive
        acc3[1] = gam[t - 1, 2] + log1m_exp( log_survival );// dead last year
          // come from dead
        acc3[2] = gam[t - 1, 3]; // dead can only be dead
        acc3 += log_obs_nothere; // should not see if dead
        gam[t, 3] = log_sum_exp_vec(acc3);

      }// end forward algorithm

      // backward pass
      // we start at the last period, and work our way back
      
      z[i,T] = categorical_rng(softmax(gam[Tp1,]));
      
      for (tt in 1:(T - 1)) {
        if(z[i,T-tt+1]==1){
          bkw[1] = log1m_exp(log_recruit);
          bkw[2] = negINF;
          bkw[3] = negINF;
        }
        if(z[i,T-tt+1]==2){
          bkw[1] = log_recruit;
          bkw[2] = log_survival;
          bkw[3] = negINF;
        }
        if(z[i,T-tt+1]==3){
          bkw[1] = negINF;
          bkw[2] = log1m_exp(log_survival);
          bkw[3] = 0;
        }
        bkw += gam[Tp1-tt,];
        z[i,T-tt] = categorical_rng(softmax(bkw));
      }
    } // end FFBS
  }// end local stuff
}

