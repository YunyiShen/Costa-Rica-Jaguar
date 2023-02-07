// Jolly-Seber model as a multistate hidden Markov model
// adapted from stan-dev/example-models/blob/master/BPA/Ch.10/js_ms.stan
// Jaguars are locked in place through out years
//--------------------------------------
// States:
// 1 not yet entered
// 2 alive
// 3 dead

functions {
  matrix log_sum_exp_helper(matrix[] x, int M, int n_grid){
    matrix[M, n_grid] res;
    for(i in 1:M){
      for(j in 1:n_grid){
        res[i,j] = log_sum_exp(x[:,i,j]);
      }
    }
    return res;
  }

  vector rowSum(matrix x){
    int nrow = rows(x);
    int ncol = cols(x);
    vector[nrow] res;
    res = rep_vector(0, nrow);
    if(nrow<ncol){
      for (i in 1:nrow){
        res[i] = sum(x[i,]);
      }
    }
    else{
      for (i in 1:ncol){
        res += x[,i];
      }
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
  matrix[n_grid, n_env] envX; // environmental variables at each primary occasion, NO NEED OF intercept!
}

transformed data {
  int Tm1 = T - 1;
  int Tp1 = T + 1;

  // deal with environmental covariates
  //real log_pixel_area = log(pixel_area);
  matrix[n_grid, n_trap] sq_dist;
  int<lower = 0> y_red[M, n_trap, T]; // reduced detection history, it is ok to ignore time order etc. as we assume secondary occasions are exchangable
  int<lower = 0> deploy_red[n_trap, T];
  int<lower = 0> y_red_sum[M, T]; // number of detection in that year
  for (j in 1:n_trap) {
    for (k in 1:T) {
      deploy_red[j, k] = sum(deploy[j, :, k]);
      for(i in 1:M){
        y_red[i,j,k] = sum(y[i,j,:,k]);
      }
    }
    for (i in 1:n_grid) {
      sq_dist[i, j] = squared_distance(grid_pts[i, ], X[j, ]);
    }
  }
  for(i in 1:M){
    for (k in 1:T) {
      y_red_sum[i,k] = sum(y_red[i, :, k]);
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
  matrix[M, n_grid] log_lik;

  // Define probabilities of state S(t+1) given S(t)
  log_init_recruit = log(psi);
  log_survival = log(phi);
  
  log_recruit = log(gamma);
      

  { // begin local scope
    real negINF = -1e20; // local negative inf, to work in log scale avoiding underflow, such a effective negInf is needed to avoid autograd error caused by exp(-Inf) 
    //real log_po[M, T]; // detection probability, given the state is alive, otherwise is -inf (emission in HMM sense), details later
    matrix[M, n_grid] acc1;
    matrix[M, n_grid] acc2[2];
    matrix[M, n_grid] acc3[2];
    matrix[M, n_grid] log_obs;
    matrix[M, n_grid] log_obs_nothere; // more a place holder to check if we have seen an individual at a trap when it should NOT be seen
    matrix[M, n_grid] gam[Tp1,3];
    int flag = 0;

    
    matrix[n_grid, n_trap] dist_tmp;
    vector[n_grid] log_probs; // log Poisson intensity, at pixels, up to a constant (we do not need it since we condition on numer of total points)
    
    dist_tmp = log_p0 - alpha1 * sq_dist;
    log_probs = log_softmax( envX * beta_env);// intensity at that year

    

    
    // Forward algorithm
    
      /*
      TODO: vectorize the loop over grids. By
      1) keep tracking of gamma for each grid, i.e. make gam[T, 3, n_grid]
      2) keep tracking of log_obs of each grid, get log_obs a vector of n_grid
      3) make calculation of log_obs vectorized
      
      TODO: consider if we can even vectorize the loop over individuals, given they are considered independent
      DONE, testing
      */
      // All individuals are in state 1 (not recruited) at t=0, we work in log scale
      gam[1, 1] = rep_matrix(0, M, n_grid);
      gam[1, 2] = rep_matrix(negINF, M, n_grid); // not sure if this is a good idea but effectively -inf, exp(-inf) can cause problem in Stan due to autograd
      gam[1, 3] = rep_matrix(negINF, M, n_grid);
      // manually set the first period's recruitment
      log_obs = rep_matrix(0, M, n_grid);
      for (j in 1:n_trap) {
          log_obs += to_vector(y_red[:, j, 1]) * (dist_tmp[:,j])'+ 
            (rep_vector(deploy_red[j, 1],M)-to_vector(y_red[:, j, 1])) * 
            (log1m_exp(dist_tmp[:,j]))'; // prob seeing the detection history at a given activaty center, for a given individual
      }
        
      log_obs_nothere = to_vector(y_red_sum[:, 1]) * rep_row_vector(negINF,n_grid);
      // change to un-entered
        // get from un-entered
      acc1 = gam[1, 1];
      acc1 += (log1m_exp(log_init_recruit) + log_obs_nothere);
      gam[2, 1] = acc1;
        
      // change to alive
        // come from un-entered
      acc2[1] = gam[1, 1];// recruited into
      acc2[1] += log_init_recruit;
        // come from alive
      acc2[2] = gam[1, 2];// survived
      acc2[2] += (log_obs + log_survival);
          
      gam[2, 2] = log_sum_exp_helper(acc2,M, n_grid);
      // change to dead
        // come from alive
      acc3[1] = gam[1, 2];// dead last year
      acc3[1] += log1m_exp( log_survival );
      acc3[1] += log_obs_nothere; // should not see if dead
        
        // come from dead
      acc3[2] = gam[1, 3]; // dead can only be dead
      acc3[2] += log_obs_nothere; // should not see if dead
          
      gam[2, 3] = log_sum_exp_helper(acc3,M, n_grid);
      // we iterate to T + 1, because we inserted a dummy period where 
      // every individual is in the "not recruited" state
      for (t in 3:(Tp1)) {
        log_obs = rep_matrix(0, M, n_grid);
          for (j in 1:n_trap) {
            log_obs += to_vector(y_red[:, j, t - 1]) * (dist_tmp[:,j])'+ 
              (rep_vector(deploy_red[j, t - 1],M)-to_vector(y_red[:, j, t - 1])) * 
              (log1m_exp(dist_tmp[:,j]))'; // prob seeing the detection history at a given activaty center
          }
        log_obs_nothere = to_vector(y_red_sum[:, t - 1]) * rep_row_vector(negINF,n_grid);
        // change to un-entered
          // get from un-entered
        acc1 = gam[t - 1, 1] ;
        acc1 += (log1m_exp(log_recruit) + log_obs_nothere);
        gam[t, 1] = acc1;
          
        // change to alive
          // come from un-entered
        acc2[1] = gam[t - 1, 1];// recruited into
        acc2[1] += log_recruit;
          // come from alive
        acc2[2] = gam[t - 1, 2] ;// survived
        acc2[2] += (log_obs + log_survival);
          
        gam[t, 2] = log_sum_exp_helper(acc2,M, n_grid);
        // change to dead
          // come from alive
        acc3[1] = gam[t - 1, 2];// dead last year
        acc3[1]+= log1m_exp( log_survival );
        acc3[1] += log_obs_nothere; // should not see if dead
          // come from dead
        acc3[2] = gam[t - 1, 3]; // dead can only be dead
        acc3[2] += log_obs_nothere; // should not see if dead
        
        gam[t, 3] = log_sum_exp_helper(acc3,M, n_grid);
      } // end time
	     
      log_lik = log_sum_exp_helper(gam[Tp1],M, n_grid) + rep_vector(1,M) * log_probs'; // prior on center location
  } // end local scope
} 

model {
  p0 ~ normal(0, 3);
  alpha1 ~ lognormal(0, 1);
  beta_env ~ normal(0, 10);
  for(i in 1:M){
    target += log_sum_exp(log_lik[i, :]);
  }
  
}

// note that we could sample z and N in the generated quantities block


generated quantities {
  int<lower=0, upper = 3> z[M, T]; // Latent state of all individuals at all primary occasions
  int<lower = 1, upper = n_grid> s[M]; // Activity center for all individuals at all primary occasions 
  { // begin local stuff

    matrix[n_grid,n_trap] log_lik_center_grid; // we need to save all traps, given we are not marginalize over grids any more 
    
    real log_po[M, T];// detection probability, 


    // local stuff for FFBS
    real negINF = -1e20; // local negative inf, to work in log scale avoiding underflow, such a effective negInf is needed to avoid autograd error caused by exp(-Inf) 
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
    

    // FFBS for latent state, see https://github.com/probcomp/metaprob/issues/17
    for (i in 1:M) {
      s[i] = categorical_rng( softmax((log_lik[i,:])')); // sample activity center for each individual 
      // All individuals are in state 1 (not recruited) at t=0, we work in log scale
      
      // forward pass
      gam[1, 1] = 0;
      gam[1, 2] = negINF; // not sure if this is a good idea but effectively -inf, exp(-inf) can cause problem in Stan due to autograd
      gam[1, 3] = negINF;
      
      log_obs = 0;
      for (j in 1:n_trap) {
          
        log_obs += sum(y[i, j, :, 1]) * dist_tmp[s[i],j]+ 
                  (sum(deploy[j, :, 1 ])-sum(y[i, j, :, 1])) * 
                    log1m_exp(dist_tmp[s[i],j]); // prob seeing the detection history at a given activaty center        
      }
      
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
      gam[2,2] = log_sum_exp(acc2);

        // change to dead
          // come from alive
      acc3[1] = gam[1, 2] + log1m_exp( log_survival );// dead last year
        // come from dead
      acc3[2] = gam[1, 3]; // dead can only be dead
      acc3 += log_obs_nothere; // should not see if dead
      gam[2,3] = log_sum_exp(acc3);

      // we iterate to T + 1, because we inserted a dummy period where 
      // every individual is in the "not recruited" state
      // forward algorithm, calculate P(z_T, o_1,...,o_T)
      for (t in 3:(Tp1)) {
        // likelihood of seeing/not seeing the individual when it is alive
        log_obs = 0;
      	for (j in 1:n_trap) {
          
          log_obs += sum(y[i, j, :, t-1]) * dist_tmp[s[i],j]+ 
                  (sum(deploy[j, :, t-1 ])-sum(y[i, j, :, t-1])) * 
                    log1m_exp(dist_tmp[s[i],j]); // prob seeing the detection history at a given activaty center        
        }
        
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

        gam[t,2] = log_sum_exp(acc2);

        // change to dead
          // come from alive
        acc3[1] = gam[t - 1, 2] + log1m_exp( log_survival );// dead last year
          // come from dead
        acc3[2] = gam[t - 1, 3]; // dead can only be dead
        acc3 += log_obs_nothere; // should not see if dead
        gam[t, 3] = log_sum_exp(acc3);

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

