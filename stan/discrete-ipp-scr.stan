functions {
  matrix log_sum_exp_helper(matrix x, matrix y,int M, int n_grid){
    matrix[M, n_grid] res;
    for(i in 1:M){
      for(j in 1:n_grid){
        res[i,j] = log_sum_exp(x[i,j], y[i,j]);
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
  int<lower = 1> M; // upper bound, number of individuals
  int<lower = 1> n_env; // number of environmental variables
  int<lower = 1> n_trap; // number of traps
  int<lower = 1> n_grid; // number of grids
  int<lower = 1> Kmax;
  matrix[n_grid, 2] grid_pts; // grid locations
  matrix[n_grid, n_env] envX; // environmental variables
  int<lower = 0, upper = 1> y[M, n_trap,Kmax]; // detection history, it is ok to ignore time order etc. as we assume secondary occasions are exchangable
  int<lower = 0, upper = 1> deploy[n_trap,Kmax]; //deployment history of all camera traps at all occasions )(primary and secondary)
  matrix[n_trap, 2] X; //trap locations
}

transformed data {
  matrix[n_grid, n_trap] sq_dist;
  int<lower = 0> y_red[M, n_trap]; // reduced detection history, it is ok to ignore time order etc. as we assume secondary occasions are exchangable
  int<lower = 0> deploy_red[n_trap];
  int<lower = 0> y_red_sum[M]; // number of detections of each individual
  for (j in 1:n_trap) {
    
    deploy_red[j] = sum(deploy[j, :]);
    for(i in 1:M){
      y_red[i,j] = sum(y[i,j,:]);
    }
    // square distance
    for (i in 1:n_grid) {
      sq_dist[i, j] = squared_distance(grid_pts[i, ], X[j, ]);
    }
  }
  for(i in 1:M){
    y_red_sum[i] = sum(y_red[i, :]);
  }
}

parameters {
  // SCR
  real<lower = 0, upper = 1> psi;
  real<lower = 0, upper = 1> p0;
  real<lower = 0> alpha1;

  // environmental things
  vector[n_env] beta_env; // regression coefficients for activity centers
}

transformed parameters {
  
  real<upper = 0> log_p0 = log(p0);
  real<upper = 0> log_psi = log(psi);
  matrix[M, n_grid] log_lik;

  {
    real negINF = -1e20;
    matrix[M, n_grid] log_obs;
    matrix[M, n_grid] log_obs_nothere;

    matrix[n_grid, n_trap] dist_tmp;
    vector[n_grid] log_probs; // log Poisson intensity, at pixels, up to a constant (we do not need it since we condition on numer of total points)
    
    dist_tmp = log_p0 - alpha1 * sq_dist; // distance sampling thing
    log_probs = log_softmax( envX * beta_env);// intensity at that year

    log_obs = rep_matrix(0, M, n_grid);

    // if the individual is alive, the detection probability 
    for (j in 1:n_trap) {
        log_obs += to_vector(y_red[:, j]) * (dist_tmp[:,j])'+ 
          (rep_vector(deploy_red[j],M)-to_vector(y_red[:, j])) * 
          (log1m_exp(dist_tmp[:,j]))'; // prob seeing the detection history at a given activaty center, for a given individual
    }
    
    // if the individual is not alive, the detection probability
    log_obs_nothere = to_vector(y_red_sum) * rep_row_vector(negINF,n_grid);
    log_obs_nothere += log1m_exp(log_psi);
    log_obs += log_psi;

    log_lik = log_sum_exp_helper(log_obs, log_obs_nothere,  M, n_grid) + rep_vector(1,M) * log_probs';// log likelihood of each individual at each grid cell

  }
}




model {
  // priors
  psi ~ beta(1, 1);
  p0 ~ beta(1, 1);
  alpha1 ~ gamma(1, 1);
  beta_env ~ normal(0, 10);
  
  // likelihood
  for(i in 1:M){
    target += log_sum_exp(log_lik[i, :]);
  }
}


generated quantities {
  int<lower=0, upper = 1> z[M]; // Latent state of all individuals at all primary occasions
  int<lower = 1, upper = n_grid> s[M]; // Activity center for all individuals at all primary occasions 

  // local calculations
  {
    matrix[M, n_grid] log_lik_center_grid; // we need to save all traps, given we are not marginalize over grids any more 
    

    
    real negINF = -1e20; // local negative inf, to work in log scale avoiding underflow, such a effective negInf is needed to avoid autograd error caused by exp(-Inf) 
    matrix[n_grid, n_trap] dist_tmp;
    matrix[M, n_grid] log_obs;
    matrix[M, n_grid] log_obs_nothere;
    vector[n_grid] log_probs; // log Poisson intensity, at pixels, up to a constant (we do not need it since we condition on numer of total points)
    

    
    dist_tmp = log_p0 - alpha1 * sq_dist; // distance sampling thing
    log_probs = log_softmax( envX * beta_env);// intensity at that year

    log_obs = rep_matrix(0, M, n_grid);

    // if the individual is alive, the detection probability 
    for (j in 1:n_trap) {
        log_obs += to_vector(y_red[:, j]) * (dist_tmp[:,j])'+ 
          (rep_vector(deploy_red[j],M)-to_vector(y_red[:, j])) * 
          (log1m_exp(dist_tmp[:,j]))'; // prob seeing the detection history at a given activaty center, for a given individual
    }
    
    // if the individual is not alive, the detection probability
    log_obs_nothere = to_vector(y_red_sum) * rep_row_vector(negINF,n_grid);
    
    // use to construct z
    log_obs_nothere = log_obs_nothere + rep_vector(1,M) * log_probs' + log1m_exp(log_psi);
    log_obs = log_obs + rep_vector(1,M) * log_probs' + log_psi;
    for (j in 1:M) {
      z[j] = bernoulli_logit_rng(log_sum_exp(log_obs[j, :])-
                                    log_sum_exp(log_obs_nothere[j, :]));
    }

    // construct s
    log_lik_center_grid = log_sum_exp_helper(log_obs, log_obs_nothere,  M, n_grid);// log likelihood of each individual at each grid cell

    for (j in 1:M) {
      s[j] = categorical_rng(softmax(to_vector(log_lik_center_grid[j, :])));
    }


  }



}
