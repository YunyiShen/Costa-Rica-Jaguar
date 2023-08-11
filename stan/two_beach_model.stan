data {
  // dimension things
  int<lower = 1> M; // number of potential individuals

  // detections
  int<lower = 0> y1[M]; // total detections on beach one
  int<lower = 0> y2[M]; // total detections on beach two

  int<lower = 1> deploy_sum[2]; // total number of camera-occasions deployed 

}

parameters{
  simplex[4] psi; // 1. not there, 2. use beach1, 3. use beach2, 4. use both
  real p_mean; // mean logit detection probability
  real <lower = 0> p_std; // variance logit detection 
  vector[M] logit_p; // logit detection probability
}

transformed parameters {
   vector[4] log_psi = log(psi);
}

model {
  psi ~ dirichlet(rep_vector(1, 4));
  p_mean ~ normal(0,10);
  p_std~gamma(1,1);
  logit_p ~ normal(p_mean, p_std);

  // calculate likelihood 
  {
    
    vector[4] log_lik; // log likelihood for each individual, condition on its state
    
    real negINF = -1e20; //effective negInf but without the warning due to autodiff
    for( i in 1:M){
      log_lik[1] = negINF * (y1[i] + y2[i]) + log_psi[1];
      log_lik[2] = binomial_logit_lpmf(y1[i] | deploy_sum[1], logit_p[i]) +
                      negINF * y2[i] + log_psi[2]; 
      log_lik[3] = negINF * y1[i] + 
                      binomial_logit_lpmf(y2[i] | deploy_sum[2], logit_p[i]) + 
                      log_psi[3];
      log_lik[4] = binomial_logit_lpmf(y1[i] | deploy_sum[1], logit_p[i]) +
                      binomial_logit_lpmf(y2[i] | deploy_sum[2], logit_p[i]) + 
                      log_psi[4];
      target += log_sum_exp(log_lik);            
    }
  }
}

generated quantities {
  vector[M] state;
  vector[3] psi_reduced; 
  {
    
    vector[4] log_lik; // log likelihood for each individual, condition on its state
    
    real negINF = -1e20; //effective negInf but without the warning due to autodiff
    psi_reduced = psi[2:4]/(1-psi[1]);
    for(i in 1:M){
      log_lik[1] = negINF * (y1[i] + y2[i]) + log_psi[1];
      log_lik[2] = binomial_logit_lpmf(y1[i] | deploy_sum[1], logit_p[i]) +
                      negINF * y2[i] + log_psi[2]; 
      log_lik[3] = negINF * y1[i] + 
                      binomial_logit_lpmf(y2[i] | deploy_sum[2], logit_p[i]) + 
                      log_psi[3];
      log_lik[4] = binomial_logit_lpmf(y1[i] | deploy_sum[1], logit_p[i]) +
                      binomial_logit_lpmf(y2[i] | deploy_sum[2], logit_p[i]) + 
                      log_psi[4];
      state[i] = categorical_rng( softmax(log_lik));          
    }
  }
}

