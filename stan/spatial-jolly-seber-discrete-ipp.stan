// Jolly-Seber model as a multistate hidden Markov model
// adapted from stan-dev/example-models/blob/master/BPA/Ch.10/js_ms.stan
//--------------------------------------
// States:
// 1 not yet entered
// 2 alive
// 3 dead

data {
  int<lower = 0> M;
  int<lower = 0> T;
  int<lower = 0> Kmax;
  int<lower = 1> n_trap;
  int<lower = 0> K[n_trap,T];
  int<lower = 0, upper = n_trap + 1> y[M, Kmax, T]; // detection history, it is ok to ignore time order etc. as we assume secondary occasions are exchangable
  matrix[n_trap, 2] X; //trap locations
  vector[2] xlim;
  vector[2] ylim;

  //delt with environmental covariates, set up a grid for Poisson intensity approximations
  int<lower = 1> n_grid;
  matrix[n_grid, 2] grid_pts;
  real pixel_area;

  // get enviromental variables at grid points
  int<lower = 1> n_env;
  matrix[n_grid, n_env + 1] envX[T]; // environmental variables at each primary occasion, INCLUDE intercept!
}

transformed data {
  int Tm1 = T - 1;
  int Tp1 = T + 1;

  // deal with environmental covariates
  real log_pixel_area = log(pixel_area);
  matrix[n_grid, n_trap] sq_dist;
  real logM = log(M);
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
  vector<lower = xlim[1], upper = xlim[2]>[T] sx[M]; // x-coordinates of activatiy centers 
  vector<lower = ylim[1], upper = ylim[2]>[T] sy[M]; // y-coordinates of activatiy centers
  //vector<lower = 1, upper = n_grid>[T] s[M]; // activity centers, on grid points
  real p0;
  real<lower = 0> alpha1;

  // deal with environmental covariates
  vector[n_env+1] beta_env; // regression coefficients for activity centers
}

transformed parameters {
  real<lower = 0, upper = 1> c_phi = 1 - phi;
  real log_p0 = log(p0);
  simplex[3] ps[3, T];
  vector[M] log_lik;

  // Define probabilities of state S(t+1) given S(t)
  ps[1, 1, 1] = 1 - psi;
  ps[1, 1, 2] = psi;
  ps[1, 1, 3] = 0;
  ps[2, 1, 1] = 0;
  ps[2, 1, 2] = 1;
  ps[2, 1, 3] = 0;
  ps[3, 1, 1] = 0;
  ps[3, 1, 2] = 0;
  ps[3, 1, 3] = 1;
  
  for (t in 2:T) {
      ps[1, t, 1] = 1 - gamma[t-1];
      ps[1, t, 2] = gamma[t-1];
      ps[1, t, 3] = 0;
      ps[2, t, 1] = 0;
      ps[2, t, 2] = phi;
      ps[2, t, 3] = c_phi;
      ps[3, t, 1] = 0;
      ps[3, t, 2] = 0;
      ps[3, t, 3] = 1;
  }

  { // begin local scope
    vector[2] s;
    vector[n_trap + 1] logits;
    real acc[3];
    vector[3] gam[Tp1];
    vector[n_trap + 1] po[T, 3];

    vector[n_grid] log_mu; // log Poisson intensity, at pixels
    vector[n_grid] log_probs;

    real min_dist = 1e20;
    real sqrdist_tmp;
    int min_dist_grid = 0;

    for (t in 1:T) {// loop over primary occasions first
      log_mu = envX[t] * beta_env + log_pixel_area;// intensity at that year
      log_probs = log_softmax(log_mu);
      for (i in 1:M) {
        s[1] = sx[i, t];
        s[2] = sy[i, t];
        for (j in 1:n_grid) {
          sqrdist_tmp = squared_distance(s, grid_pts[j, ])
          if(sqrdist_tmp < min_dist){
            min_dist = sqrdist_tmp;
            min_dist_grid = i;
          }
        }
        loglik[i] += log_probs[min_dist_grid];// add log prob of being at that pixel
        for (j in 1:n_trap) {
          po[t, 1, j] = 0; // not recruited
          po[t, 3, j] = 0; // dead
          logits[j] = log_p0 - alpha1 * squared_distance(s, X[j, ]);// distance sampling 
        }
        logits[n_trap + 1] = 0;
        po[t, 1, n_trap + 1] = 1;
        po[t, 3, n_trap + 1] = 1;
        po[t, 2, ] = softmax(logits);
      }
    }
    // Forward algorithm
    for (i in 1:M) {
      // All individuals are in state 1 (not recruited) at t=0
      gam[1, 1] = 1;
      gam[1, 2] = 0;
      gam[1, 3] = 0;
  
      // we iterate to T + 1, because we inserted a dummy period where 
      // every individual is in the "not recruited" state
      for (t in 2:(Tp1)) {
        for (k in 1:3) {
          for (j in 1:3) {
            acc[j] = gam[t - 1, j] * ps[j, t - 1, k];
            // the loop below accounts for differences in K across years
            for (occasion in 1:K[t-1]) {
              acc[j] *= po[t-1, k, y[i, occasion, t - 1]];
            }
          }
          gam[t, k] = sum(acc);
        }
      }
      log_lik[i] = log(sum(gam[Tp1]));
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
  array[M, T] int<lower=1, upper=3> z; // Latent state
  
  // Generate z[]
  for (i in 1 : M) {
    z[i, 1] = 1;
    for (t in 2 : Tp1) {
      z[i, t] = categorical_rng(ps[z[i, t - 1], i, t - 1]);
    }
  }
  
}