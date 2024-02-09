// autoregressive DSEM with manifest variables and random interepts
data {
  int<lower=1> N; 	// number of observational units
  int<lower=1> N_obs; 	// observations in total: N * TP
  int<lower=1> N_obs_id[N]; // number of observations for each unit
  vector[N_obs] y; 	// D*N_obs array of observations

  // handling of missing values
  int n_miss;             // total number of missings across D
  int pos_miss[n_miss];   // array of missings' positions

  // - time-invariant variables:
  int n_cov;           // number of covariates - minimum of 1 for intercepts
  int n_cov_bs;
  int n_cov_mat[n_cov_bs, 2];
  matrix[N, n_cov] W;  // predictors of individual parameters

  // priors
  real prior_sd_R_loc;
  real prior_sd_R_scale;
  real prior_gamma_loc;
  real prior_gamma_scale;
  real prior_sigma_loc;
  real prior_sigma_scale;
  matrix[1,n_cov_bs] prior_b_re_pred_loc;
  matrix[1,n_cov_bs] prior_b_re_pred_scale;

}

transformed data{
  int n_cov_bs_use;
  if(n_cov==1){
    n_cov_bs_use = 0;
  } else {
    n_cov_bs_use = n_cov_bs;
  }
}

parameters {
  vector[N] b_free;                      // person-specific parameter
  real sd_R;                             // random effect SD
  vector[n_miss] y_impute;               // vector to store imputed values
  row_vector[1] gammas;                  // fixed effect (intercepts)
  real<lower=0> sigma;
  real ar;
  vector[n_cov_bs_use] b_re_pred;
//  vector[n_fixed] b_fix;
}

transformed parameters {
  matrix[N, 1] bmu;                      // expected values of person-specific parameters
  vector[N] b;

  // regression weights of time-invariant variables as predictors of random effects
  matrix[n_cov, 1] b_re_pred_mat;
  b_re_pred_mat[1,] = gammas;
  if(n_cov>1){
     // start with setting all values to zero
     for(i in 2:n_cov){
       b_re_pred_mat[i,1] = 0;
    }
    // replacement with parameters to estimate
    for(i in 1:n_cov_bs){
    b_re_pred_mat[n_cov_mat[i,1],n_cov_mat[i,2]] = b_re_pred[i];
  }
 }
  // calculate population means (intercepts) of person-specific parameters
  bmu = W * b_re_pred_mat;

  // create array of (person-specific) parameters to use in model
  b = b_free;
  }


model {
  int pos = 1;       // initialize position indicator
  int obs_id = 1;    // declare local variable to store variable number of obs per person
  vector[N_obs] y_merge;

  y_merge = y;      // add observations
  if(n_miss>0){
    // add imputed values for missings on each indicator
    y_merge[pos_miss] = y_impute;
  }

  // priors
  sd_R ~ cauchy(prior_sd_R_loc, prior_sd_R_scale);
  sigma ~ cauchy(prior_sigma_loc, prior_sigma_scale);
  // ar ~ normal(prior_ar_loc, prior_ar_scale);
  ar ~ normal(0, 2);

  // priors on fixed effects
  gammas ~ normal(prior_gamma_loc,prior_gamma_scale);

if(n_cov > 1){
    b_re_pred ~ normal(to_vector(prior_b_re_pred_loc[1,]),to_vector(prior_b_re_pred_scale[1,]));
  }

  // random intercepts
  b_free ~ normal(to_vector(bmu[,1]), sd_R);

  for (pp in 1:N) {
    // store number of observations per person
    obs_id = N_obs_id[pp];

    // local variable declaration: array of predicted values
    {
    vector[obs_id-1] mus;
    // for now:
    // create latent mean centered versions of observations
    vector[obs_id] y_cen;
    y_cen = y_merge[pos:(pos+obs_id-1)] - b[pp];

    // use build predictor matrix to calculate latent time-series means
    mus =  b[pp] + ar * y_cen[1:(obs_id-1)];

    // sampling statement
    y_merge[(pos+1):(pos+(obs_id-1))] ~ normal(mus, sigma);
    }

    // update index variables
    pos = pos + obs_id;
  }

}
