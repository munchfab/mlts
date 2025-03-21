// autoregressive DSEM with manifest variables
data {
  int<lower=1> N; 	// number of observational units
  int<lower=1> D; 	// number of time-varying constructs
  int<lower=1, upper=3> maxLag; // maximum lag
  int<lower=1> N_obs; 	// observations in total: N * TP
  int<lower=1> n_pars;
  int<lower=D> n_random;   // number of random effects
  int n_fixed;
  array[1,n_fixed] int is_fixed;
  array[n_random] int is_random;  // which parameters to model person-specific
  array[N] int<lower=1> N_obs_id; // number of observations for each unit
  array[D] vector[N_obs] y; 	    // D*N_obs array of observations

  // handling of missing values
  int n_miss;                      // total number of missings across D
  array[D] int n_miss_D;                 // missings per D
  array[D,max(n_miss_D)] int pos_miss_D; // array of missings' positions

  // model adaptations based on user inputs:
  // - fixing parameters to constant values:
  // - innovation variances
  array[D] int<lower=0,upper=1> innos_rand;
  int n_innos_fix;
  array[D] int innos_fix_pos;
  array[D] int innos_pos;

  // - dynamic model specification per D
  array[D] int<lower=1> N_pred;     // Number of predictors per dimension
  array[D,D*maxLag] int<lower=0> D_pred;    // matrix to index predictors to use per dimension
  array[D,D*maxLag] int<lower=0> Lag_pred;  // matrix to index lag of used predictors
  array[D] int Dpos1;  // index positions of danymic effect parameters
  array[D] int Dpos2;

  // - time-invariant variables:
  // covariates as predictors of random effects
  int<lower=1> n_cov;           // number of covariates - minimum of 1 for intercepts
  int n_cov_bs;
  array[n_cov_bs, 2] int n_cov_mat;
  matrix[N, n_cov] W;  // predictors of individual parameters

  // outcome prediction
  int n_out;                 // number of outcome variables
  array[n_out,1] int n_out_bs;     // number of predictors per outcome
  int n_out_bs_max;          // number of predictors per outcome
  int n_out_bs_sum;          // number of predictors per outcome
  array[n_out,n_out_bs_max] int n_out_b_pos; // index positions
  int n_z;              // number of additional time-invariant as outcome predictors
  matrix[N, n_z] Z;     // observations of Z
  array[n_out] vector[N] out;        // outcome

  // priors
  matrix[n_random,2] prior_gamma;
  matrix[n_random,2] prior_sd_R;
  real prior_LKJ;
  matrix[n_fixed,2] prior_b_fix;
  matrix[n_innos_fix,2] prior_sigma;
  matrix[n_cov_bs,2] prior_b_re_pred;
  matrix[n_out,2] prior_alpha_out;
  matrix[n_out_bs_sum,2] prior_b_out;
  matrix[n_out,2] prior_sigma_out;

  // - covariances of innovations:
  int n_inno_covs; // number of potential innovation covs to include
  int n_obs_cov;   // total number of residuals
  array[1,n_inno_covs] int inno_cov_pos;
  array[2] int<lower=-1,upper=1> inno_cov_load;
}


parameters {
  array[N] vector[n_random] b_free;            // person-specific parameter
  vector<lower=0>[n_random] sd_R;        // random effect SD
  vector<lower=0>[n_innos_fix] sigma;    // SDs of fixed innovation variances
  cholesky_factor_corr[n_random] L;           // cholesky factor of random effects correlation matrix
  vector[n_miss] y_impute;               // vector to store imputed values
  row_vector[n_random] gammas;           // fixed effect (intercepts)
  vector[n_cov_bs] b_re_pred;            // regression coefs of RE prediction
  vector[n_fixed] b_fix;
  vector[n_out] alpha_out;               // outcome precition intercepts
  vector<lower=0>[n_out] sigma_out;      // residual SD(s) of outcome(s)
  vector[n_out_bs_sum] b_out_pred;       // regression coefs of out prediction
  array[n_inno_covs] vector[n_obs_cov] eta_cov;
}

transformed parameters {
  matrix[N, n_random] bmu;     // gammas of person-specific parameters
  matrix[N,n_pars] b;
  array[D] vector<lower = 0>[N] sd_noise;
  array[n_inno_covs] vector<lower = 0>[N] sd_inncov;
  matrix[n_cov, n_random] b_re_pred_mat = rep_matrix(0, n_cov, n_random);

 // REs regressed on covariates
  b_re_pred_mat[1,] = gammas;
  if(n_cov>1){
     for(i in 1:n_cov_bs){
     b_re_pred_mat[n_cov_mat[i,1],n_cov_mat[i,2]] = b_re_pred[i];
    }
  }
  // calculate population means (intercepts) of person-specific parameters
  bmu = W * b_re_pred_mat;

  // create array of (person-specific) parameters to use in model
  for(i in 1:n_random){
    b[,is_random[i]] = to_vector(b_free[,i]);
  }
  if(n_fixed>0){
    for(i in 1:n_fixed){
      b[,is_fixed[1,i]] = rep_vector(b_fix[i],N);
    }
  }

  // transformation of log-innovation variances if modeled as person-specific
  for(i in 1:D){
    if(innos_rand[i] == 0){
      sd_noise[i,] = rep_vector(sigma[innos_fix_pos[i]],N);
    } else {
      sd_noise[i,] = sqrt(exp(b[,innos_pos[i]]));
    }
  }
  // transform log innovation covarainces
  if(n_inno_covs > 0){
      for(i in 1:n_inno_covs){
        sd_inncov[i,1:N] = sqrt(exp(to_vector(b[,inno_cov_pos[1,i]])));
      }
    }
 }

model {
  int pos = 1;       // initialize position indicator
  int p_miss = 1;    // running counter variable to index positions on y_impute
  int pos_cov = 1;   // covariance position
  int obs_id = 1;    // declare local variable to store variable number of obs per person
  array[D] vector[N_obs] y_merge;
  matrix[n_random, n_random] SIGMA = diag_pre_multiply(sd_R, L);

  y_merge = y;      // add observations
  if(n_miss>0){
    for(i in 1:D){
    // add imputed values for missings on each indicator
    y_merge[i,pos_miss_D[i,1:n_miss_D[i]]] = segment(y_impute, p_miss, n_miss_D[i]);
    p_miss = p_miss + n_miss_D[i];    // update counter for next indicator i+1
  }
  }

  // (Hyper-)Priors
  gammas ~ normal(prior_gamma[,1],prior_gamma[,2]);
  sd_R ~ cauchy(prior_sd_R[,1], prior_sd_R[,2]);
  L ~ lkj_corr_cholesky(prior_LKJ);

  if(n_innos_fix>0){
    sigma ~ cauchy(prior_sigma[,1], prior_sigma[,2]);
  }

  if(n_cov > 1){
    b_re_pred ~ normal(prior_b_re_pred[,1], prior_b_re_pred[,2]);
  }
  if(n_out > 0){
    alpha_out ~ normal(prior_alpha_out[,1], prior_alpha_out[,2]);
    b_out_pred ~ normal(prior_b_out[,1], prior_b_out[,2]);
    sigma_out ~ cauchy(prior_sigma_out[,1], prior_sigma_out[,2]);
  }

  if(n_fixed > 0){
    b_fix ~ normal(prior_b_fix[,1],prior_b_fix[,2]);
  }

  for (pp in 1:N) {
    // store number of observations per person
    obs_id = (N_obs_id[pp]);

    // individual parameters from (multivariate) normal distribution
    if(n_random == 1){
      b_free[pp,1] ~ normal(bmu[pp,1], sd_R[1]);
    } else {
      b_free[pp, 1:n_random] ~ multi_normal_cholesky(bmu[pp, 1 : n_random], SIGMA);
    }

    array[n_inno_covs] vector[obs_id-maxLag] eta_cov_id;
    if(n_inno_covs>0){
      for(i in 1:n_inno_covs){
          eta_cov_id[i,] = segment(eta_cov[i,], pos_cov, (obs_id-maxLag));
          eta_cov_id[i,] ~ normal(0, sd_inncov[i,pp]);
      }
    }

    // local variable declaration: array of predicted values
    {
    array[D] vector[obs_id-maxLag] mus;

    // create latent mean centered versions of observations
    array[D] vector[obs_id] y_cen;

    for(d in 1:D){ // start loop over dimensions
      y_cen[d,] = y_merge[d,pos:(pos+obs_id-1)] - b[pp,d];
    }

    for(d in 1:D){ // start loop over dimensions
      // build prediction matrix for specific dimensions
      {
      matrix[(obs_id-maxLag),N_pred[d]+n_inno_covs] b_mat; // adjust for non-fully crossed models

      for(nd in 1:N_pred[d]){ // start loop over number of predictors in each dimension
          int lag_use = Lag_pred[d,nd];
          b_mat[,nd] = y_cen[D_pred[d, nd],(1+maxLag-lag_use):(obs_id-lag_use)];
      }
      if(n_inno_covs>0){
          b_mat[,(N_pred[d]+1)] = eta_cov_id[1,]; // add innovation covariance factor scores
          // use build predictor matrix to calculate latent time-series means
          vector[N_pred[d]+1] b_use;
          b_use[1:N_pred[d]] = to_vector(b[pp, Dpos1[d]:Dpos2[d]]);
          b_use[N_pred[d]+1] = inno_cov_load[d];
          mus[d,] =  b[pp,d] + b_mat * b_use;
        } else{
          mus[d,] =  b[pp,d] + b_mat * to_vector(b[pp, Dpos1[d]:Dpos2[d]]);
        }
      }

      // sampling statement
      y_merge[d,(pos+maxLag):(pos+(obs_id-1))] ~ normal(mus[d,], sd_noise[d,pp]);
      }
     }
    // update index variables
    pos = pos + obs_id;
    pos_cov = pos_cov + obs_id - maxLag;
  } // end loop over subjects

  // outcome prediction: get expectations of outcome values
  if(n_out > 0){
    int k = 1;
    matrix[N,n_random+n_z] b_z = append_col(b[,is_random],Z);
    for(i in 1:n_out){
      int n_bs = n_out_bs[i,1];      // number of predictors for each outcome
      out[i,] ~ normal(alpha_out[i] + b_z[,n_out_b_pos[i,1:n_bs]] * segment(b_out_pred,k,n_bs), sigma_out[i]);
      k = k + n_bs; // update index
    }
  }
}


generated quantities{
  matrix[n_random,n_random] bcorr; // random coefficients correlation matrix
  bcorr = multiply_lower_tri_self_transpose(L);
}
