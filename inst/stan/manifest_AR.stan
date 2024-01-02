// autoregressive DSEM with manifest variables
data {
  int<lower=1> N; // number of observational units
  int<lower=1> TP; // number of time points
  int<lower=1> N_obs; // observations in total: N * TP
  array[N] int<lower=1> N_obs_id; // number of observations for each unit
  vector[N_obs] y; // vector of data points
  // 
  // conditionals:
  int<lower=0, upper=1> logv_is_random; // constant or random innovation variances? 0 = no, 1 = yes
  // int<lower=0, upper=1> re_as_outcome; individual parameters as outcome? 0 = no, 1 = yes
  int n_cov; // number of covariates - minimum of 1 for intercepts
  matrix[N, n_cov] W;  // predictors of individual parameters
  // individual effect as predicors of time-invariant outcome
  int n_out; // number of outcomes 
  int n_out_pred; // number of predictors (including intercepts when out_is_std = 0)
  int n_out_pred_b; // number of additional person parameters to use as predictors 
  array[1, n_out_pred - n_out_pred_b] int out_pred_which; // index of RE to use as predictor (1 = MU, 2 = AR, 3 = LOGV)
  array[n_out] vector[N] out;
  matrix[N, n_out_pred_b] out_pred_b;
  int<lower=0, upper=1> out_is_std; // were outcome variables standardized? 0 = no, 1 = yes
  // 
  // addon for missing 
  int<lower=0> N_miss; // has to be > 2 at the moment or next line fails 
  array[1, N_miss] int pos_miss; // positions of missings
} 
transformed data{
  int n_random; // number of random effects  
  int n_inno; // number of values for innovation variances 
  // select the number of random parameters depending on whether logv is constant
  if (logv_is_random == 0) {
    n_random = 2;
    n_inno = 1;
    } else {
      n_random = 3;
      n_inno = N;
    }
   // standardize outcome
 }
parameters {
  array[N] vector[n_random] b; // person-specific parameter
  vector<lower=0>[n_random] sigma; // random effect SD
  cholesky_factor_corr[n_random] L; // cholesky factor of random effects correlation matrix
  vector[N_miss] y_impute;
  matrix[n_cov, n_random] btw_pred;
  // prediction parameter regressing outcome on individual parameter 
  matrix[n_out, n_out_pred] bs;
  real<lower=0> sigma_out[n_out];
  real alpha_out[n_out];
} 
transformed parameters {
  // population means of person-specific parameters 
  matrix[N, n_random] bmu;                 
  // vector of innovation variances - in case of random innovation variances 
  // retransform by exponentiation of log innovations sampled from MVN
  vector<lower = 0>[n_inno] sd_noise;
  // vector to store observed data and imputed values in vase of missing obs
  vector[N_obs] y_merge;
  // (standardized) individual parameters as predictors
  matrix[N,n_out_pred] b_pred;
  // calculate population means of person-specific parameters
  bmu = W * btw_pred;
  if(logv_is_random == 1) {
    sd_noise = sqrt(exp(to_vector(b[, 3])));
  } else {
    sd_noise[1] = sigma[3];
  }
  // impute misings only when N_miss > 2
  y_merge = y;
  if (N_miss > 0) {  
    y_merge[pos_miss[1, ]] = y_impute; // include missings at relevant position
  }
  // prepare predictor matrix if outcome is provided 
  if(n_out > 0) {
    for (nn in 1 : (n_out_pred - n_out_pred_b)) {
      if (out_is_std == 1) {
        b_pred[, nn] = (to_vector(b[, out_pred_which[1,nn]]) - mean(b[, out_pred_which[1, nn]])) / sd(b[, out_pred_which[1, nn]]);
      } else {
        b_pred[, nn] = to_vector(b[, out_pred_which[1, nn]]);
      }
    }
    for (nn in 1 : n_out_pred_b) {
      b_pred[, nn + n_out_pred - n_out_pred_b] = out_pred_b[, nn];
    }
  }
}
model {
  int pos = 1; // initialize position indicator
  matrix[n_random, n_random] SIGMA = diag_pre_multiply(sigma[1 : n_random], L);
  // priors on elements in the random effect covariance matrix
  L ~ lkj_corr_cholesky(1);  // maybe include option to specify priors in fit function?
  // priors on random effects variances
  sigma ~ cauchy(0, 10);
  // priors on fixed effects
  for (kk in 1 : n_cov) {
    for (nn in 1 : n_random) {
      if (nn == 1) {
        btw_pred[kk, nn] ~ normal(0, 10);
      } else {
      btw_pred[kk, nn] ~ normal(0, 1);
      }
    }
  }
  // priors on prediction parameter 
  if (n_out > 0) {
    if (out_is_std) {
      to_vector(bs) ~ normal(0, 5);
      sigma_out ~ cauchy(0, 1);
      alpha_out ~ normal(0, 1);
    } else {
      to_vector(bs) ~ normal(0, 10);
      sigma_out ~ cauchy(0, 10);
      alpha_out ~ normal(0, 100);
    }
  }
  for (pp in 1:N) {
    // individual parameters from multivariate normal distribution
    b[pp] ~ multi_normal_cholesky(bmu[pp, 1 : n_random], SIGMA); 
    vector[(N_obs_id[pp])] mus;
    mus = b[pp, 1] + b[pp, 2] * (y_merge[pos : (pos + (N_obs_id[pp] - 1))] - b[pp, 1]);
    // if innovation variances are fixed across subjects 
    if (logv_is_random == 1) {
      y_merge[(pos + 1) : (pos + N_obs_id[pp])] ~ normal(mus, sd_noise[pp]);
    } else {
      y_merge[(pos + 1):(pos + N_obs_id[pp])] ~ normal(mus, sd_noise[1]);
    }
    // update index variables
    pos = pos + N_obs_id[pp] + 1; 
    }
  // conditional: prediction of time-invariant outcome 
  if (n_out > 0) {
    for (oo in 1:n_out) {
    out[oo,] ~ normal(b_pred * bs[oo,]' + alpha_out[oo], sigma_out[oo]);
    //out[oo,] ~ normal(b_pred * to_vector(bs[oo,]), sigma_out[oo]);
    }
  }
}
generated quantities{
  matrix[n_random, n_random] bcorr; // random coefficients correlation matrix
  matrix[n_random, n_random] bcov; // random coefficients covariance matrix
  array[n_out] vector[N] y_hat; // predictions
  real R2_out[n_out]; // explained outcome variance 
  // create random coefficient matrices
  bcorr = multiply_lower_tri_self_transpose(L);
  bcov = quad_form_diag(bcorr, sigma[1:n_random]);
  // conditional: prediction of time-invariant outcome
  if (n_out > 0) {
    for (oo in 1 : n_out) {
    y_hat[oo, 1 : N] = alpha_out[oo] + b_pred * bs[oo, ]';
    R2_out[oo] = ((sd(out[oo,])^2) - (sigma_out[oo]^2)) / (sd(out[oo, ])^2);
    // R2_out[oo] = (sd(y_hat[oo,])^2)/(sd(out[oo,])^2);
    }
  }
  // posterior predictive checks for y
  int pos = 1;
  array[N_obs] real y_rep;
  for (pp in 1:N) {
    if (logv_is_random == 1) {
      y_rep[pos] = b[pp, 1];
      y_rep[(pos + 1):(pos + N_obs_id[pp])] = normal_rng(b[pp, 1] + b[pp, 2] * (y_merge[pos : (pos + (N_obs_id[pp] - 1))] - b[pp, 1]), sd_noise[pp]);
    } else {
      y_rep[pos] = b[pp, 1];
      y_rep[(pos + 1):(pos + N_obs_id[pp])] = normal_rng(b[pp, 1] + b[pp, 2] * (y_merge[pos : (pos + (N_obs_id[pp] - 1))] - b[pp, 1]), sd_noise[1]);
    }
    pos = pos + N_obs_id[pp] + 1; 
  }
}
