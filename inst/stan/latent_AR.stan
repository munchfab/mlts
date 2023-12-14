// autoregressive DSEM with latent variables
data {
  int<lower=1> N;
  int<lower=1> TP;
  int<lower=1> N_obs;
  array[N] int<lower=1> N_obs_id;
  int<lower=2> N_ind;
  array[N_ind] vector[N_obs] y;
  // 
  // conditionals: 
  int logv_is_random; // constant or random innovation variances? 0 = no, 1 = yes
  // int re_as_outcome; // individual parameters as outcome? 0 = no, 1 = yes
  int n_cov; // number of covariates - minimum of 1 for intercepts
  matrix[N, n_cov] W; // predictors of individual parameters
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
  array[N_ind] int<lower=0> N_miss; // has to be > 2 at the moment or next line fails 
  array[N_ind, max(N_miss)] int pos_miss; // positions of missings
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
  vector[n_random] b[N]; // person-specific parameter
  vector<lower=0>[n_random] sigma; // random effect SD
  cholesky_factor_corr[n_random] L; // cholesky factor of random effects correlation matrix
  vector[sum(N_miss)] y_impute;
  matrix[n_cov, n_random] btw_pred;
  // prediction parameter regressing outcome on individual parameter 
  matrix[n_out, n_out_pred] bs;
  real<lower=0> sigma_out[n_out];
  real alpha_out[n_out];
  // addon for measurement model
  vector[N_ind - 1] lamB;
  vector[N_ind - 1] lamW;
  vector<lower=0>[N_ind - 1] sigmaB;
  vector<lower=0>[N_ind] sigmaW;
  vector[N_ind - 1] item_int;
  array[N_ind - 1] vector[N] yB;
  vector[N_obs] etaW;
}
transformed parameters {
  // population means of person-specific parameters 
  matrix[N, n_random] bmu;
  // vector of innovation variances - in case of random innovation variances 
  // retransform by exponentiation of log innovations sampled from MVN
  vector<lower = 0>[n_inno] sd_noise;
  // vector to store observed data and imputed values in vase of missing obs
  // array[N_ind] vector[N_obs] y_merge;
  // (standardized) individual parameters as predictors 
  matrix[N, n_out_pred] b_pred;
  // calculate population means of person-specific parameters 
  bmu = W * btw_pred;
  if (logv_is_random == 1) {
    sd_noise = sqrt(exp(to_vector(b[, 3])));
  } else {
    sd_noise[1] = sigma[3];
  }
  // prepare predictor matrix if outcome is provided 
  if(n_out > 0) {
    for (nn in 1 : (n_out_pred - n_out_pred_b)) {
      if (out_is_std == 1) {
        b_pred[, nn] = (to_vector(b[, out_pred_which[1, nn]]) - mean(b[, out_pred_which[1,nn]])) / sd(b[, out_pred_which[1,nn]]);
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
  int pos = 1;
  int p_miss = 1;
  matrix[n_random, n_random] SIGMA = diag_pre_multiply(sigma[1 : n_random], L);
  // vector to store observed data and imputed values in vase of missing obs
  array[N_ind] vector[N_obs] y_merge;
  // impute misings only when N_miss > 2
  y_merge = y;
  if (max(N_miss) > 0) {  
    for (i in 1 : N_ind) {
      y_merge[i, pos_miss[i, 1 : N_miss[i]]] = y_impute[p_miss : (p_miss - 1 + N_miss[i])];
      p_miss = p_miss + N_miss[i];
    }
  }
  // priors on elements in the random effect covariance matrix
  L ~ lkj_corr_cholesky(1);
  // priors on random effects variances
  sigma ~ cauchy(0, 10);
  // Priors on fixed effects
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
  // priors on measurement model parameters 
  lamB ~ normal(0, 1);
  lamW ~ normal(0, 1);
  sigmaB ~ cauchy(0, 1);
  sigmaW ~ cauchy(0, 1);
  item_int ~ normal(0, 10);
  // if innovation variances are fixed across subjects 
  // if (logv_is_random == 1) {
  for (pp in 1 : N) {
    array[N_ind] vector[N_obs_id[pp] + 1] yW;
    // individual parameters from multivariate normal distribution 
    b[pp] ~ multi_normal_cholesky(bmu[pp,1 : n_random], SIGMA);
    // decomposing between and within-part  
    for (i in 1 : N_ind) {
      if (i == 1) {
        yW[i, ] = y_merge[i, pos : (pos + (N_obs_id[pp]))] - b[pp, 1];
        yW[i, ] ~ normal(etaW[pos : (pos + (N_obs_id[pp]))], sigmaW[i]);
      } else {
        yB[i - 1, pp] ~ normal(item_int[i - 1] + lamB[i - 1] * b[pp, 1], sigmaB[i - 1]);
        yW[i, ] = y_merge[i, pos : (pos + (N_obs_id[pp]))] - yB[i - 1, pp];
        yW[i, ] ~ normal(lamW[i - 1] * etaW[pos : (pos + (N_obs_id[pp]))], sigmaW[i]);
      }
    }
    // local - why?
    {
    vector[(N_obs_id[pp])] mus = b[pp, 2] * (etaW[pos : (pos + (N_obs_id[pp] - 1))]);
    // mus = b[pp, 2] * (etaW[pos : (pos + (N_obs_id[pp] - 1))]);
    if (logv_is_random == 1) {
      etaW[(pos + 1) : (pos + N_obs_id[pp])] ~ normal(mus, sd_noise[pp]);
    } else {
      etaW[(pos + 1) : (pos + N_obs_id[pp])] ~ normal(mus, sd_noise[1]);
    }
    }
    // update index variables
    pos = pos + N_obs_id[pp] + 1;
    } // end loop over subjects 
  // conditional: prediction of time-invariant outcome 
  if (n_out > 0) {
    for(oo in 1 : n_out) {
      out[oo, ] ~ normal(b_pred * bs[oo, ]' + alpha_out[oo], sigma_out[oo]);
      //out[oo,] ~ normal(b_pred * to_vector(bs[oo,]), sigma_out[oo]);
    }
  }
}
generated quantities{
  matrix[n_random, n_random] bcorr; // random coefficients correlation matrix
  matrix[n_random, n_random] bcov; // random coefficients covariance matrix
  array[n_out] vector[N] yhat;
  real R2_out[n_out]; // explained outcome variance 
  bcorr = multiply_lower_tri_self_transpose(L);
  bcov = quad_form_diag(bcorr, sigma[1 : n_random]);
  // conditional: prediction of time-invariant outcome 
  if (n_out > 0) {
    for(oo in 1 : n_out) {
      yhat[oo, 1 : N] = alpha_out[oo] + b_pred * bs[oo, ]';
      R2_out[oo] = ((sd(out[oo, ])^2)-(sigma_out[oo]^2))/(sd(out[oo, ])^2);
      //R2_out[oo] = (sd(yhat[oo,])^2)/(sd(out[oo,])^2);
    }
  }
}
