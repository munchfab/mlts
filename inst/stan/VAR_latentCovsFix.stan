// autoregressive DSEM with manifest variables
data {
  int<lower=1> N; 	        // number of observational units
  int<lower=1> D;           // number of latent constructs
  int<lower=1> D_cen;
  array[D] int<lower=1> D_np;     // number of indicators per construct
  int<lower=1> n_p; 	      // number of manifest indicators
  array[n_p] int D_perP;          // indicate dimension per indicator
  array[n_p] int is_SI;           // indicate if single-indicator per construct
  array[D] int D_pos_is_SI;       // indicate position of single-indicator per construct
  int<lower=1, upper=3> maxLag; // maximum lag
  int<lower=1> N_obs; 	    // observations in total: N * TP
  int<lower=1> n_pars;
  int<lower=D_cen> n_random;    // number of random effects
  int n_fixed;
  array[1,n_fixed] int is_fixed;
  array[n_random] int is_random;  // which parameters to model person-specific
  array[N] int<lower=1> N_obs_id; // number of observations for each unit
  array[n_p] vector[N_obs] y;     // observations as array of vectors

  // handling of missing values
  int n_miss;                        // total number of missings across indicators
  array[n_p] int n_miss_p;                 // missings per indicator
  array[n_p,max(n_miss_p)] int pos_miss_p; // array of missings' positions
  array[n_p] int n_obs_p;
  array[n_p, max(n_obs_p)] int pos_obs_p;

  // censored models
  real censL_val;
  int n_censL;                         // total number of obs at LB across D
  array[n_p] int n_censL_p;                  // obs at LB per D
  array[n_p,max(n_censL_p)] int pos_censL_p; // array of obs at LBs' positions
  real censR_val;
  int n_censR;                         // total number of obs at LB across D
  array[n_p] int n_censR_p;                  // obs at LB per D
  array[n_p,max(n_censR_p)] int pos_censR_p; // array of obs at LBs' positions

  // model adaptations based on user inputs:
  // - fixing parameters to constant values:
  // - innovation variances
  array[D_cen] int<lower=0,upper=1> innos_rand;
  int n_innos_fix;
  array[D_cen] int innos_fix_pos;
  array[D_cen] int innos_pos;

  // - dynamic model specification per D
  array[D] int<lower=0> N_pred;     // Number of predictors per dimension
  array[D,max(N_pred)] int<lower=0> D_pred;    // matrix to index predictors to use per dimension
  array[D,max(N_pred)] int<lower=0> Lag_pred;  // matrix to index lag of used predictors
  array[D] int Dpos1;  // index positions of danymic effect parameters
  array[D] int Dpos2;
  int n_int;
  array[D,max(N_pred)] int D_pred2;    // matrix to index predictors to use per dimension
  array[D,max(N_pred)] int Lag_pred2;  // matrix to index lag of used predictors


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

  // indexing information on constraints
  int n_etaW_free;
  int n_loadBfree;
  int n_loadB_equalW;
  int n_loadWfree;
  int n_alphafree;
  int n_sigmaBfree;
  int n_sigmaWfree;
  array[n_loadBfree] int pos_loadBfree; // positions in relation to the 1:n_p indicators
  array[n_loadB_equalW] int pos_loadB_equalW;
  array[n_loadWfree] int pos_loadWfree;
  array[n_alphafree] int pos_alphafree;
  array[n_sigmaBfree] int pos_sigmaBfree;
  array[n_sigmaWfree] int pos_sigmaWfree;
  // index random manifest indicator means
  int n_YB_free;        // number of indicators for which mu (YB) is not determined by random item mean
  array[n_p] int YB_free_pos; //
  array[n_p] int mu_is_etaB;  //
  array[n_p] int mu_etaB_pos; // indicate whether to use etaB or random item mean
  // get SDs for standardized results
  int<lower=0,upper=1> standardized;

  // priors
  matrix[n_random,2] prior_gamma;
  matrix[n_random,2] prior_sd_R;
  real prior_LKJ;
  matrix[n_innos_fix,2] prior_sigma;
  matrix[n_fixed,2] prior_b_fix;
  matrix[n_cov_bs,2] prior_b_re_pred;
  matrix[n_out,2] prior_alpha_out;
  matrix[n_out_bs_sum,2] prior_b_out;
  matrix[n_out,2] prior_sigma_out;
  matrix[n_alphafree,2] prior_alpha;
  matrix[n_loadBfree,2] prior_loadB;
  matrix[n_loadWfree,2] prior_loadW;
  matrix[n_sigmaBfree,2] prior_sigmaB;
  matrix[n_sigmaWfree,2] prior_sigmaW;

  array[D] int<lower=0,upper=1> is_wcen;
  array[n_p] int<lower=0,upper=1> p_is_wcen;
  array[n_p] int<lower=0,upper=n_p> p_is_wcen_pos;
  array[D] int<lower=0,upper=D> D_cen_pos;
  array[D] int<lower=0,upper=n_p> Dp_cen_pos;
}

transformed data{
  int n_SD_etaW_all;
  int n_SD_etaW_i;
  n_SD_etaW_all = standardized == 1 ? n_etaW_free : 0;
  n_SD_etaW_i = standardized == 1 ? N : 0;
}

parameters {
  array[N] vector[n_random] b_free;            // person-specific parameter
  vector<lower=0>[n_random] sd_R;        // random effect SD
  vector<lower=0>[n_innos_fix] sigma;    // SDs of fixed innovation variances
  cholesky_factor_corr[n_random] L;      // cholesky factor of random effects correlation matrix
  cholesky_factor_corr[D_cen] L_inno;    // cholesky factor of prediction errors
  vector[n_miss] y_impute;               // vector to store imputed values
  vector<upper=censL_val>[n_censL] y_impute_censL;
  vector<lower=censR_val>[n_censR] y_impute_censR;
  row_vector[n_random] gammas;           // fixed effect (intercepts)
  vector[n_cov_bs] b_re_pred;            // regression coefs of RE prediction
  vector[n_fixed] b_fix;
  vector[n_out] alpha_out;               // outcome precition intercepts
  vector<lower=0>[n_out] sigma_out;      // residual SD(s) of outcome(s)
  vector[n_out_bs_sum] b_out_pred;       // regression coefs of out prediction
  // measurement model parameters
  vector[n_loadBfree] loadB_free;
  vector[n_loadWfree] loadW_free;
  vector[n_alphafree] alpha_free;
  vector<lower=0>[n_sigmaBfree] sigmaB_free;
  vector<lower=0>[n_sigmaWfree] sigmaW_free;
  array[n_etaW_free] vector[N_obs] etaW_free;
  array[n_YB_free] vector[N] YB_free;
}

transformed parameters {
  matrix[N, n_random] bmu;     // gammas of person-specific parameters
  matrix[N,n_pars] b;
  array[N] vector[D_cen] sd_noise;
  matrix[n_cov, n_random] b_re_pred_mat = rep_matrix(0, n_cov, n_random);

  vector[n_p] loadB = rep_vector(1, n_p); // measurement model parameters
  vector[n_p] loadW = rep_vector(1, n_p);
  vector[n_p] alpha = rep_vector(0, n_p);
  vector[n_p] sigmaB = rep_vector(0, n_p);
  vector[n_p] sigmaW = rep_vector(0, n_p);

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
  for(i in 1:D_cen){
    if(innos_rand[i] == 0){
      sd_noise[,i] = to_array_1d(rep_vector(sigma[innos_fix_pos[i]],N));
    } else {
      sd_noise[,i] = to_array_1d(sqrt(exp(b[,innos_pos[i]])));
    }
  }

  // replace values for parameters to estimate
  loadW[pos_loadWfree] = loadW_free;
  loadB[pos_loadBfree] = loadB_free;
  loadB[pos_loadB_equalW] = loadW[pos_loadB_equalW];
  alpha[pos_alphafree] = alpha_free;
  sigmaB[pos_sigmaBfree] = sigmaB_free;
  sigmaW[pos_sigmaWfree] = sigmaW_free;
 }

model {
  int pos = 1;       // initialize position indicator
  int p_miss = 1;    // running counter variable to index positions on y_impute
  int p_censL = 1;
  int p_censR = 1;
  int obs_id = 1;    // declare local variable to store variable number of obs per person
  matrix[n_random, n_random] SIGMA = diag_pre_multiply(sd_R, L);
  matrix[D_cen, D_cen] SIGMA_inno = diag_pre_multiply(sd_noise[1,], L_inno);
  array[n_p] vector[N_obs] y_merge;
  array[sum(p_is_wcen)] vector[N_obs] Ymus;
  array[n_p] vector[N] YB;

  y_merge = y;          // add observations
  // add imputed values for missings on each indicator
  for(i in 1:n_p){
    if(n_miss_p[i]>0){
      y_merge[i,pos_miss_p[i,1:n_miss_p[i]]] = segment(y_impute, p_miss, n_miss_p[i]);
      p_miss = p_miss + n_miss_p[i];    // update counter for next indicator i+1
    }
  }

  // replace values at censor thresholds
  for(i in 1:n_p){
    if(n_censL_p[i]>0){
    // add imputed values for observations at floor (threshold for censoring)
    y_merge[i,pos_censL_p[i,1:n_censL_p[i]]] = segment(y_impute_censL, p_censL, n_censL_p[i]);
    p_censL = p_censL + n_censL_p[i];    // update counter for next indicator i+1
    }
    if(n_censR_p[i]>0){
    // add imputed values for observations at ceiling (threshold for censoring)
    y_merge[i,pos_censR_p[i,1:n_censR_p[i]]] = segment(y_impute_censR, p_censR, n_censR_p[i]);
    p_censR = p_censR + n_censR_p[i];    // update counter for next indicator i+1
    }
  }

  // (Hyper-)Priors
  target += normal_lpdf(gammas | prior_gamma[,1],prior_gamma[,2]);
  target += cauchy_lpdf(sd_R | prior_sd_R[,1], prior_sd_R[,2]);
  target += lkj_corr_cholesky_lpdf(L | prior_LKJ);
  target += lkj_corr_cholesky_lpdf(L_inno | prior_LKJ);

  if(n_innos_fix>0){
    target += cauchy_lpdf(sigma | prior_sigma[,1], prior_sigma[,2]);
  }

  if(n_cov > 1){
    target += normal_lpdf(b_re_pred | prior_b_re_pred[,1], prior_b_re_pred[,2]);
  }
  if(n_out > 0){
    target += normal_lpdf(alpha_out | prior_alpha_out[,1], prior_alpha_out[,2]);
    target += normal_lpdf(b_out_pred | prior_b_out[,1], prior_b_out[,2]);
    target += cauchy_lpdf(sigma_out | prior_sigma_out[,1], prior_sigma_out[,2]);
  }

  if(n_fixed > 0){
    target += normal_lpdf(b_fix | prior_b_fix[,1],prior_b_fix[,2]);
  }

  // priors on measurement model parameter
  target += normal_lpdf(alpha_free | prior_alpha[,1], prior_alpha[,2]);
  target += normal_lpdf(loadB_free | prior_loadB[,1], prior_loadB[,2]);
  target += normal_lpdf(loadW_free | prior_loadW[,1], prior_loadW[,2]);
  target += cauchy_lpdf(sigmaB_free | prior_sigmaB[,1], prior_sigmaB[,2]);
  target += cauchy_lpdf(sigmaW_free | prior_sigmaW[,1], prior_sigmaW[,2]);

  // indicator between-part
  for(i in 1:n_p){
    if(p_is_wcen[i] == 0){
      YB[i,] = rep_vector(0,N);
    } else if(mu_is_etaB[i] == 1){
      YB[i,] = b[,mu_etaB_pos[i]];
    } else {
      YB[i,] = YB_free[YB_free_pos[i],];
    }
  }

  for (pp in 1:N) {
    // store number of observations per person
    obs_id = (N_obs_id[pp]);
    int pos_etaW_free = 1;    // running counter variable to index positition on etaW_free
    array[obs_id - maxLag] vector[D_cen] etaW_use;
    array[D] vector[obs_id] etaW_id;
    for(d in 1:D){
      if(is_wcen[d] == 0){
        etaW_id[d,] = segment(y_merge[Dp_cen_pos[d],], pos, obs_id);
      } else if(D_np[d] == 1){
        etaW_id[d,] = segment(y_merge[D_pos_is_SI[d],], pos, obs_id) - YB[D_pos_is_SI[d],pp];
      } else {
        etaW_id[d,] = segment(etaW_free[pos_etaW_free,],pos, obs_id);
        pos_etaW_free = pos_etaW_free + 1;
      }
    }

    // individual parameters from (multivariate) normal distribution
    if(n_random == 1){
      target += normal_lpdf(b_free[pp,1] | bmu[pp,1], sd_R[1]);
    } else {
      target += multi_normal_cholesky_lpdf(b_free[pp, 1:n_random] | bmu[pp, 1 : n_random], SIGMA);
    }

    // dynamic process
    array[obs_id-maxLag] vector[D_cen] mus;

    for(d in 1:D){ // start loop over dimensions

      if(is_wcen[d] == 1){

      // build prediction matrix for specific dimensions
      matrix[(obs_id-maxLag),N_pred[d]] b_mat;  // adjust for non-fully crossed models

      for(nd in 1:N_pred[d]){ // start loop over number of predictors in each dimension
        int lag_use = Lag_pred[d,nd];
        if(D_pred2[d,nd] == -99){
          b_mat[,nd] = etaW_id[D_pred[d, nd],(1+maxLag-lag_use):(obs_id-lag_use)];
        } else {
          int lag_use2 = Lag_pred2[d,nd];
          b_mat[,nd] = etaW_id[D_pred[d, nd],(1+maxLag-lag_use):(obs_id-lag_use)] .*
                       etaW_id[D_pred2[d, nd],(1+maxLag-lag_use2):(obs_id-lag_use2)];
        }
      }

      // use build predictor matrix to calculate latent time-series means
      mus[,D_cen_pos[d]] = to_array_1d(b_mat * to_vector(b[pp, Dpos1[d]:Dpos2[d]]));
      etaW_use[,D_cen_pos[d]] = to_array_1d(segment(etaW_id[d,], (1+maxLag), (obs_id-maxLag)));
      }
    }

    // sampling statement
    if(n_innos_fix == D_cen){
      target += multi_normal_cholesky_lpdf(etaW_use | mus, SIGMA_inno);
    } else {
      target += multi_normal_cholesky_lpdf(etaW_use | mus, diag_pre_multiply(sd_noise[pp,], L_inno));
      }

    // expected indicator scores
    for(i in 1:n_p){
      if(p_is_wcen[i] == 1){
        if(is_SI[i] == 0){
          Ymus[p_is_wcen_pos[i],(pos):(pos-1+obs_id)] = YB[i,pp] + loadW[i] * etaW_id[D_perP[i],];
        } else {
          Ymus[p_is_wcen_pos[i],(pos):(pos-1+obs_id)] = rep_vector(0, obs_id);
        }
      }
    }

    pos = pos + obs_id;
    } // end loop over subjects

    // sampling statements
    for(i in 1:n_p){
      if(p_is_wcen[i] == 1){
        if(mu_is_etaB[i] == 0){
          target += normal_lpdf(YB[i,] | alpha[i] + loadB[i]*b[,mu_etaB_pos[i]], sigmaB[i]);
          }
        if(is_SI[i] == 0){
          target += normal_lpdf(y_merge[i,pos_obs_p[i,1:n_obs_p[i]]] | Ymus[p_is_wcen_pos[i],pos_obs_p[i,1:n_obs_p[i]]], sigmaW[i]);
        }
      }
    }

  // outcome prediction: get expectations of outcome values
  if(n_out > 0){
    int k = 1;
    matrix[N,n_random+n_z] b_z = append_col(b[,is_random],Z);
    for(i in 1:n_out){
      int n_bs = n_out_bs[i,1];      // number of predictors for each outcome
      target += normal_lpdf(out[i,] | alpha_out[i] + b_z[,n_out_b_pos[i,1:n_bs]] * segment(b_out_pred,k,n_bs), sigma_out[i]);
      k = k + n_bs; // update index
    }
  }
}

generated quantities{
  matrix[n_random,n_random] bcorr; // random coefficients correlation matrix
  matrix[D_cen,D_cen] bcorr_inn; // random coefficients correlation matrix
  vector[n_SD_etaW_all] SD_etaW;
  array[n_SD_etaW_all] vector[n_SD_etaW_i] SD_etaW_i;
  bcorr = multiply_lower_tri_self_transpose(L);
  bcorr_inn = multiply_lower_tri_self_transpose(L_inno);
  if(standardized == 1){
    for(i in 1:n_SD_etaW_all){
      SD_etaW[i] = sd(etaW_free[i,]);
      {
      int pos = 1;
      for(p in 1:N){
        SD_etaW_i[i,p] = sd(segment(etaW_free[i,], pos,N_obs_id[p]));
        pos = pos + N_obs_id[p];
        }
      }
    }
  }
}
