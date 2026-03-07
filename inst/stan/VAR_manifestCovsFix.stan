// autoregressive DSEM with manifest variables
data {
  int<lower=1> N; 	// number of observational units
  int<lower=1> G;   // number of groups
  int<lower=1> D; 	// number of time-varying constructs
  int<lower=1> D_cen;
  int<lower=0, upper=3> maxLag; // maximum lag
  int<lower=1> N_obs; 	// observations in total: N * TP
  int<lower=1> n_pars;
  int<lower=1> n_random;   // number of random effects

  // new
  int<lower=0> n_mvn1;
  int<lower=0> n_mvn2;
  int<lower=0> n_iid;
  array[n_iid] int pos_iid;
  array[n_mvn1] int pos_mvn1;
  array[n_mvn2] int pos_mvn2;

  int n_fixed;
  array[1,n_fixed] int is_fixed;
  array[n_random] int is_random;  // which parameters to model person-specific
  array[N] int<lower=1> N_obs_id; // number of observations for each unit
  array[D] vector[N_obs] y; 	    // D*N_obs array of observations

  // handling of missing values
  int n_miss;                      // total number of missings across D
  array[D] int n_miss_D;                 // missings per D
  array[D,max(n_miss_D)] int pos_miss_D; // array of missings' positions

  // censored models
  real censL_val;
  int n_censL;                      // total number of obs at LB across D
  array[D] int n_censL_D;                 // obs at LB per D
  array[D,max(n_censL_D)] int pos_censL_D; // array of obs at LBs' positions
  real censR_val;
  int n_censR;                      // total number of obs at LB across D
  array[D] int n_censR_D;                 // obs at LB per D
  array[D,max(n_censR_D)] int pos_censR_D; // array of obs at LBs' positions

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
  array[D,max(N_pred)] int D_pred2;
  array[D,max(N_pred)] int Lag_pred2;

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
  matrix[n_innos_fix,2] prior_sigma;
  matrix[n_fixed,2] prior_b_fix;
  matrix[n_cov_bs,2] prior_b_re_pred;
  matrix[n_out,2] prior_alpha_out;
  matrix[n_out_bs_sum,2] prior_b_out;
  matrix[n_out,2] prior_sigma_out;

  // group-specific?
  array[G] int N_G;               // number of clusters by group
  array[N] int g_id;              // index group per cluster
  array[G,max(N_G)] int g_id_pos; // index group by G*N_G array

  // - covariances of innovations:
  // int<lower=1> n_inno_covs; // number of potential innovation covs to include
  // int<lower=0,upper=1> inno_cov0;    // fixed to zero
  // int<lower=0,upper=1> inno_cov_fix; // fixed to zero
  // addon for time-varying exogenous variables
  array[D] int<lower=0,upper=1> is_wcen;
  array[D] int<lower=0,upper=D> D_cen_pos;

  array[D] int<lower=0,upper=1> is_rdsem;
  array[D] int<lower=0> N_pred_rdsem;
  array[D,max(N_pred_rdsem)] int<lower=0> D_pred_rdsem;
  array[D] int Dpos1_rdsem;
}


parameters {
  array[N] vector[n_random] b_free;            // person-specific parameter
  array[G] vector<lower=0>[n_random] sd_R; // random effect SD
  array[G] vector<lower=0>[n_innos_fix] sigma;    // SDs of fixed innovation variances
  array[G] cholesky_factor_corr[n_mvn1] L;      // cholesky factor of random effects correlation matrix
  array[G] cholesky_factor_corr[n_mvn2] L2;     // cholesky factor of random effects correlation matrix
  array[G] cholesky_factor_corr[D_cen] L_inno;    // cholesky factor of prediction errors
  array[G] row_vector[n_random] gammas;           // fixed effect (intercepts)
  array[G] vector[n_cov_bs] b_re_pred;            // regression coefs of RE prediction
  array[G] vector[n_fixed] b_fix;
  array[G] vector[n_out] alpha_out;               // outcome precition intercepts
  array[G] vector<lower=0>[n_out] sigma_out;      // residual SD(s) of outcome(s)
  array[G] vector[n_out_bs_sum] b_out_pred;       // regression coefs of out prediction
  vector[n_miss] y_impute;               // vector to store imputed values
  vector<upper=censL_val>[n_censL] y_impute_censL;
  vector<lower=censR_val>[n_censR] y_impute_censR;
}

transformed parameters {
  matrix[N, n_random] bmu;     // gammas of person-specific parameters
  matrix[N,n_pars] b;
  array[N] vector[D_cen] sd_noise;
  array[G] matrix[n_cov, n_random] b_re_pred_mat;

 // REs regressed on covariates
  for(g in 1:G){
    b_re_pred_mat[g] = rep_matrix(0, n_cov, n_random);
    b_re_pred_mat[g,1,] = gammas[g,];
    if(n_cov>1){
      for(i in 1:n_cov_bs){
      b_re_pred_mat[g,n_cov_mat[i,1],n_cov_mat[i,2]] = b_re_pred[g,i];
      }
    }
    // calculate population means (intercepts) of person-specific parameters
    bmu[g_id_pos[g,1:N_G[g]],] = W[g_id_pos[g,1:N_G[g]],] * b_re_pred_mat[g];
  }

  // create array of (person-specific) parameters to use in model
  for(i in 1:n_random){
    b[,is_random[i]] = to_vector(b_free[,i]);
  }
  if(n_fixed>0){
    for(i in 1:n_fixed){
      for(g in 1:G){
        b[g_id_pos[g,1:N_G[g]],is_fixed[1,i]] = rep_vector(b_fix[g,i],N_G[g]);
      }
    }
  }

  // transformation of log-innovation variances if modeled as person-specific
  for(i in 1:D_cen){
    if(innos_rand[i] == 0){
      for(g in 1:G){
        sd_noise[g_id_pos[g,1:N_G[g]],i] = to_array_1d(rep_vector(sigma[g,innos_fix_pos[i]],N_G[g]));
        }
    } else {
      sd_noise[,i] = to_array_1d(sqrt(exp(b[,innos_pos[i]])));
    }
  }
 }

model {
  int pos = 1;       // initialize position indicator
  int p_miss = 1;    // running counter variable to index positions on y_impute
  int p_censL = 1;
  int p_censR = 1;
  int obs_id = 1;    // declare local variable to store variable number of obs per person
  array[D] vector[N_obs] y_merge;
  array[G] matrix[n_mvn1, n_mvn1] SIGMA;
  array[G] matrix[n_mvn2, n_mvn2] SIGMA2;

  for(g in 1:G){
    if(n_mvn1 > 0){
      SIGMA[g] = diag_pre_multiply(sd_R[g,pos_mvn1], L[g,]);
    }
    if(n_mvn2 > 0){
      SIGMA2[g] = diag_pre_multiply(sd_R[g,pos_mvn2], L2[g,]);
    }
  }

  y_merge = y;      // add observations

  for(i in 1:D){
    if(n_miss_D[i]>0){
    // add imputed values for missings on each indicator
    y_merge[i,pos_miss_D[i,1:n_miss_D[i]]] = segment(y_impute, p_miss, n_miss_D[i]);
    p_miss = p_miss + n_miss_D[i];    // update counter for next indicator i+1
    }
  }
  // replace values at censor thresholds
  for(i in 1:D){
    if(n_censL_D[i]>0){
    // add imputed values for observations at floor (threshold for censoring)
    y_merge[i,pos_censL_D[i,1:n_censL_D[i]]] = segment(y_impute_censL, p_censL, n_censL_D[i]);
    p_censL = p_censL + n_censL_D[i];    // update counter for next indicator i+1
    }
    if(n_censR_D[i]>0){
    // add imputed values for observations at ceiling (threshold for censoring)
    y_merge[i,pos_censR_D[i,1:n_censR_D[i]]] = segment(y_impute_censR, p_censR, n_censR_D[i]);
    p_censR = p_censR + n_censR_D[i];    // update counter for next indicator i+1
    }
  }

  // (Hyper-)Priors
  for(g in 1:G){
    target += normal_lpdf(gammas[g,] | prior_gamma[,1],prior_gamma[,2]);
    target += cauchy_lpdf(sd_R[g,] | prior_sd_R[,1], prior_sd_R[,2]);
    if(n_mvn1>0){
      target += lkj_corr_cholesky_lpdf(L[g,] | prior_LKJ);
    }
    if(n_mvn2>0){
      target += lkj_corr_cholesky_lpdf(L2[g,] | prior_LKJ);
    }
    target += lkj_corr_cholesky_lpdf(L_inno[g,] | prior_LKJ);

    if(n_innos_fix>0){
      target += cauchy_lpdf(sigma[g,] | prior_sigma[,1], prior_sigma[,2]);
    }

    if(n_cov > 1){
      target += normal_lpdf(b_re_pred[g,] | prior_b_re_pred[,1], prior_b_re_pred[,2]);
    }
    if(n_out > 0){
      target += normal_lpdf(alpha_out[g,] | prior_alpha_out[,1], prior_alpha_out[,2]);
      target += normal_lpdf(b_out_pred[g,] | prior_b_out[,1], prior_b_out[,2]);
      target += cauchy_lpdf(sigma_out[g,] | prior_sigma_out[,1], prior_sigma_out[,2]);
    }

    if(n_fixed > 0){
      target += normal_lpdf(b_fix[g,] | prior_b_fix[,1],prior_b_fix[,2]);
    }
  }

  for (pp in 1:N) {
    // store number of observations per person
    obs_id = (N_obs_id[pp]);
    // obtain group
    int pp_g = g_id[pp];

    array[obs_id - maxLag] vector[D_cen] y_use;

  // individual parameters from (multivariate) normal distribution
    if(n_iid > 0){
      for(jj in pos_iid){
         target += normal_lpdf(b_free[pp,jj] | bmu[pp,jj], sd_R[pp_g, jj]);
      }
    }
    if(n_mvn1 > 0) {
      target += multi_normal_cholesky_lpdf(b_free[pp, pos_mvn1] | bmu[pp, pos_mvn1], SIGMA[pp_g]);
    }
    if(n_mvn2 > 0) {
      target += multi_normal_cholesky_lpdf(b_free[pp, pos_mvn2] | bmu[pp, pos_mvn2], SIGMA2[pp_g]);
    }

    // local variable declaration: array of predicted values
    array[obs_id-maxLag] vector[D_cen] mus;

    // create latent mean centered versions of observations
    array[D] vector[obs_id] y_cen;

    for(d in 1:D){ // start loop over dimensions
      if(is_wcen[d] == 1){
        y_cen[d,] = y_merge[d,pos:(pos+obs_id-1)] - b[pp,D_cen_pos[d]];
      } else {
        y_cen[d,] = y_merge[d,pos:(pos+obs_id-1)];
      }
    }

    if( sum(is_rdsem) > 0 ){
      for (d in 1:D) {
        if (is_rdsem[d] == 1 ) {
          for(k in 1:N_pred_rdsem[d]){
             y_cen[d,] = y_cen[d,] - b[pp,(Dpos1_rdsem[d]+(k-1))] * y_cen[D_pred_rdsem[d,k],];
          }
        }
      }
    }

    for(d in 1:D){ // start loop over dimensions

      if(is_wcen[d] == 1){
        if(N_pred[d] > 0){
          // build prediction matrix for specific dimensions
          matrix[(obs_id-maxLag),N_pred[d]] b_mat; // adjust for non-fully crossed models

          for(nd in 1:N_pred[d]){ // start loop over number of predictors in each dimension
              int lag_use = Lag_pred[d,nd];
              if(D_pred2[d,nd] == -99){
                  b_mat[,nd] = y_cen[D_pred[d, nd],(1+maxLag-lag_use):(obs_id-lag_use)];
              } else {
                int lag_use2 = Lag_pred2[d,nd];
                b_mat[,nd] = y_cen[D_pred [d,nd],(1+maxLag-lag_use ):(obs_id-lag_use)] .*
                             y_cen[D_pred2[d,nd],(1+maxLag-lag_use2):(obs_id-lag_use2)];
                }
              }
          mus[,D_cen_pos[d]] =  to_array_1d(b_mat * to_vector(b[pp, Dpos1[d]:Dpos2[d]]));
        } else {
          mus[,D_cen_pos[d]] =  to_array_1d(rep_vector(0,(obs_id-maxLag)));
        }

      y_use[,D_cen_pos[d]] = to_array_1d(segment(y_cen[d,], (1+maxLag), (obs_id-maxLag)));
      }
    }

    // sampling statement
    if(n_innos_fix == D_cen){
      target += multi_normal_cholesky_lpdf(y_use | mus, diag_pre_multiply(sd_noise[pp,], L_inno[pp_g,]));
    } else {
      target += multi_normal_cholesky_lpdf(y_use | mus, diag_pre_multiply(sd_noise[pp,], L_inno[pp_g,]));
    }

    // update index variables
    pos = pos + obs_id;
  }

  // outcome prediction: get expectations of outcome values
  if(n_out > 0){
    for(g in 1:G){
      int k = 1;
      matrix[N_G[g],n_random+n_z] b_z = append_col(b[g_id_pos[g,1:N_G[g]],],Z[g_id_pos[g,1:N_G[g]],]);
      for(i in 1:n_out){
        int n_bs = n_out_bs[i,1];      // number of predictors for each outcome
        target += normal_lpdf(out[i,g_id_pos[g,1:N_G[g]]] | alpha_out[g,i] + b_z[,n_out_b_pos[i,1:n_bs]] * segment(b_out_pred[g,],k,n_bs), sigma_out[g,i]);
        k = k + n_bs; // update index
        }
    }
  }

}


generated quantities{
  array[G] matrix[n_mvn1,n_mvn1] bcorr;  // random coefficients correlation matrix
  array[G] matrix[n_mvn2,n_mvn2] bcorr2; // random coefficients correlation matrix
  array[G] matrix[D_cen,D_cen] bcorr_inn; // random coefficients correlation matrix

  for(g in 1:G){
    if(n_mvn1 > 0){
      bcorr[g] = multiply_lower_tri_self_transpose(L[g]);
    }
    if(n_mvn2 > 0){
      bcorr2[g] = multiply_lower_tri_self_transpose(L2[g]);
    }
    bcorr_inn[g,] = multiply_lower_tri_self_transpose(L_inno[g,]);
  }
}
