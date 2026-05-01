// autoregressive DSEM with manifest variables
data {
  int<lower=1> N; 	        // number of observational units
  int<lower=1> G;   // number of groups
  int<lower=1> D;           // number of latent constructs
  int<lower=1> D_cen;
  array[D] int<lower=1> D_np;     // number of indicators per construct
  int<lower=1> n_p; 	      // number of manifest indicators
  array[n_p] int D_perP;          // indicate dimension per indicator
  array[n_p] int is_SI;           // indicate if single-indicator per construct
  array[D] int D_pos_is_SI;       // indicate position of single-indicator per construct
  int<lower=0, upper=3> maxLag; // maximum lag
  int<lower=1> N_obs; 	    // observations in total: N * TP
  int<lower=1> n_pars;
  int<lower=D_cen> n_random;    // number of random effects

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

  // - covariances of innovations:
  int n_inno_covs; // number of potential innovation covs to include
  int n_obs_cov;   // total number of residuals
  array[1,n_inno_covs] int inno_cov_pos;
  array[2] int<lower=-1,upper=1> inno_cov_load;

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

  array[D] int<lower=0,upper=1> is_rdsem;
  array[D] int<lower=0> N_pred_rdsem;
  array[D,max(N_pred_rdsem)] int<lower=0> D_pred_rdsem;
  array[D] int Dpos1_rdsem;
}

transformed data{
  int n_SD_etaW_all;
  int n_SD_etaW_i;
  n_SD_etaW_all = standardized == 1 ? n_etaW_free : 0;
  n_SD_etaW_i = standardized == 1 ? N : 0;
}

parameters {
  array[N] vector[n_random] b_free;            // person-specific parameter
  array[G] vector<lower=0>[n_random] sd_R;        // random effect SD
  array[G] vector<lower=0>[n_innos_fix] sigma;    // SDs of fixed innovation variances
  array[G] cholesky_factor_corr[n_mvn1] L;      // cholesky factor of random effects correlation matrix
  array[G] cholesky_factor_corr[n_mvn2] L2;     // cholesky factor of random effects correlation matrix
  vector[n_miss] y_impute;               // vector to store imputed values
  vector<upper=censL_val>[n_censL] y_impute_censL;
  vector<lower=censR_val>[n_censR] y_impute_censR;
  array[G] row_vector[n_random] gammas;           // fixed effect (intercepts)
  array[G] vector[n_cov_bs] b_re_pred;            // regression coefs of RE prediction
  array[G] vector[n_fixed] b_fix;
  array[G] vector[n_out] alpha_out;               // outcome precition intercepts
  array[G] vector<lower=0>[n_out] sigma_out;      // residual SD(s) of outcome(s)
  array[G] vector[n_out_bs_sum] b_out_pred;       // regression coefs of out prediction
  array[n_inno_covs] vector[n_obs_cov] eta_cov;

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
  matrix[N,n_pars-n_innos_fix] b;
  array[D_cen] vector[N] sd_noise;
  array[n_inno_covs] vector[N] sd_inncov;
  array[G] matrix[n_cov, n_random] b_re_pred_mat;
  vector[n_p] loadB = rep_vector(1, n_p); // measurement model parameters
  vector[n_p] loadW = rep_vector(1, n_p);
  vector[n_p] alpha = rep_vector(0, n_p);
  vector[n_p] sigmaB = rep_vector(0, n_p);
  vector[n_p] sigmaW = rep_vector(0, n_p);

 // REs regressed on covariates
  b_re_pred_mat[1] = rep_matrix(0, n_cov, n_random);
  b_re_pred_mat[1,1,] = gammas[1,];
  if(n_cov>1){
     for(i in 1:n_cov_bs){
     b_re_pred_mat[1,n_cov_mat[i,1],n_cov_mat[i,2]] = b_re_pred[1,i];
    }
  }
  // calculate population means (intercepts) of person-specific parameters
  bmu = W * b_re_pred_mat[1,,];

  // create array of (person-specific) parameters to use in model
  for(i in 1:n_random){
    b[,is_random[i]] = to_vector(b_free[,i]);
  }
  if(n_fixed>0){
    for(i in 1:n_fixed){
      b[,is_fixed[1,i]] = rep_vector(b_fix[1,i],N);
    }
  }

  // transformation of log-innovation variances if modeled as person-specific
  for(i in 1:D_cen){
    if(innos_rand[i] == 0){
      sd_noise[i,] = rep_vector(sigma[1,innos_fix_pos[i]],N);
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
  int pos_cov = 1;   // covariance position
  int p_miss = 1;    // running counter variable to index positions on y_impute
  int p_censL = 1;
  int p_censR = 1;
  int obs_id = 1;    // declare local variable to store variable number of obs per person
  array[n_p] vector[N_obs] y_merge;
  array[max(p_is_wcen_pos)] vector[N_obs] Ymus;
  array[n_p] vector[N] YB;
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

  // add imputed values for missings on each indicator
  y_merge = y;
  for(i in 1:n_p){
    if(n_miss_p[i]>0){
      y_merge[i,pos_miss_p[i,1:n_miss_p[i]]] = segment(y_impute, p_miss, n_miss_p[i]);
      p_miss = p_miss + n_miss_p[i];
    }
  }

  // replace values at censor thresholds
  for(i in 1:n_p){
    if(n_censL_p[i]>0){
    // add imputed values for observations at floor (threshold for censoring)
    y_merge[i,pos_censL_p[i,1:n_censL_p[i]]] = segment(y_impute_censL, p_censL, n_censL_p[i]);
    p_censL = p_censL + n_censL_p[i];
    }
    if(n_censR_p[i]>0){
    // add imputed values for observations at ceiling (threshold for censoring)
    y_merge[i,pos_censR_p[i,1:n_censR_p[i]]] = segment(y_impute_censR, p_censR, n_censR_p[i]);
    p_censR = p_censR + n_censR_p[i];
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

    int pp_g = 1;

    // store number of observations per person
    obs_id = (N_obs_id[pp]);
    int pos_etaW_free = 1;    // running counter variable to index positition on etaW_free
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

    if( sum(is_rdsem) > 0 ){
      for (d in 1:D) {
        if (is_rdsem[d] == 1 ) {
          for(k in 1:N_pred_rdsem[d]){
             etaW_id[d,] = etaW_id[d,] - b[pp,(Dpos1_rdsem[d]+(k-1))] * etaW_id[D_pred_rdsem[d,k],];
          }
        }
      }
    }


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

    // dynamic process
    array[D_cen] vector[obs_id-maxLag] mus;

    array[n_inno_covs] vector[obs_id-maxLag] eta_cov_id;
    if(n_inno_covs > 0){
       for(i in 1:n_inno_covs){
          eta_cov_id[i,] = segment(eta_cov[i,], pos_cov, (obs_id-maxLag));
          target += normal_lpdf(eta_cov_id[i,] | 0, sd_inncov[i,pp]);
          }
      }

    for(d in 1:D){ // start loop over dimensions
      if(is_wcen[d] == 1){

      if(N_pred[d]>0){
        // build prediction matrix for specific dimensions
        int n_cols = (n_inno_covs>0 && d<3) ? N_pred[d]+n_inno_covs : N_pred[d];
        matrix[(obs_id-maxLag),n_cols] b_mat;
        vector[n_cols] b_use;
        b_use[1:N_pred[d]] = to_vector(b[pp, Dpos1[d]:Dpos2[d]]);

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

        if(n_inno_covs>0&&d<3){
          b_use[N_pred[d]+1] = inno_cov_load[d];
          b_mat[,(N_pred[d]+1)] = eta_cov_id[1,]; // add innovation covariance factor scores
          }

          // use build predictor matrix to calculate latent time-series means
          if(D_np[d] == 1){
              //mus[D_cen_pos[d],] = YB[D_pos_is_SI[d],pp] + b_mat * b_use;
              mus[D_cen_pos[d],] = b_mat * b_use;
              } else {
              mus[D_cen_pos[d],] = b_mat * b_use;
              }

        } else {
          mus[D_cen_pos[d],] = rep_vector(0,(obs_id-maxLag));
          }

        target += normal_lpdf(etaW_id[d,(1+maxLag):obs_id] | mus[D_cen_pos[d],], sd_noise[D_cen_pos[d],pp]);
      }
    } // end loop over dimensions

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
    // update index variables
    pos = pos + obs_id;
    pos_cov = pos_cov + obs_id - maxLag;
    } // end loop over subjects

    // sampling statements
    for(i in 1:n_p){
      if(p_is_wcen[i] == 1){
        if(mu_is_etaB[i] == 0){
          target += normal_lpdf(YB[i,] | alpha[i] + loadB[i]*b[,mu_etaB_pos[i]], sigmaB[i]);
        }
        if(is_SI[i] == 0){
          target += normal_lpdf(y_merge[i, pos_obs_p[i,1:n_obs_p[i]]] | Ymus[p_is_wcen_pos[i],pos_obs_p[i,1:n_obs_p[i]]], sigmaW[i]);
        }
      }
    }

  // outcome prediction: get expectations of outcome values
  if(n_out > 0){
    int k = 1;
    matrix[N,n_random+n_z] b_z = append_col(b[,is_random],Z);
    for(i in 1:n_out){
      int n_bs = n_out_bs[i,1];      // number of predictors for each outcome
      target += normal_lpdf(out[i,] | alpha_out[1,i] + b_z[,n_out_b_pos[i,1:n_bs]] * segment(b_out_pred[1],k,n_bs), sigma_out[1,i]);
      k = k + n_bs; // update index
    }
  }
}

generated quantities{
  array[G] matrix[n_mvn1,n_mvn1] bcorr;  // random coefficients correlation matrix
  array[G] matrix[n_mvn2,n_mvn2] bcorr2; // random coefficients correlation matrix
  vector[n_SD_etaW_all] SD_etaW;
  array[n_SD_etaW_all] vector[n_SD_etaW_i] SD_etaW_i;
  for(g in 1:G){
    if(n_mvn1 > 0){
      bcorr[g] = multiply_lower_tri_self_transpose(L[g]);
    }
    if(n_mvn2 > 0){
      bcorr2[g] = multiply_lower_tri_self_transpose(L2[g]);
    }
  }

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
