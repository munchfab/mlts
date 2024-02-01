// autoregressive DSEM with manifest variables
data {
  int<lower=1> N; 	// number of observational units
  int<lower=2> D; 	// number of time-varying constructs
  int<lower=1> N_obs; 	// observations in total: N * TP
  int<lower=1> n_pars;
  int<lower=D> n_random;   // number of random effects
  int n_fixed;
  int is_fixed[1,n_fixed];
  int is_random[n_random];  // which parameters to model person-specific
  int<lower=1> N_obs_id[N]; // number of observations for each unit
  vector[N_obs] y[D]; 	    // D*N_obs array of observations

  // handling of missing values
  int n_miss;                      // total number of missings across D
  int n_miss_D[D];                 // missings per D
  int pos_miss_D[D,max(n_miss_D)]; // array of missings' positions

  // model adaptations based on user inputs:
  // - fixing parameters to constant values:
  // - innovation variances
  int<lower=0,upper=1> innos_rand[D];
  int n_innos_fix;
  int innos_fix_pos[D];
  int innos_pos[D];
  
  // - time-invariant variables:
  int n_cov;           // number of covariates - minimum of 1 for intercepts
  matrix[N, n_cov] W;  // predictors of individual parameters

  // - dynamic model specification per D
  int<lower=1> N_pred[D];     // Number of predictors per dimension
  int<lower=0> D_pred[D,D];   // matrix to index predictors to use per dimension
  int Dpos1[D];  // index positions of danymic effect parameters
  int Dpos2[D];

  // - covariances of innovations:
  int n_inno_covs; // number of potential innovation covs to include
//  int<lower=0,upper=1> inno_cov0;    // fixed to zero
//  int<lower=0,upper=1> inno_cov_fix; // fixed to zero
  int n_obs_cov;   // total number of residuals 
  int inno_cov_pos[1,n_inno_covs];
}

parameters {
  vector[n_random] b_free[N];              // person-specific parameter
  vector<lower=0>[n_random] sd_R;          // random effect SD   
  vector<lower=0>[D-sum(innos_rand)] sigma;// SDs of fixed innovation variances 
  cholesky_factor_corr[n_random] L;        // cholesky factor of random effects correlation matrix
  vector[n_miss] y_impute;                 // vector to store imputed values
  row_vector[n_random] gammas;             // fixed effect (intercepts)
  matrix[n_random,n_cov - 1] btw_free;
  vector[n_fixed] b_fix;
  
  vector[n_obs_cov] eta_cov[n_inno_covs];
}

transformed parameters {
  matrix[N, n_random] bmu;     // gammas of person-specific parameters
  vector[n_pars] b[N];
  // array of innovation variances - in case of random innovation variances
  // retransform by exponentiation of log innovations sampled from MVN
  vector<lower = 0>[N] sd_noise[D];
  vector<lower = 0>[N] sd_inncov[n_inno_covs];

  // regression weights of time-invariant variables as predictors of random effects
  matrix[n_cov, n_random] btw_pred;
  btw_pred[1,] = gammas;
  if(n_cov>1){
   btw_pred[2,] = to_row_vector(btw_free[,1]);
 }
  // calculate population means (intercepts) of person-specific parameters
  bmu = W * btw_pred;
  
  // create array of (person-specific) parameters to use in model
  for(i in 1:n_random){
    b[,is_random[i]] = b_free[,i];
  }
  if(n_fixed>0){
    for(i in 1:n_fixed){
      b[,is_fixed[1,i]] = to_array_1d(rep_vector(b_fix[i],N));
    }
  }
  
  // transformation of log-innovation variances if modeled as person-specific
  for(i in 1:D){
    if(innos_rand[i] == 0){
      sd_noise[i,] = rep_vector(sigma[innos_fix_pos[i]],N);
 //     b[,innos_pos[i]] = to_array_1d(sd_noise[i,]); // probably redundant
    } else {
      sd_noise[i,] = sqrt(exp(to_vector(b[,innos_pos[i]])));
    }
  }
  
  // transform log innovation covarainces 
  for(i in 1:n_inno_covs){
    sd_inncov[i,1:N] = sqrt(exp(to_vector(b[,inno_cov_pos[1,i]])));
    }
 }

model {
  int pos = 1;       // initialize position indicator
  int pos_cov = 1;   // covariance position
  int p_miss = 1;    // running counter variable to index positions on y_impute
  int obs_id = 1;    // declare local variable to store variable number of obs per person
  matrix[n_random, n_random] SIGMA = diag_pre_multiply(sd_R, L);
  vector[N_obs] y_merge[D];

  y_merge = y;      // add observations
  if(n_miss>0){
    for(i in 1:D){
    // add imputed values for missings on each indicator
    y_merge[i,pos_miss_D[i,1:n_miss_D[i]]] = segment(y_impute, p_miss, n_miss_D[i]);
    p_miss = p_miss + n_miss_D[i];    // update counter for next indicator i+1
  }
  }

  // priors on elements in the random effect covariance matrix
  L ~ lkj_corr_cholesky(1);  // maybe include option to specify priors in fit function?

  // priors on (random effects) variances
  sigma ~ cauchy(0, 2.5);

  // priors on fixed effects
  for(i in 1:n_random){
    gammas[i] ~ normal(0,10);
    if(n_cov > 1){
    btw_free[i,] ~ normal(0,2);
    }
  }

  if(n_fixed > 0){
    b_fix ~ normal(0,2);
  }


  for (pp in 1:N) {
    // store number of observations per person
    obs_id = (N_obs_id[pp]);

    // individual parameters from multivariate normal distribution
    b_free[pp, 1:n_random] ~ multi_normal_cholesky(bmu[pp, 1 : n_random], SIGMA);

    // local variable declaration: array of predicted values
    vector[obs_id-1] mus[D];
    
    // for now:
    // create latent mean centered versions of observations
    vector[obs_id] y_cen[D];

    for(d in 1:D){ // start loop over dimensions
    y_cen[d,] = y_merge[d,pos:(pos+obs_id-1)] - b[pp,d];
    }
    
    // inno covs 
    vector[obs_id-1] eta_cov_id[n_inno_covs];
    for(i in 1:n_inno_covs){
      eta_cov_id[i,] = segment(eta_cov[i,], pos_cov, (obs_id-1));
      eta_cov_id[i,] ~ normal(0, sd_inncov[i,pp]);
    } 
    
    
    for(d in 1:D){ // start loop over dimensions
      // build prediction matrix for specific dimensions
      {
      matrix[(obs_id-1),N_pred[d]] b_mat; // adjust for non-fully crossed models
      for(nd in 1:N_pred[d]){ // start loop over number of predictors in each dimension
        b_mat[,nd] = y_cen[D_pred[d, nd],1:(obs_id-1)];
      }

      // use build predictor matrix to calculate latent time-series means
      mus[d,] =  b[pp,d] + b_mat * b[pp, Dpos1[d]:Dpos2[d]] + eta_cov_id[1,];
      }
      // sampling statement
      y_merge[d,(pos+1):(pos+(obs_id-1))] ~ normal(mus[d,], sd_noise[d,pp]);
    }
    // update index variables
    pos = pos + obs_id;
    pos_cov = pos_cov + obs_id - 1;
  }

}

generated quantities{
  matrix[n_random,n_random] bcorr; // random coefficients correlation matrix
//  matrix[n_random,n_random] bcov; // random coefficients covariance matrix
  // create random coefficient matrices
  bcorr = multiply_lower_tri_self_transpose(L);
//  bcov = quad_form_diag(bcorr, sd_R);
}
