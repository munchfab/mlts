// y_merge --> array der Länge D mit Vektoren der Länge n_obs: beobachtete Werte
// D --> Anzahl der time varying constructs
// n_miss_D --> array der Länge D mit Anzahl der missings pro D
// pos_miss_D --> 2 dimensionales array unterschiedlicher Spaltenzahl mit den
  // Positionen der missings
// y_impute --> PARAMETER mit dem missing Werten

array[] vector missings_and_censoring (array[] vector y,
                                        array[] int n_miss_D,
                                        array[,] int pos_miss_D,
                                        vector y_impute) {
  int p_miss = 1; // running counter variable zum indexieren von y_impute
  int N_obs = num_elements(y[1]);
  int D = size(y);

  array[D] vector[N_obs] y_out = y;

  for(i in 1:D){
    if(n_miss_D[i]>0){
    // add imputed values for missings on each indicator
    y_out[i,pos_miss_D[i,1:n_miss_D[i]]] = segment(y_impute, p_miss, n_miss_D[i]);
    p_miss = p_miss + n_miss_D[i];    // update counter for next indicator i+1
    }
  }
  return y_out;
}

