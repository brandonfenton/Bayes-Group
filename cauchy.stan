data{
  // Number of observations, number of levels of each factor, number of
  // Combinations of factors
	int num_obs;
	int m; // Number of Binomial trials
    real<lower=0> cauchysd; // Scale for cauchy priors
	int num_b;
	int num_v;
	int num_n;
	int num_i;
	int num_bv;
	int num_bn;
	int num_bi;
	int num_vn;
	int num_vi;
	int num_ni;
	int num_bvn;
	int num_bvi;
	int num_bni;
	int num_vni;

	// Send in data (actually levels)
	int y[num_obs];
	int blk[num_obs];
	int vty[num_obs];
	int nit[num_obs];
	int ino[num_obs];
	int bv_int[num_obs];
	int bn_int[num_obs];
	int bi_int[num_obs];
	int vn_int[num_obs];
	int vi_int[num_obs];
	int ni_int[num_obs];
	int bvn_int[num_obs];
	int bvi_int[num_obs];
	int bni_int[num_obs];
	int vni_int[num_obs];
}

parameters{
  real mu; // Overall mean
  vector[num_b] b;
  vector[num_v] v;
  vector[num_bv] b_v;
  vector[num_i] i;
  vector[num_n] n;
  vector[num_bi] b_i;
  vector[num_bn] b_n;
  vector[num_ni] n_i;
  vector[num_vi] v_i;
  vector[num_vn] v_n;
  vector[num_bni] b_n_i;
  vector[num_bvi] b_v_i;
  vector[num_bvn] b_v_n;
  vector[num_vni] v_n_i;

  // Standard deviations of the batches
  real <lower=0> sigma_b;    // Batch 1
  real <lower=0> sigma_v;    // Batch 2
  real <lower=0> sigma_bv;   // Batch 3
  real <lower=0> sigma_i;    // Batch 4
  real <lower=0> sigma_n;    // Batch 5
  real <lower=0> sigma_bi;   // Batch 6
  real <lower=0> sigma_bn;   // Batch 7
  real <lower=0> sigma_ni;   // Batch 8
  real <lower=0> sigma_vi;   // Batch 9
  real <lower=0> sigma_vn;   // Batch 10
  real <lower=0> sigma_bni;  // Batch 11
  real <lower=0> sigma_bvi;  // Batch 12
  real <lower=0> sigma_bvn;  // Batch 13
  real <lower=0> sigma_vni;  // Batch 14
}

transformed parameters{
	vector[num_obs] l_p;
	vector[num_obs] p;

	// logit(p) = Xb
	for (o in 1:num_obs) {
  	l_p[o] <- mu + b[blk[o]] + v[vty[o]] + n[nit[o]] + i[ino[o]] +
      b_v[bv_int[o]] + b_n[bn_int[o]] + b_i[bi_int[o]] + v_n[vn_int[o]] +
      v_i[vi_int[o]] + n_i[ni_int[o]] + b_v_n[bvn_int[o]] +
      b_v_i[bvi_int[o]] + b_n_i[bni_int[o]] + v_n_i[vni_int[o]];
		p[o] <- inv_logit(l_p[o]);
	}
}

model{
  // likelihood
	y ~ binomial(m,p);

	// priors on effects
	mu ~ normal(0, 10);  //batch of 1
	b ~ normal(0, sigma_b); //batch of 6
	v ~ normal(0, sigma_v); //batch of 5
	n ~ normal(0, sigma_n); //batch of 3
	i ~ normal(0, sigma_i); //batch of 2
	b_v ~ normal(0, sigma_bv); //batch of 30
	b_n ~ normal(0, sigma_bn); //batch of 18
	b_i ~ normal(0, sigma_bi); //batch of 12
	v_n ~ normal(0, sigma_vn); // batch of 15
	v_i ~ normal(0, sigma_vi); //batch of 10
	n_i ~ normal(0, sigma_ni); //batch of 6
	b_v_n ~ normal(0, sigma_bvn); //batch of 90
	b_v_i ~ normal(0, sigma_bvi); //batch of 60
	b_n_i ~ normal(0, sigma_bni); //batch of 36
	v_n_i ~ normal(0, sigma_vni); // batch of 30

	// priors on standard deviations
 	sigma_b ~ cauchy(0, cauchysd);
	sigma_v ~ cauchy(0, cauchysd);
	sigma_n ~ cauchy(0, cauchysd);
	sigma_i ~ cauchy(0, cauchysd);
	sigma_bv ~ cauchy(0, cauchysd);
	sigma_bn ~ cauchy(0, cauchysd);
	sigma_bi ~ cauchy(0, cauchysd);
	sigma_vn ~ cauchy(0, cauchysd);
  sigma_vi ~ cauchy(0, cauchysd);
	sigma_ni ~ cauchy(0, cauchysd);
	sigma_bvn ~ cauchy(0, cauchysd);
	sigma_bvi ~ cauchy(0, cauchysd);
	sigma_bni ~ cauchy(0, cauchysd);
	sigma_vni ~ cauchy(0, cauchysd);
}
