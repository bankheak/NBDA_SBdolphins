//stan
data {
    int<lower=0> K;                // Number of trials
    int<lower=0> Q;                // Number of individuals in each trial
    int<lower=1> P;                // Number of unique individuals
    array[K] int<lower=0> N;       // Number of individuals that learned during observation period
    array[K] int<lower=0> N_c;     // Number of right-censored individuals
    array[K, Q] int<lower=-1> ind_id; // IDs of individuals
    array[K] int<lower=1> T;       // Maximum time periods
    int<lower=1> T_max;            // Max timesteps reached
    array[K,P] int t;     // Time of acquisition for each individual
    array[K, T_max] real<lower=0> D; // Scaled durations
    int<lower=1> N_networks;
    array[N_networks, K, T_max] matrix[P, P] A;  // network matrices
    array[K] matrix[T_max, P] Z;   // Knowledge state * cue matrix
    array[K] matrix[T_max, P] Zn;   // Knowledge state
    array[K,T_max] matrix[P, 3] ILV_age_group;
    vector[P] ILV_sex;
    int<lower=0> N_veff;
}
parameters {
    real log_lambda_0_mean;  // Log baseline learning rate
    vector[N_networks] log_s_prime_mean;
    vector[3] beta_ILVi_age_group;
    real beta_ILVi_sex;
    vector[3] beta_ILVs_age_group;
    real beta_ILVs_sex;
}
transformed parameters {
    vector<lower=0>[N_networks] s_prime;
    real<lower=0> lambda_0;
    array[K, T_max] vector[P] age_group_i;
    vector[P] sex_i;
    array[K, T_max] vector[P] age_group_s;
    vector[P] sex_s;
    s_prime = exp(log_s_prime_mean);
    lambda_0 = exp(log_lambda_0_mean);
    sex_i = ILV_sex * beta_ILVi_sex;
    sex_s = ILV_sex * beta_ILVs_sex;
    for (trial in 1:K) {
        for (timestep in 1:T_max) age_group_i[trial][timestep] = ILV_age_group[trial][timestep] * beta_ILVi_age_group;
        for (timestep in 1:T_max) age_group_s[trial][timestep] = ILV_age_group[trial][timestep] * beta_ILVs_age_group;
    }
}
model {
    log_lambda_0_mean ~ normal(-4, 0.5);
    log_s_prime_mean ~ normal(-4, 0.5);
    beta_ILVi_age_group ~ normal(0,1);
    beta_ILVi_sex ~ normal(0,1);
    beta_ILVs_age_group ~ normal(0,1);
    beta_ILVs_sex ~ normal(0,1);
    for (trial in 1:K) {
        for (n in 1:N[trial]) {
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];
            if (learn_time > 0) {
                for (time_step in 1:learn_time) {
                    real ind_term = exp(age_group_i[trial,time_step,id] + sex_i[id]);
                    real net_effect = 0;
                    for (network in 1:N_networks) {
                        net_effect += s_prime[network] * dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
                    }
                    real soc_term = net_effect* exp(age_group_s[trial,time_step,id] + sex_s[id]);
                    real base_rate = lambda_0 * ind_term + soc_term;
                    real lambda = fmax(base_rate, 1e-6) * D[trial, time_step];
                    target += -lambda;
                    if (time_step == learn_time) {
                        target += log( (lambda_0 * ind_term + soc_term));
                    }
                }
            }
        }
        if (N_c[trial] > 0) {
            for (c in 1:N_c[trial]) {
                int id = ind_id[trial, N[trial] + c];
                for (time_step in 1:T[trial]) {
                    real ind_term = exp(age_group_i[trial,time_step,id] + sex_i[id]);
                    real net_effect = 0;
                    for (network in 1:N_networks) {
                        net_effect += s_prime[network] * dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
                    }
                    real soc_term = net_effect* exp(age_group_s[trial,time_step,id] + sex_s[id]);
                    real base_rate = lambda_0 * ind_term + soc_term;
                    real lambda = fmax(base_rate, 1e-6) * D[trial, time_step];
                    target += -lambda;
                }
            }
        }
    }
}
