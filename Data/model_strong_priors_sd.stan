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
    array[K] int<lower=0> time_max; //Duration of obs period for each trial
    array[K, T_max] int<lower=0> D_int; // integer durations
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
generated quantities {
    vector<lower=0>[N_networks] s = s_prime ./ lambda_0;
    matrix[K, Q] log_lik_matrix = rep_matrix(0.0, K, Q);           // LL for each observation
    //for %ST
    int count_ST = 0;
    vector[N_networks] psocn_sum = rep_vector(0.0, N_networks);
    for (trial in 1:K) {
        for (n in 1:N[trial]) {
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];
            if (learn_time > 0){
                real cum_hazard = 0; //set val before adding
                for (time_step in 1:learn_time) {
                    real ind_term = exp(age_group_i[trial,time_step,id] + sex_i[id]);
                    real net_effect = 0;
                    for (network in 1:N_networks) {
                        net_effect += s_prime[network] * dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
                    }
                    real soc_term = net_effect* exp(age_group_s[trial,time_step,id] + sex_s[id]);
                    real base_rate = lambda_0 * ind_term + soc_term;
                    real lambda = fmax(base_rate, 1e-6) * D[trial, time_step];
                    cum_hazard += lambda; // accumulate hazard
                    //if it learn_time, record the ll
                    if (time_step == learn_time){
                        log_lik_matrix[trial, n] = log( (lambda_0 * ind_term + soc_term)) - cum_hazard;
                        for (network in 1:N_networks) {
                            real Tn = dot_product(A[network, trial, time_step][id, ], Z[trial][time_step, ]);
                            psocn_sum[network] += (s_prime[network] * D[trial, time_step] * exp(age_group_s[trial,time_step,id] + sex_s[id])  * Tn) / lambda;
                        }
                        count_ST += 1;
                    }
                }
            }
        }
        // Contributions of censored individuals
        if (N_c[trial] > 0) {
            for (c in 1:N_c[trial]) {
                int id = ind_id[trial, N[trial] + c];
                int censor_time = T[trial]; // Censoring time (end of observation)
                // compute cumulative hazard up to the censoring time
                real cum_hazard = 0;
                for (time_step in 1:censor_time) {
                    real ind_term = exp(age_group_i[trial,time_step,id] + sex_i[id]);
                    real net_effect = 0;
                    for (network in 1:N_networks) {
                        net_effect += s_prime[network] * dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
                    }
                    real soc_term = net_effect* exp(age_group_s[trial,time_step,id] + sex_s[id]);
                    real base_rate = lambda_0 * ind_term + soc_term;
                    real lambda = fmax(base_rate, 1e-6) * D[trial, time_step];
                    cum_hazard += lambda; // accumulate hazard
                }
                // Compute per-individual log likelihood
                log_lik_matrix[trial, N[trial] + c] = -cum_hazard;
            }
        }
    }
    vector[N_networks] percent_ST;

if (count_ST > 0) {
    percent_ST = psocn_sum / count_ST;
} else {
    percent_ST = rep_vector(0.0, N_networks);
}
    matrix[K, Q] acquisition_time;         // simulated acquisition times
    for (trial in 1:K) {
        for (n in 1:Q) {
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];
            // if demonstrator, skip simulation
            if (learn_time < 0) {
                acquisition_time[trial, n] = 0;
                continue;
            }
            real cum_hazard = 0; //set val before adding
            real threshold = -log(uniform_rng(0, 1));
            int global_time = 1;
            acquisition_time[trial, n] = time_max[trial];
            for (time_step in 1:T[trial]) {
                for (micro_time in 1:D_int[trial, time_step]){
                    real ind_term = exp(age_group_i[trial,time_step,id] + sex_i[id]);
                    real net_effect = 0;
                    for (network in 1:N_networks) {
                        net_effect += s_prime[network] * dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
                    }
                    real soc_term = net_effect* exp(age_group_s[trial,time_step,id] + sex_s[id]);
                    real base_rate = lambda_0 * ind_term + soc_term;
                    real lambda = fmax(base_rate, 1e-6);
                    cum_hazard += lambda;
                    if (cum_hazard > threshold) {
                        acquisition_time[trial, n] = global_time;
                        break;  // exit inner loop
                    }
                    global_time += 1;
                }
                if (cum_hazard > threshold) break;  // exit outer loop
            }
        }
    }
    // Flatten log_lik_matrix into log_lik
    array[K * Q] real log_lik;
    int idx = 1;
    for (trial in 1:K) {
        for (n in 1:Q) {
            log_lik[idx] = log_lik_matrix[trial, n];
            idx += 1;
        }
    }
}
