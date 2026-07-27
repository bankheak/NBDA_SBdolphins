data {
    int<lower=0> K;
    int<lower=0> Q;
    int<lower=1> P;
    array[K] int<lower=0> N;
    array[K] int<lower=0> N_c;
    array[K, Q] int<lower=-1> ind_id;
    array[K] int<lower=1> T;
    int<lower=1> T_max;
    array[K,P] int t;
    array[K, T_max] real<lower=0> D;
    int<lower=1> N_networks;

    // Precomputed network effects
    array[N_networks, K, T_max] vector[P] E;

    array[K] matrix[T_max, P] Z;
    array[K] matrix[T_max, P] Zn;

    array[K,T_max] matrix[P, 3] ILV_age_group;
    vector[P] ILV_sex;

    int<lower=0> N_veff;
    array[K] int<lower=0> time_max;
    array[K, T_max] int<lower=0> D_int;
}

parameters {
    real log_lambda_0_mean;
    vector[N_networks] log_s_prime_mean;

    vector[3] beta_ILVi_age_group;
    real beta_ILVi_sex;

    vector[3] beta_ILVs_age_group;
    real beta_ILVs_sex;
}

transformed parameters {
    vector<lower=0>[N_networks] s_prime = exp(log_s_prime_mean);
    real<lower=0> lambda_0 = exp(log_lambda_0_mean);

    array[K, T_max] vector[P] age_group_i;
    array[K, T_max] vector[P] age_group_s;

    vector[P] sex_i = ILV_sex * beta_ILVi_sex;
    vector[P] sex_s = ILV_sex * beta_ILVs_sex;

    for (trial in 1:K) {
        for (timestep in 1:T_max) {

            age_group_i[trial][timestep] =
                ILV_age_group[trial][timestep] *
                beta_ILVi_age_group;

            age_group_s[trial][timestep] =
                ILV_age_group[trial][timestep] *
                beta_ILVs_age_group;
        }
    }
}

model {

    // Priors
    log_lambda_0_mean ~ normal(-4, 0.5);
    log_s_prime_mean ~ normal(-4, 0.5);

    beta_ILVi_age_group ~ normal(0, 1);
    beta_ILVi_sex ~ normal(0, 1);

    beta_ILVs_age_group ~ normal(0, 1);
    beta_ILVs_sex ~ normal(0, 1);

    for (trial in 1:K) {

        // Learned individuals
        for (n in 1:N[trial]) {

            int id = ind_id[trial, n];
            int learn_time = t[trial, id];

            if (learn_time > 0) {

                for (time_step in 1:learn_time) {

                    real ind_term =
                        exp(age_group_i[trial,time_step,id]
                            + sex_i[id]);

                    real net_effect = 0;

                    for (network in 1:N_networks) {
                        net_effect +=
                            s_prime[network] *
                            E[network, trial, time_step][id];
                    }

                    real soc_term =
                        net_effect *
                        exp(age_group_s[trial,time_step,id]
                            + sex_s[id]);

                    real base_rate =
                        lambda_0 * ind_term +
                        soc_term;

                    real lambda =
                        fmax(base_rate, 1e-6) *
                        D[trial, time_step];

                    target += -lambda;

                    if (time_step == learn_time)
                        target += log(base_rate);
                }
            }
        }

        // Censored individuals
        for (c in 1:N_c[trial]) {

            int id = ind_id[trial, N[trial] + c];

            for (time_step in 1:T[trial]) {

                real ind_term =
                    exp(age_group_i[trial,time_step,id]
                        + sex_i[id]);

                real net_effect = 0;

                for (network in 1:N_networks) {
                    net_effect +=
                        s_prime[network] *
                        E[network, trial, time_step][id];
                }

                real soc_term =
                    net_effect *
                    exp(age_group_s[trial,time_step,id]
                        + sex_s[id]);

                real base_rate =
                    lambda_0 * ind_term +
                    soc_term;

                real lambda =
                    fmax(base_rate, 1e-6) *
                    D[trial, time_step];

                target += -lambda;
            }
        }
    }
}

generated quantities {

    vector<lower=0>[N_networks] s =
        s_prime ./ lambda_0;

    matrix[K, Q] log_lik_matrix =
        rep_matrix(0.0, K, Q);

    int count_ST = 0;
    vector[N_networks] psocn_sum =
        rep_vector(0.0, N_networks);

    for (trial in 1:K) {

        // Learned individuals
        for (n in 1:N[trial]) {

            int id = ind_id[trial, n];
            int learn_time = t[trial, id];

            if (learn_time > 0) {

                real cum_hazard = 0;

                for (time_step in 1:learn_time) {

                    real ind_term =
                        exp(age_group_i[trial,time_step,id]
                            + sex_i[id]);

                    real net_effect = 0;

                    for (network in 1:N_networks) {
                        net_effect +=
                            s_prime[network] *
                            E[network, trial, time_step][id];
                    }

                    real soc_term =
                        net_effect *
                        exp(age_group_s[trial,time_step,id]
                            + sex_s[id]);

                    real rate =
                        lambda_0 * ind_term +
                        soc_term;

                    real lambda =
                        fmax(rate, 1e-6) *
                        D[trial, time_step];

                    cum_hazard += lambda;

                    if (time_step == learn_time) {

                        log_lik_matrix[trial, n] =
                            log(rate) - cum_hazard;

                        for (network in 1:N_networks) {

                            real Tn =
                                E[network, trial, time_step][id];

                            real contribution =
                                (s_prime[network]
                                 * D[trial, time_step]
                                 * exp(age_group_s[trial,time_step,id]
                                       + sex_s[id])
                                 * Tn) /
                                fmax(lambda, 1e-6);

                            psocn_sum[network] +=
                                contribution;
                        }

                        count_ST += 1;
                    }
                }
            }
        }

        // Censored individuals
        for (c in 1:N_c[trial]) {

            int id = ind_id[trial, N[trial] + c];
            real cum_hazard = 0;

            for (time_step in 1:T[trial]) {

                real ind_term =
                    exp(age_group_i[trial,time_step,id]
                        + sex_i[id]);

                real net_effect = 0;

                for (network in 1:N_networks) {
                    net_effect +=
                        s_prime[network] *
                        E[network, trial, time_step][id];
                }

                real soc_term =
                    net_effect *
                    exp(age_group_s[trial,time_step,id]
                        + sex_s[id]);

                real lambda =
                    fmax(lambda_0 * ind_term + soc_term,
                         1e-6) *
                    D[trial, time_step];

                cum_hazard += lambda;
            }

            log_lik_matrix[trial,
                           N[trial] + c] =
                -cum_hazard;
        }
    }

    vector[N_networks] percent_ST;

    if (count_ST > 0)
        percent_ST =
            psocn_sum / count_ST;
    else
        percent_ST =
            rep_vector(0.0, N_networks);

    array[K * Q] real log_lik;

    int idx = 1;

    for (trial in 1:K) {
        for (n in 1:Q) {

            log_lik[idx] =
                log_lik_matrix[trial, n];

            idx += 1;
        }
    }
}