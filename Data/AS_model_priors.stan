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

    array[K,T_max] real<lower=0> D;

    array[K,T_max] matrix[P,3] ILV_age_group;
    vector[P] ILV_sex;

    array[K] int<lower=0> time_max;
    array[K,T_max] int<lower=0> D_int;
}

parameters {

    real log_lambda_0_mean;

    vector[3] beta_ILVi_age_group;

    real beta_ILVi_sex;
}

transformed parameters {

    real<lower=0> lambda_0 =
        exp(log_lambda_0_mean);

    array[K,T_max] vector[P] age_group_i;

    vector[P] sex_i =
        ILV_sex * beta_ILVi_sex;

    for (trial in 1:K) {

        for (timestep in 1:T_max) {

            age_group_i[trial][timestep] =
                ILV_age_group[trial][timestep] *
                beta_ILVi_age_group;
        }
    }
}

model {

    log_lambda_0_mean ~ normal(-4,0.5);

    beta_ILVi_age_group ~ normal(0,1);

    beta_ILVi_sex ~ normal(0,1);

    for (trial in 1:K) {

        // learned individuals

        for (n in 1:N[trial]) {

            int id =
                ind_id[trial,n];

            int learn_time =
                t[trial,id];

            if (learn_time > 0) {

                for (time_step in 1:learn_time) {

                    real ind_term =
                        exp(
                            age_group_i[trial,time_step,id]
                            + sex_i[id]
                        );

                    real base_rate =
                        lambda_0 * ind_term;

                    real lambda =
                        fmax(base_rate,1e-6) *
                        D[trial,time_step];

                    target += -lambda;

                    if (time_step == learn_time)
                        target += log(base_rate);
                }
            }
        }

        // censored individuals

        for (c in 1:N_c[trial]) {

            int id =
                ind_id[trial,N[trial] + c];

            for (time_step in 1:T[trial]) {

                real ind_term =
                    exp(
                        age_group_i[trial,time_step,id]
                        + sex_i[id]
                    );

                real lambda =
                    fmax(
                        lambda_0 * ind_term,
                        1e-6
                    ) *
                    D[trial,time_step];

                target += -lambda;
            }
        }
    }
}

generated quantities {

    matrix[K,Q] log_lik_matrix =
        rep_matrix(0.0,K,Q);

    for (trial in 1:K) {

        for (n in 1:N[trial]) {

            int id =
                ind_id[trial,n];

            int learn_time =
                t[trial,id];

            if (learn_time > 0) {

                real cum_hazard = 0;

                for (time_step in 1:learn_time) {

                    real ind_term =
                        exp(
                            age_group_i[trial,time_step,id]
                            + sex_i[id]
                        );

                    real rate =
                        lambda_0 * ind_term;

                    real lambda =
                        fmax(rate,1e-6) *
                        D[trial,time_step];

                    cum_hazard += lambda;

                    if (time_step == learn_time) {

                        log_lik_matrix[trial,n] =
                            log(rate) -
                            cum_hazard;
                    }
                }
            }
        }

        for (c in 1:N_c[trial]) {

            int id =
                ind_id[trial,N[trial] + c];

            real cum_hazard = 0;

            for (time_step in 1:T[trial]) {

                real ind_term =
                    exp(
                        age_group_i[trial,time_step,id]
                        + sex_i[id]
                    );

                real lambda =
                    fmax(
                        lambda_0 * ind_term,
                        1e-6
                    ) *
                    D[trial,time_step];

                cum_hazard += lambda;
            }

            log_lik_matrix[
                trial,
                N[trial] + c
            ] = -cum_hazard;
        }
    }

    array[K * Q] real log_lik;

    int idx = 1;

    for (trial in 1:K) {

        for (n in 1:Q) {

            log_lik[idx] =
                log_lik_matrix[trial,n];

            idx += 1;
        }
    }
}