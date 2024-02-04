functions {
    real custom_function_lpdf(real x, real a, real b) {
        return lgamma(a + x) - lgamma(a) - lgamma(x + 1) + a * log(b) - (a + x) * log(b + 1);
    }
}

data {
    int<lower=1> N;     // Number of samples
    array[N] int x;     // Data
}

parameters {
    real<lower=0> a;    // Shape parameter
    real<lower=0> b;    // Scale parameter
}

model {
    // Non-informative priors
    a ~ uniform(0, 1000);
    b ~ uniform(0, 1000);
    
    // Likelihood
    for (i in 1:N) {
        target += custom_function_lpdf(x[i] | a, b);
    }
}
