data {
    int N;
    array[N] int x;
}

parameters {
    real<lower=0> a;
    real<lower=0> b;
}

model {
    // Non-informative prior
    target += log(1 / (a * b));
    // Likelihood
    for (i in 1:N){
        target += lgamma(a + x[i]) - lgamma(a) - lgamma(x[i] + 1) + a * log(b) - (a + x[i]) * log(b + 1);
    }
}
