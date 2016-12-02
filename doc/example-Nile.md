
## Nile

This is a short ($n = 100$) univariate time series of the annual flow volumes of
the Nile River at Aswan between 1871 and 1970.
This series is described in @DurbinKoopman2012 and had been analyzed by @Cobb1978 and @Balke1993, and in numerous time series textbooks.
A notable feature of the series is a seeming structural break in 1899, around the time of the completion of the Aswan dam.


```r
data("Nile", package = "datasets")
Nile_ <- data_frame(year = year(as.Date(Nile)),
                    flow = as.numeric(Nile),
                    obs = seq_along(Nile))
```


```r
ggplot(Nile_, aes(x = year, y = flow)) +
  geom_point() +
  geom_line() +
  ylab("Annual Flow") + xlab("")
```

<img src="example-Nile_files/figure-html/Nile_plot-1.png" width="672" />

### Local Level Model

The Nile data can be modeled as a local level model,
$$
\begin{aligned}[t]
y_t &= \mu_t + \varepsilon_t & \varepsilon_t & \sim N(0, \sigma_{\varepsilon}^2) \\
\mu_{t + 1} &= \mu_t + \eta_t &
\eta_t & \sim N(0, \sigma^2_{\eta})
\end{aligned}
$$


```
functions {
  #include ssm.stan
}
data {
  int<lower = 1> n;
  vector[1] y[n];
  vector<lower = 0.0>[1] a1;
  cov_matrix[1] P1;
  real<lower = 0.0> y_scale;
}
transformed data {
  // system matrices
  matrix[1, 1] T;
  matrix[1, 1] Z;
  matrix[1, 1] R;
  vector[1] c;
  vector[1] d;
  int m;
  int p;
  int q;
  int filter_sz;
  m = 1;
  p = 1;
  q = 1;
  T[1, 1] = 1.0;
  Z[1, 1] = 1.0;
  R[1, 1] = 1.0;
  c[1] = 0.0;
  d[1] = 0.0;
  filter_sz = ssm_filter_size(m, p);
}
parameters {
  real<lower = 0.0> sigma_eta;
  real<lower = 0.0> sigma_epsilon;
}
transformed parameters {
  matrix[1, 1] H;
  matrix[1, 1] Q;
  H[1, 1] = pow(sigma_epsilon, 2);
  Q[1, 1] = pow(sigma_eta * sigma_epsilon, 2);
}
model {
  y ~ ssm_constant_lpdf(d, Z, H, c, T, R, Q, a1, P1);
  sigma_epsilon ~ cauchy(0.0, y_scale);
  sigma_eta ~ cauchy(0.0, 1.0);
}
generated quantities {
  vector[filter_sz] filtered[n];
  vector[1] alpha[n];
  vector[1] eta[n];
  vector[1] eps[n];
  // filtering
  filtered = ssm_filter(y,
    rep_array(d, 1),
    rep_array(Z, 1),
    rep_array(H, 1),
    rep_array(c, 1),
    rep_array(T, 1),
    rep_array(R, 1),
    rep_array(Q, 1), a1, P1);
  // sampling states
  alpha = ssm_simsmo_states_rng(filtered,
    rep_array(d, 1),
    rep_array(Z, 1),
    rep_array(H, 1),
    rep_array(c, 1),
    rep_array(T, 1),
    rep_array(R, 1),
    rep_array(Q, 1), a1, P1);
  eps = ssm_simsmo_eps_rng(filtered,
    rep_array(d, 1),
    rep_array(Z, 1),
    rep_array(H, 1),
    rep_array(c, 1),
    rep_array(T, 1),
    rep_array(R, 1),
    rep_array(Q, 1), a1, P1);
  eta = ssm_simsmo_eta_rng(filtered,
    rep_array(d, 1),
    rep_array(Z, 1),
    rep_array(H, 1),
    rep_array(c, 1),
    rep_array(T, 1),
    rep_array(R, 1),
    rep_array(Q, 1), a1, P1);
}
```


```r
local_level_mod <- ssm_stan_model("local_level.stan")
```


```r
nile_1_data <- within(list(), {
  y <- matrix(Nile_$flow)
  n <- nrow(y)
  a1 <- array(0, 1)
  P1 <- matrix(10 ^ 7)
  y_scale <- sd(Nile_$flow)
})
nile_1_samples <-
  sampling(local_level_mod,
           chains = 1,
           iter = 500,
           data = nile_1_data)
#> 
#> SAMPLING FOR MODEL 'local_level' NOW (CHAIN 1).
#> 
#> Chain 1, Iteration:   1 / 500 [  0%]  (Warmup)
#> Chain 1, Iteration:  50 / 500 [ 10%]  (Warmup)
#> Chain 1, Iteration: 100 / 500 [ 20%]  (Warmup)
#> Chain 1, Iteration: 150 / 500 [ 30%]  (Warmup)
#> Chain 1, Iteration: 200 / 500 [ 40%]  (Warmup)
#> Chain 1, Iteration: 250 / 500 [ 50%]  (Warmup)
#> Chain 1, Iteration: 251 / 500 [ 50%]  (Sampling)
#> Chain 1, Iteration: 300 / 500 [ 60%]  (Sampling)
#> Chain 1, Iteration: 350 / 500 [ 70%]  (Sampling)
#> Chain 1, Iteration: 400 / 500 [ 80%]  (Sampling)
#> Chain 1, Iteration: 450 / 500 [ 90%]  (Sampling)
#> Chain 1, Iteration: 500 / 500 [100%]  (Sampling)
#>  Elapsed Time: 6.77371 seconds (Warm-up)
#>                6.50536 seconds (Sampling)
#>                13.2791 seconds (Total)
```

Now, summarize the MCMC samples using the `summary` function on the `stanfit` object.
Additionally, I use the `tidy_stan_summary` function to make the results of `summary` easier to work with.
This converts the results of `summary` from a list of matrices to a list of data frames, and also parses the parameter names so that it is easier to select particular parameter values by name. I also will only use only the summary statistics for the combined chains.

```r
nile_1_summary <- tidy_stan_summary(summary(nile_1_samples))[["all"]] %>%
  left_join(Nile_, by = c("dim_1" = "obs"))
```
The estimated variances of the observation and state variances,

```r
filter(nile_1_summary, parameter %in% c("H", "Q")) %>%
  select(parname, mean, se_mean, p2.5, p97.5, n_eff, Rhat)
#> # A tibble: 2 × 7
#>   parname  mean se_mean  p2.5 p97.5 n_eff  Rhat
#>     <chr> <dbl>   <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1  H[1,1] 15031     301  9607 20772  79.3  1.02
#> 2  Q[1,1]  2025     130   433  4996 103.6  1.01
```
are similar to the MLE estimates producted by `StructTS`,

```r
StructTS(Nile_$flow, type = "level")
#> 
#> Call:
#> StructTS(x = Nile_$flow, type = "level")
#> 
#> Variances:
#>   level  epsilon  
#>    1469    15099
```
However, since the Bayesian estimates are means, the MLE estimates are modes,
and the posterior distribution of the variances are right skewed, the means are
larger than the posterior modes.



```r
str_keep <- function(string, pattern) {
  string[str_detect(string, pattern)]
}

ggplot(filter(nile_1_summary, parameter == "alpha"),
       aes(x = year,
           ymin = mean - 2 * sd,
           ymax = mean + 2 * sd)) +
  geom_ribbon(alpha = 0.3) +
  geom_line(aes(y = mean)) +
  geom_point(aes(y = flow)) +
  ylab("Annual river flow") +
  xlab("Observation") +
  theme_minimal()
```

<img src="example-Nile_files/figure-html/nile_1_states-1.png" width="672" />



**TODO** Diagnostics. What are the relevant Bayesian analogs?


### Local level with known intervention (intercept)


```r
nile_2_mod <- ssm_stan_model("local_level_reg.stan")
```


```r
nile_2_data <- nile_1_data
nile_2_data[["x"]] <- matrix(as.integer(Nile_$year > 1899))
nile_2_data[["k"]] <- ncol(nile_2_data[["x"]])
nile_2_samples <- sampling(nile_2_mod, chains = 1, iter = 500,
                           data = nile_2_data)
#> 
#> SAMPLING FOR MODEL 'local_level_reg' NOW (CHAIN 1).
#> 
#> Chain 1, Iteration:   1 / 500 [  0%]  (Warmup)
#> Chain 1, Iteration:  50 / 500 [ 10%]  (Warmup)
#> Chain 1, Iteration: 100 / 500 [ 20%]  (Warmup)
#> Chain 1, Iteration: 150 / 500 [ 30%]  (Warmup)
#> Chain 1, Iteration: 200 / 500 [ 40%]  (Warmup)
#> Chain 1, Iteration: 250 / 500 [ 50%]  (Warmup)
#> Chain 1, Iteration: 251 / 500 [ 50%]  (Sampling)
#> Chain 1, Iteration: 300 / 500 [ 60%]  (Sampling)
#> Chain 1, Iteration: 350 / 500 [ 70%]  (Sampling)
#> Chain 1, Iteration: 400 / 500 [ 80%]  (Sampling)
#> Chain 1, Iteration: 450 / 500 [ 90%]  (Sampling)
#> Chain 1, Iteration: 500 / 500 [100%]  (Sampling)
#>  Elapsed Time: 29.3301 seconds (Warm-up)
#>                4.77555 seconds (Sampling)
#>                34.1056 seconds (Total)
#> The following numerical problems occured the indicated number of times after warmup on chain 1
#>                                                                                                      count
#> Exception thrown at line 2525: Exception thrown at line 1071: Exception thrown at line 279: multiply    14
#> When a numerical problem occurs, the Hamiltonian proposal gets rejected.
#> See http://mc-stan.org/misc/warnings.html#exception-hamiltonian-proposal-rejected
#> If the number in the 'count' column is small, do not ask about this message on stan-users.
```


```r
nile_2_summary <- tidy_stan_summary(summary(nile_2_samples))[["all"]] %>%
  left_join(Nile_, by = c("dim_1" = "obs"))
```



```r
ggplot(filter(nile_2_summary, parameter == "mu"),
       aes(x = year,
           ymin = mean - 2 * sd,
           ymax = mean + 2 * sd)) +
  geom_ribbon(alpha = 0.3) +
  geom_line(aes(y = mean)) +
  geom_point(aes(y = flow)) +
  ylab("Annual river flow") +
  xlab("Observation") +
  theme_minimal()
```

<img src="example-Nile_files/figure-html/nile_2_states-1.png" width="672" />


### Local Level with known intervention (variance)


```r
nile_3_mod <- ssm_stan_model("local_level_interven.stan")
```


```r
nile_3_data <- nile_1_data
nile_3_data[["s"]] <- ifelse(Nile_$year == 1899, 10, 1)
nile_3_samples <- sampling(nile_3_mod, chains = 1, iter = 500,
                           data = nile_3_data)
#> 
#> SAMPLING FOR MODEL 'local_level_interven' NOW (CHAIN 1).
#> 
#> Chain 1, Iteration:   1 / 500 [  0%]  (Warmup)
#> Chain 1, Iteration:  50 / 500 [ 10%]  (Warmup)
#> Chain 1, Iteration: 100 / 500 [ 20%]  (Warmup)
#> Chain 1, Iteration: 150 / 500 [ 30%]  (Warmup)
#> Chain 1, Iteration: 200 / 500 [ 40%]  (Warmup)
#> Chain 1, Iteration: 250 / 500 [ 50%]  (Warmup)
#> Chain 1, Iteration: 251 / 500 [ 50%]  (Sampling)
#> Chain 1, Iteration: 300 / 500 [ 60%]  (Sampling)
#> Chain 1, Iteration: 350 / 500 [ 70%]  (Sampling)
#> Chain 1, Iteration: 400 / 500 [ 80%]  (Sampling)
#> Chain 1, Iteration: 450 / 500 [ 90%]  (Sampling)
#> Chain 1, Iteration: 500 / 500 [100%]  (Sampling)
#>  Elapsed Time: 5.4583 seconds (Warm-up)
#>                4.1302 seconds (Sampling)
#>                9.5885 seconds (Total)
#> The following numerical problems occured the indicated number of times after warmup on chain 1
#>                                                                                                      count
#> Exception thrown at line 2523: Exception thrown at line 1071: Exception thrown at line 279: multiply    14
#> When a numerical problem occurs, the Hamiltonian proposal gets rejected.
#> See http://mc-stan.org/misc/warnings.html#exception-hamiltonian-proposal-rejected
#> If the number in the 'count' column is small, do not ask about this message on stan-users.
```


```r
nile_3_summary <- tidy_stan_summary(summary(nile_3_samples))[["all"]] %>%
  left_join(Nile_, by = c("dim_1" = "obs"))
```



```r
ggplot(filter(nile_3_summary, parameter == "alpha"),
       aes(x = year,
           ymin = mean - 2 * sd,
           ymax = mean + 2 * sd)) +
  geom_ribbon(alpha = 0.3) +
  geom_line(aes(y = mean)) +
  geom_point(aes(y = flow)) +
  ylab("Annual river flow") +
  xlab("Observation") +
  theme_minimal()
```

<img src="example-Nile_files/figure-html/nile_3_states-1.png" width="672" />

### Local Level model with Sparse State Disturbances


```r
nile_4_mod <- ssm_stan_model("local_level_tvvar.stan")
```


```r
nile_4_data <- nile_1_data
nile_4_data[["s"]] <- 1 / nrow(Nile_)
nile_4_samples <- sampling(nile_4_mod, chains = 1, iter = 500,
                           data = nile_4_data)
#> 
#> SAMPLING FOR MODEL 'local_level_tvvar' NOW (CHAIN 1).
#> 
#> Chain 1, Iteration:   1 / 500 [  0%]  (Warmup)
#> Chain 1, Iteration:  50 / 500 [ 10%]  (Warmup)
#> Chain 1, Iteration: 100 / 500 [ 20%]  (Warmup)
#> Chain 1, Iteration: 150 / 500 [ 30%]  (Warmup)
#> Chain 1, Iteration: 200 / 500 [ 40%]  (Warmup)
#> Chain 1, Iteration: 250 / 500 [ 50%]  (Warmup)
#> Chain 1, Iteration: 251 / 500 [ 50%]  (Sampling)
#> Chain 1, Iteration: 300 / 500 [ 60%]  (Sampling)
#> Chain 1, Iteration: 350 / 500 [ 70%]  (Sampling)
#> Chain 1, Iteration: 400 / 500 [ 80%]  (Sampling)
#> Chain 1, Iteration: 450 / 500 [ 90%]  (Sampling)
#> Chain 1, Iteration: 500 / 500 [100%]  (Sampling)
#>  Elapsed Time: 31.9269 seconds (Warm-up)
#>                12.5386 seconds (Sampling)
#>                44.4654 seconds (Total)
#> The following numerical problems occured the indicated number of times after warmup on chain 1
#>                                                                                                      count
#> Exception thrown at line 2524: Exception thrown at line 1071: Exception thrown at line 279: multiply    11
#> When a numerical problem occurs, the Hamiltonian proposal gets rejected.
#> See http://mc-stan.org/misc/warnings.html#exception-hamiltonian-proposal-rejected
#> If the number in the 'count' column is small, do not ask about this message on stan-users.
```


```r
nile_4_summary <- tidy_stan_summary(summary(nile_4_samples))[["all"]] %>%
  left_join(Nile_, by = c("dim_1" = "obs"))
```



```r
ggplot(filter(nile_4_summary, parameter == "alpha"),
       aes(x = year,
           ymin = mean - 2 * sd,
           ymax = mean + 2 * sd)) +
  geom_ribbon(alpha = 0.3) +
  geom_line(aes(y = mean)) +
  geom_point(aes(y = flow)) +
  ylab("Annual river flow") +
  xlab("Observation") +
  theme_minimal()
```

<img src="example-Nile_files/figure-html/nile_4_states-1.png" width="672" />
