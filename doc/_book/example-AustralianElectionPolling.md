
## Australian Election Polling

Example from Simon Jackman's textbook (Ex 9.3) and 2005 paper.

The data consist of 239 polls by the Australian pollsters (Newspoll, Nielsen, Morgan, and Galaxy) betwen the October 9, 2004, and November 24, 2007 Australian Federal elections.

```r
data("AustralianElectionPolling", package = "pscl")
glimpse(AustralianElectionPolling)
#> Observations: 239
#> Variables: 14
#> $ ALP         <dbl> 39.5, 39.0, 38.0, 36.0, 33.0, 36.5, 39.0, 37.0, 34...
#> $ Lib         <dbl> 44.5, 44.0, 46.0, 46.5, 47.0, 45.5, 46.0, 47.0, 52...
#> $ Nat         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,...
#> $ Green       <dbl> 8.5, 8.5, 6.0, 9.0, 8.0, 9.5, 6.0, 7.5, 8.0, 7.0, ...
#> $ FamilyFirst <dbl> 2.0, 1.5, 0.0, 2.5, 0.0, 2.0, 0.0, 2.0, 0.0, 0.0, ...
#> $ Dems        <dbl> 2.0, 2.0, 0.0, 1.5, 0.0, 1.5, 0.0, 1.5, 2.0, 0.0, ...
#> $ OneNation   <dbl> 1.0, 1.0, 0.0, 1.0, 0.0, 1.5, 0.0, 1.0, 1.0, 0.0, ...
#> $ DK          <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,...
#> $ sampleSize  <dbl> 1451, 2090, 1150, 1451, 1130, 2136, 1132, 2010, 14...
#> $ org         <fctr> Morgan, F2F, Morgan, F2F, Newspoll, Morgan, F2F, ...
#> $ startDate   <date> 2004-10-30, 2004-11-13, 2004-11-19, 2004-11-27, 2...
#> $ endDate     <date> 2004-11-07, 2004-11-21, 2004-11-21, 2004-12-05, 2...
#> $ source      <chr> "", "http://www.roymorgan.com/news/polls/2004/3808...
#> $ remark      <chr> "", "face-to-face", "", "face-to-face", "", "face-...
```

We will estimate the the proportion of the population that will vote the Australian Labor Party (ALP) in a House of Representatives election for each day between the elections. 
We model this as 

$$
\begin{aligned}[t]
y_{j,t} &= \delta_{j} + \xi_t + \varepsilon_{j,t} & \varepsilon_{j,t} & \sim N\left(0, \frac{y_{j,t} (1 - y_{j,t})}{r_{j,t}} \right) \\
\xi_t & \xi_{t - 1} + \omega_t & \omega_{t} & \sim N(0, \sigma_{\omega})
\end{aligned}
$$
The observation variances are treated as known, and derived from the sample size and results of the poll.
Jackman 2012 (p. 475) puts a prior on the state variance to ensure that the daily changes in thestate are not too large, with a 95% interval of $\pm 2%$, implying $\sigma_{\omega} = 0.01$.

Polls have between 500 to 2,600 respondents, with a median of 1,156.
While most of the polls use complex survey methods such as post-stratification weights to adjust for non-response, for the purposes of this analysis we will use the stated sample size of the poll as its effective size, and use that to construct the measurement error. There are 1,142 days in the period under consideration.

The unique polling firms and number of polls are,

```r
AustralianElectionPolling %>% group_by(org) %>% tally()
#> # A tibble: 5 x 2
#>             org     n
#>          <fctr> <int>
#> 1        Galaxy    10
#> 2   Morgan, F2F    94
#> 3      Newspoll    78
#> 4       Nielsen    42
#> 5 Morgan, Phone    15
```
The face-to-face and phone polls by Morgan are treated separately.
Each polling firm will be given it's own "house effect" term to allow them to 
be systematically biased.
However, we will require that the average of the house effects is zero; 
meaning that while each polling firm may be biased, on average they are not.
$$
\delta_j \sim N(0, d)
$$
Jackman uses a fixed value of $d$ that gives a 95% interval spanning from -0.15 to .15. This corresponds to $d = 0.075$.

In the 2004 election the Australian Labor Party (ALP) won 37.6% percent of
first preferences, and in 2007 they won 43.4%.
These will be treated as coming from a source with no house effect and an infinitely large sample size.


```raustralianelectionpolling
library("dplyr")
data("AustralianElectionPolling", package = "pscl")
glimpse(AustralianElectionPolling)
```

Unique polling organizations (Morgan face-to-face and phone polls treated separately):

```r
AustralianElectionPolling %>% group_by(org) %>% tally()
#> # A tibble: 5 x 2
#>             org     n
#>          <fctr> <int>
#> 1        Galaxy    10
#> 2   Morgan, F2F    94
#> 3      Newspoll    78
#> 4       Nielsen    42
#> 5 Morgan, Phone    15
```

Date of the poll is the median day (rounded down) for when they are in the field.


```r
aus_elections <-
  bind_rows(
    data_frame(date = as.Date("2004-10-09"), ALP = 37.6, sd = .00025),
    data_frame(date = as.Date("2007-11-24"), ALP = 43.4, sd = .00025)
  )

aus_polling <-
  AustralianElectionPolling %>%
  group_by(org, startDate, endDate) %>%
  # There are two Nielson polls on 2007-11-19 to 21.
  summarize(ALP = weighted.mean(ALP, sampleSize),
            sampleSize = sum(sampleSize)) %>%
  ungroup() %>%
  mutate(date = floor_date(as.Date(startDate) + difftime(endDate, startDate, unit = "days"), unit = "day"),
         y = ALP / 100,
         sd = sqrt(y * (1 - y) / sampleSize),
         org = as.character(org)) %>%
  select(date, org, y, sd) %>%
  bind_rows(select(mutate(aus_elections, y = ALP / 100, org = "Election"), -ALP)) %>%
  mutate(time = as.integer(difftime(date, as.Date("2004-10-09"), unit = "days")) + 1)

datelist <-
  data_frame(date = seq(min(aus_polling[["date"]]), max(aus_polling[["date"]]), by = "days"),
             time = min(aus_polling[["time"]]):max(aus_polling[["time"]]))

polling_values <-
  aus_polling %>%
  select(date, time, org, y) %>%
  spread(org, y) %>%
  full_join(select(datelist, date), by = "date") %>%
  arrange(date) %>%
  select(-time)

polling_sd <-
  aus_polling %>%
  select(date, time, org, sd) %>%
  spread(org, sd) %>%
  full_join(select(datelist, date), by = c("date")) %>%
  arrange(date) %>%  
  select(-time)


aus_elec_data <-
  within(list(), {
    y <- as.matrix(select(polling_values, -date))
    missing <- ssm_process_na(y)
    p_t <- missing$n
    y_idx <- missing$idx
    rm(missing)
    y[is.na(y)] <- 0
    n <- nrow(y)
    p <- ncol(y)
    sigma_eps <- as.matrix(select(polling_sd, -date)) %>%
      apply(2, function(.) if_else(is.na(.), median(., na.rm = TRUE), .))
    sigma_delta <- .0075 # prior on house effects. p. 479. Jackman text.
    zeta <- 0.1
    a1 <- array(0.5)
    P1 <- matrix(0.25, 1, 1)
  })
```


```r
polling_mod <- ssm_stan_model("polling.stan")
#> In file included from filee49614c5cd4f.cpp:8:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/stan/math/rev/core.hpp:42:
#> /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/stan/math/rev/core/set_zero_all_adjoints.hpp:14:17: warning: unused function 'set_zero_all_adjoints' [-Wunused-function]
#>     static void set_zero_all_adjoints() {
#>                 ^
#> In file included from filee49614c5cd4f.cpp:8:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/stan/math/rev/core.hpp:43:
#> /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/stan/math/rev/core/set_zero_all_adjoints_nested.hpp:17:17: warning: 'static' function 'set_zero_all_adjoints_nested' declared in header file should be declared 'static inline' [-Wunneeded-internal-declaration]
#>     static void set_zero_all_adjoints_nested() {
#>                 ^
#> In file included from filee49614c5cd4f.cpp:8:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/stan/math/rev/mat.hpp:8:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/stan/math/prim/mat.hpp:55:
#> /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/stan/math/prim/mat/fun/autocorrelation.hpp:19:14: warning: function 'fft_next_good_size' is not needed and will not be emitted [-Wunneeded-internal-declaration]
#>       size_t fft_next_good_size(size_t N) {
#>              ^
#> In file included from filee49614c5cd4f.cpp:8:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/stan/math/rev/mat.hpp:8:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/stan/math/prim/mat.hpp:36:
#> /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/stan/math/prim/mat/err/check_positive_ordered.hpp:39:67: warning: unused typedef 'size_type' [-Wunused-local-typedef]
#>       typedef typename index_type<Matrix<T_y, Dynamic, 1> >::type size_type;
#>                                                                   ^
#> In file included from filee49614c5cd4f.cpp:8:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/stan/math/rev/mat.hpp:8:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/stan/math/prim/mat.hpp:232:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/stan/math/prim/arr.hpp:33:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/StanHeaders/include/stan/math/prim/arr/functor/integrate_ode_rk45.hpp:13:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/BH/include/boost/numeric/odeint.hpp:61:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/BH/include/boost/numeric/odeint/util/multi_array_adaption.hpp:29:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/BH/include/boost/multi_array.hpp:21:
#> In file included from /Library/Frameworks/R.framework/Versions/3.3/Resources/library/BH/include/boost/multi_array/base.hpp:28:
#> /Library/Frameworks/R.framework/Versions/3.3/Resources/library/BH/include/boost/multi_array/concept_checks.hpp:42:43: warning: unused typedef 'index_range' [-Wunused-local-typedef]
#>       typedef typename Array::index_range index_range;
#>                                           ^
#> /Library/Frameworks/R.framework/Versions/3.3/Resources/library/BH/include/boost/multi_array/concept_checks.hpp:43:37: warning: unused typedef 'index' [-Wunused-local-typedef]
#>       typedef typename Array::index index;
#>                                     ^
#> /Library/Frameworks/R.framework/Versions/3.3/Resources/library/BH/include/boost/multi_array/concept_checks.hpp:53:43: warning: unused typedef 'index_range' [-Wunused-local-typedef]
#>       typedef typename Array::index_range index_range;
#>                                           ^
#> /Library/Frameworks/R.framework/Versions/3.3/Resources/library/BH/include/boost/multi_array/concept_checks.hpp:54:37: warning: unused typedef 'index' [-Wunused-local-typedef]
#>       typedef typename Array::index index;
#>                                     ^
#> filee49614c5cd4f.cpp:166:20: warning: unused typedef 'fun_scalar_t__' [-Wunused-local-typedef]
#>     typedef double fun_scalar_t__;
#>                    ^
#> filee49614c5cd4f.cpp:204:20: warning: unused typedef 'fun_scalar_t__' [-Wunused-local-typedef]
#>     typedef double fun_scalar_t__;
#>                    ^
#> filee49614c5cd4f.cpp:633:20: warning: unused typedef 'fun_scalar_t__' [-Wunused-local-typedef]
#>     typedef double fun_scalar_t__;
#>                    ^
#> filee49614c5cd4f.cpp:674:20: warning: unused typedef 'fun_scalar_t__' [-Wunused-local-typedef]
#>     typedef double fun_scalar_t__;
#>                    ^
#> filee49614c5cd4f.cpp:716:20: warning: unused typedef 'fun_scalar_t__' [-Wunused-local-typedef]
#>     typedef double fun_scalar_t__;
#>                    ^
#> filee49614c5cd4f.cpp:765:20: warning: unused typedef 'fun_scalar_t__' [-Wunused-local-typedef]
#>     typedef double fun_scalar_t__;
#>                    ^
#> filee49614c5cd4f.cpp:1651:20: warning: unused typedef 'fun_scalar_t__' [-Wunused-local-typedef]
#>     typedef double fun_scalar_t__;
#>                    ^
#> filee49614c5cd4f.cpp:1706:20: warning: unused typedef 'fun_scalar_t__' [-Wunused-local-typedef]
#>     typedef double fun_scalar_t__;
#>                    ^
#> filee49614c5cd4f.cpp:2446:20: warning: unused typedef 'fun_scalar_t__' [-Wunused-local-typedef]
#>     typedef double fun_scalar_t__;
#>                    ^
#> filee49614c5cd4f.cpp:3498:20: warning: unused typedef 'fun_scalar_t__' [-Wunused-local-typedef]
#>     typedef double fun_scalar_t__;
#>                    ^
#> filee49614c5cd4f.cpp:3757:20: warning: unused typedef 'fun_scalar_t__' [-Wunused-local-typedef]
#>     typedef double fun_scalar_t__;
#>                    ^
#> filee49614c5cd4f.cpp:4021:20: warning: unused typedef 'fun_scalar_t__' [-Wunused-local-typedef]
#>     typedef double fun_scalar_t__;
#>                    ^
#> filee49614c5cd4f.cpp:4502:20: warning: unused typedef 'fun_scalar_t__' [-Wunused-local-typedef]
#>     typedef double fun_scalar_t__;
#>                    ^
#> filee49614c5cd4f.cpp:4555:20: warning: unused typedef 'fun_scalar_t__' [-Wunused-local-typedef]
#>     typedef double fun_scalar_t__;
#>                    ^
#> 22 warnings generated.
polling_samples <-
  sampling(polling_mod, data = aus_elec_data, iter = 2000, chains = 1,
           init = list(list(delta = rep(0, aus_elec_data[["p"]] - 1),
                            sigma_eta = 0.005,
                            labmda = rep(1 / aus_elec_data[["n"]], aus_elec_data[["n"]]))))
#> 
#> SAMPLING FOR MODEL 'polling' NOW (CHAIN 1).
#> 
#> Chain 1, Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1, Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1, Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1, Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1, Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1, Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1, Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1, Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1, Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1, Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1, Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1, Iteration: 2000 / 2000 [100%]  (Sampling)
#>  Elapsed Time: 211.938 seconds (Warm-up)
#>                167.657 seconds (Sampling)
#>                379.595 seconds (Total)
polling_summary <- tidy_stan_summary(summary(polling_samples))
```


```r
ggplot(aus_polling, aes(x = date, y = y,
                        ymin = y - 2 * sd,
                        ymax = y + 2 * sd,
                        colour = org)) +
  geom_pointrange() +
  theme_minimal()
```

<img src="example-AustralianElectionPolling_files/figure-html/unnamed-chunk-5-1.png" width="672" />


```r
polling_samples <- polling_summary[["all"]] %>%
  left_join(datelist, by = c("dim_1" = "time"))

ggplot() +
  geom_ribbon(data = filter(polling_samples, parameter == "alpha"),
              aes(x = date, ymin = mean - 2 * sd, ymax = mean + 2 * sd), alpha = 0.3) +
  geom_line(data = filter(polling_samples, parameter == "alpha"),
            aes(x = date, y = mean)) +
  geom_point(data = aus_polling, aes(x = date, y = y, colour = org)) +
  theme_minimal()
```

<img src="example-AustralianElectionPolling_files/figure-html/capm-state-plot-1.png" width="672" />


```r
filter(polling_summary[["all"]], parameter == "delta") %>%
  select(parname, mean, se_mean, p2.5, p97.5, n_eff, Rhat)
#> # A tibble: 5 x 7
#>    parname      mean  se_mean     p2.5    p97.5 n_eff  Rhat
#>      <chr>     <dbl>    <dbl>    <dbl>    <dbl> <dbl> <dbl>
#> 1 delta[1] -0.013136 0.000181 -0.02256 -0.00477   690 1.000
#> 2 delta[2]  0.018446 0.000179  0.01123  0.02576   425 1.001
#> 3 delta[3]  0.000314 0.000199 -0.00841  0.00898   506 1.000
#> 4 delta[4]  0.004033 0.000169 -0.00323  0.01107   462 0.999
#> 5 delta[5]  0.002159 0.000173 -0.00521  0.00938   480 1.000
```

