
## CAPM

This is a simplified example of CAPM model from @Petris2010a.
A simplified capital-asset pricicing model (CAPM) model is a model of the excess
return (return relative to a risk-free asset) of stocks. In this example,
the excess return of a stock is proportional to the market return (return of a portfolio of stocks representing the whole market).


Let $y_{j,t}$ be the excess return of stock $j$ in period $t$.
The model that we will estimate is,
$$
\begin{aligned}[t]
\vec{y}_{j,t} &= \beta_{j,t} x_t + \varepsilon_{j,t} & \vec{\varepsilon}_t & \sim N(\vec{0}, \mat{\Sigma}_{\varepsilon} \\
\beta_{j,t} &= \beta_{j,t} + \eta_{j,t} & \vec{\eta}_{t} & \sim N(\vec{0}, \mat{\Sigma}_{\eta})
\end{aligned}
$$
The proportionality constants of the stocks are stock- and time-specific, with
each $\vec{\beta}_j$ following a random walk.[^alpha]
However, both the observation disturbances ($\mat{\Sigma}_{\varepsilon}$ and state disturbances ($\mat{\Sigma}_{\eta}$) allow for correlation between the stocks.
In the model estimated here, the intercepts for each stock are assumed to be zero.

The data, originally from @Berndt1991a, are the monthly returns of four stocks (Mobil, IBM, Weyer, and Citicorp) from January 1978 to December 1987.
The interest rate on 30-day Treasury Bill is used as the risk free rate (`rkfree`), and the value-weighted average returns of all stocks listed on the New York and American Stock Exchanges as the market return (`market`).

```r
data("capm", package = "StanStateSpace")
glimpse(capm)
#> Observations: 120
#> Variables: 7
#> $ date   <date> 1978-01-01, 1978-02-01, 1978-03-01, 1978-04-01, 1978-0...
#> $ mobil  <dbl> -0.046, -0.017, 0.049, 0.077, -0.011, -0.043, 0.028, 0....
#> $ ibm    <dbl> -0.029, -0.043, -0.063, 0.130, -0.018, -0.004, 0.092, 0...
#> $ weyer  <dbl> -0.116, -0.135, 0.084, 0.144, -0.031, 0.005, 0.164, 0.0...
#> $ citcrp <dbl> -0.115, -0.019, 0.059, 0.127, 0.005, 0.007, 0.032, 0.08...
#> $ market <dbl> -0.045, 0.010, 0.050, 0.063, 0.067, 0.007, 0.071, 0.079...
#> $ rkfree <dbl> 0.00487, 0.00494, 0.00526, 0.00491, 0.00513, 0.00527, 0...
```
In the analysis, the excess returns (return minus risk-free return) are used.
So, the first step is to subtract the risk-free returns (`rkfree`) from each of the stock return columns and the market return.

```r
capm <- capm %>%
  mutate_each(funs(. - rkfree), one_of(c("mobil", "ibm", "weyer", "citcrp",
                                         "market"))) %>%
  mutate(obs = row_number())
ggplot(gather(select(capm, -rkfree), asset, rate, -date, -obs),
       aes(x = date, y = rate)) +
  geom_line() +
  facet_grid(asset ~ .) +
  theme_minimal()
```

<img src="example-capm_files/figure-html/capm-data-clean-1.png" width="672" />



```r
capm_mod <- ssm_stan_model("capm.stan")

capm_data <-
  within(list(), {
    y <- as.matrix(capm[ , c("mobil", "ibm", "weyer", "citcrp")])
    x <- capm[["market"]]
    n <- nrow(capm)
    m <- ncol(y)
    y_scale <- apply(y, 2, sd)
    a1 <- rep(0, m)
    P1 <- diag(1e6, m, m)
  })

capm_samples <- sampling(capm_mod, chains = 1, iter = 1000,
                         data = capm_data)
```


```r
capm_samples <- tidy_stan_summary(summary(capm_samples))[["all"]] %>%
  left_join(mutate(capm, obs = row_number()) %>% select(date, obs),
            by = c("dim_1" = "obs"))

capm_betas <- filter(capm_samples, parameter == "beta") %>%
         mutate(stock = recode(dim_2,
                                `1` = "mobil",
                                `2` = "ibm",
                                `3` = "weyer",
                                `4` = "citcrp")) %>%
  left_join(gather(capm, stock, rate, -market, -date, -obs),
            by = c("date", "stock"))

ggplot(capm_betas,
       aes(x = date,
           ymin = mean - 2 * sd,
           ymax = mean + 2 * sd)) +
  geom_ribbon(alpha = 0.3) +
  geom_line(aes(y = mean)) +
  facet_grid(stock ~ .) +
  theme_minimal()
```

<img src="example-capm_files/figure-html/capm-state-plot-1.png" width="672" />

**TODO** Show estimated covariance matrices
