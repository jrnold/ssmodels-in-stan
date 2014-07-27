library("dplyr")
library("boot")
library("rstan")
data("coal")

fill_na <- function(x, fill = 0) {
    x[is.na(x)] <- fill
    x
}

allyears <- data.frame(year = seq(1851, 1962, by = 1))
coal_by_year <-
    (mutate(coal,
            year = floor(date))
     %>% group_by(year)
     %>% summarise(disasters = length(year))
     %>% merge(allyears, all.y = TRUE)
     %>% mutate(disasters = fill_na(disasters)))

m <- stan_model("stan/chib_coal.stan")

# Priors from Chib (98)
# 3 state model
# mu ~ Gamma(3, 1)
# p_{11}, p_{22} ~ Beta(5, 1)
m1_data <-
    within(list(), {
        y <- coal_by_year$disasters
        n <- length(y)
        m <- 2
        mu_a <- 2
        mu_b <- 1
        rho_a <- 5
        rho_b <- 1
    })

m1_est <- sampling(m, data = m1_data, chains = 1)

m1_data <-
    within(list(), {
        y <- coal_by_year$disasters
        n <- length(y)
        m <- 3
        mu_a <- 3
        mu_b <- 1
        rho_a <- 5
        rho_b <- 1
    })

m1_est <- sampling(m, data = m1_data, chains = 1)


ret <- m1_est
Mode <- function(x) {
    tbl <- table(x)
    as.integer(names(tbl)[which.max(tbl)])
}

stateprob <- extract(ret, "stateprob")[[1]]
extract(ret, "viterbi")[[1]][1, ]
apply(extract(ret, "p")[[1]], 2:3, mean)
apply(extract(ret, "mu")[[1]], 2, median)
