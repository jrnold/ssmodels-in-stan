library("rstan")
earthquakes <- read.delim("../data/hmm-with-r/earthquakes.txt", header = FALSE)
names(earthquakes) <- c("year", "quakes")

m <- stan_model("hmm.stan")
standata <-
    within(list(), {
        y <- earthquakes[["quakes"]]
        n <- length(y)
        m <- 3
    })

Gamma <- matrix(0.01, standata$m, standata$m)
Gamma <- Gamma + diag(1 - rowSums(Gamma))
init <-
    list(list(Gamma = Gamma,
              lambda = quantile(standata$y,
                  prob = seq(0, 1, length = standata$m + 2)[2:(standata$m + 1)])))
gammamm <- function(x) {
    xmean <- mean(x)
    xvar <- var(x)
    beta <- (xmean / xvar)^(1/3)
    alpha <- xmean / beta
    list(alpha = alpha, beta = beta)
}
mm <- gammamm(init[[1]][["lambda"]])
init[[1]][["lambda_a"]] <- mm$alpha
init[[1]][["lambda_b"]] <- mm$beta

ret <- optimizing(m, data = standata, init = init[[1]])
ret$par[grepl("lambda", names(ret$par))]
ret$par[grepl("delta", names(ret$par))]

ret <- sampling(m, data = standata, chains=1, iter = 5000,
                init = init)

lambda <- apply(extract(ret, "lambda")[[1]], 2, mean)
print(lambda)

apply(extract(ret, "delta")[[1]], 2, median)
median(extract(ret, "llik")[[1]])
apply(extract(ret, "Gamma")[[1]], 2:3, median)

Mode <- function(x) {
    tbl <- table(x)
    as.integer(names(tbl)[which.max(tbl)])
}

stateprob <- extract(ret, "stateprob")[[1]]
maxstate <-
    apply(apply(stateprob, 1, function(x) apply(x, 1, which.max)),
          1, Mode)

str(extract(ret, "states")[[1]])


