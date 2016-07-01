training_data_list = dget("training_data_list.txt")

fit.stan2 <- stan (file = "linear_trend_ar2_model.stan",
                   data = list (N = length (training_data_list),
                                y = matrix (c (training_data_list), nrow = 1)),
                   # pars = c ( "V", "W"),
                   chains = 4,
                   iter = 2000, warmup = 1000, thin = 1)