remove(list = ls())
library(keras)
library(vrtest)
library("readxl")
setwd("C:/Users/jbastos/1_research/99_GitHub/martingales/code")
source("TS_gen.R")

convnet = function(size) {
  model = keras_model_sequential() 
  model %>% 
    layer_conv_1d(filters = 32, kernel_size = 12,  activation = "relu",  input_shape=c(size, 1)) %>%
    layer_max_pooling_1d(pool_size = 4) %>%
    layer_flatten() %>% 
    layer_dense(units = 16, activation = 'relu') %>%
    layer_dropout(rate = 0.25) %>% 
    layer_dense(units = 1, activation = 'sigmoid')
    model %>% compile(loss = 'binary_crossentropy',
                    optimizer = optimizer_adam(lr = 0.00001),
                    metrics = c('accuracy'))
}

# ==================================================================================================
#    Wright (2001) data
# ==================================================================================================

data(exrates)
M = nrow(exrates) - 1
N = 100000

H0 = matrix(0, nrow = N, ncol = M)
H1 = matrix(0, nrow = N, ncol = M)
for(i in 1:N) {
  H0[i,] = AR(size = M, phi = 0, omega = 0, alpha = 0, beta = 1, st = FALSE, inov = 1)
  H1[i,] = AR(size = M, phi = runif(1, 0.01, 0.15), omega = 0, alpha = 0, beta = 1, st = FALSE, inov = 1)
}
xT = rbind(H0[1:N,], H1[1:N,])
yT = c(rep(0, N), rep(1, N))
xT  = array_reshape(xT, c(dim(xT), 1))
m = convnet(M)
history = m %>% fit(xT, yT, batch_size = 32, epochs = 50)
#save_model_hdf5(m, file = "convnets/Wright2001.h5")

# Scores, critical values and p-values

m = load_model_hdf5(file = "../convnets/Wright2001.h5")

H0 = matrix(0, nrow = N, ncol = M)
for(i in 1:N) {
  H0[i,] = AR(size = M, phi = 0, omega = 0, alpha = 0, beta = 1, st = FALSE, inov = 1)
}
H0_dist = m %>% predict(array_reshape(H0, c(dim(H0), 1)))

for(fx in (names(exrates))) {
  y = exrates[,fx]
  r = log(y[2:length(y)]) - log(y[1:(length(y)-1)])
  r_convnet = array_reshape(scale(r), c(1,M,1))
  score  = m %>% predict(r_convnet, verbose = 0)
  pvalue = sum(H0_dist > rep(score, N)) / N
  cat(fx, score, pvalue, AutoBoot.test(r, nboot = 5000, wild = "Normal")$pval, "\n")
}

# ==================================================================================================
#    Federal reserve exchange rates
# ==================================================================================================

M = 1095 # weekly data
#M = 251  # monthly data
N = 100000
H0 = matrix(0, nrow = N, ncol = M)
H1 = matrix(0, nrow = N, ncol = M)

for(i in 1:N) {
  H0[i,] = AR(size = M, phi = 0, omega = 0, alpha = 0, beta = 1, st = FALSE, inov = 1)
  H1[i,] = AR(size = M, phi = runif(1, 0.01, 0.15), omega = 0, alpha = 0, beta = 1, st = FALSE, inov = 1)
}
xT = rbind(H0[1:N,], H1[1:N,])
yT = c(rep(0, N), rep(1, N))
xT  = array_reshape(xT, c(dim(xT), 1))
m = convnet(M)
history = m %>% fit(xT, yT, batch_size = 32, epochs = 50)

#save_model_hdf5(m, file = "convnets/FRNY_weekly.h5")
#save_model_hdf5(m, file = "convnets/FRNY_monthly.h5")

# Scores, critical values and p-values

#m = load_model_hdf5(file = "../convnets/FRNY_weekly.h5")
m = load_model_hdf5(file = "../convnets/FRNY_monthly.h5")

#M = 1095 # weekly data
M = 251  # monthly data
N = 100000

H0 = matrix(0, nrow = N, ncol = M)
for(i in 1:N) {
  H0[i,] = AR(size = M, phi = 0, omega = 0, alpha = 0, beta = 1, st = FALSE, inov = 1)
}
H0_dist = m %>% predict(array_reshape(H0, c(dim(H0), 1)))

#df = data.frame(read_excel("../data/FX_rates.xlsx", sheet = "weekly", col_names = TRUE, na = "NA"))
df = data.frame(read_excel("../data/FX_rates.xlsx", sheet = "monthly", col_names = TRUE, na = "NA"))

for(fx in names(df[-1])) {
  y = df[,fx]
  r = log(y[2:length(y)]) - log(y[1:(length(y)-1)])
  r = r[(length(r) - M + 1):length(r)]
  r_convnet = array_reshape(scale(r), c(1,M,1))
  score  = m %>% predict(r_convnet, verbose = 0)
  pvalue = sum(H0_dist > rep(score, N)) / N
  cat(fx, score, pvalue, AutoBoot.test(r, nboot = 5000, wild = "Normal")$pval, "\n")
}

