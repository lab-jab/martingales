remove(list = ls())
library(keras)
library(vrtest)
library("readxl")

setwd("path_to/code")

source("TS_gen.R")

# ==================================================================================================
#    Auxiliary functions
# ==================================================================================================

np_crit_value <- function(scores, alpha, delta) {
  size = length(scores)
  violation_rates <- pbinom(q=0:(size-1), size = size, prob = 1 - alpha, lower.tail = FALSE)
  np_order = which(violation_rates <= delta)[1]
  return(sort(scores)[np_order])
}

convnet <- function(size) {
  model <- keras_model_sequential() 
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
M <- nrow(exrates) - 1
N <- 200000

H0 <- matrix(0, nrow = N, ncol = M)
H1 <- matrix(0, nrow = N, ncol = M)
for(i in 1:N) {
  H0[i,] <- AR(size = M, phi = 0, omega = runif(1, 0., 1E-3),
               alpha = runif(1, 0.85, 0.99), beta = 0, st = FALSE, inov = 1)
  H1[i,] <- AR(size = M, phi = runif(1, 0.01, 0.15), omega = runif(1, 0., 1E-3),
               alpha = runif(1, 0.85, 0.99), beta = 0, st = FALSE, inov = 1)
}
xT <- rbind(H0[1:N,], H1[1:N,])
yT <- c(rep(0, N), rep(1, N))
xT <- array_reshape(xT, c(dim(xT), 1))

m <- convnet(M)
history <- m %>% fit(xT, yT, batch_size = 32, epochs = 100)
save_model_hdf5(m, file = "convnets/Wright2001.h5")

# Calculate distribution for null hypothesis
m <- load_model_hdf5(file = "convnets/Wright2001.h5")
H0 <- matrix(0, nrow = N, ncol = M)
for(i in 1:N) {
  H0[i,] <- AR(size = M, phi = 0, omega = runif(1, 0., 1E-3),
               alpha = runif(1, 0.85, 0.99), beta = 0, st = FALSE, inov = 1)
}
scores_H0 <- m %>% predict(array_reshape(H0, c(dim(H0), 1)))

# Calculate critical values
crit_value_010 <- np_crit_value(scores_H0, 0.10, 0.1)
crit_value_005 <- np_crit_value(scores_H0, 0.05, 0.1)
crit_value_001 <- np_crit_value(scores_H0, 0.01, 0.1)
cat(crit_value_010, crit_value_005, crit_value_001)

# Calculate scores (test statistic)
data(exrates)
for(fx in names(exrates)) {
  y <- exrates[,fx]
  r <- log(y[2:length(y)]) - log(y[1:(length(y)-1)])
  r <- r[(length(r) - M + 1):length(r)]
  r_convnet <- array_reshape(scale(r), c(1,M,1))
  score <- m %>% predict(r_convnet, verbose = 0)
  cat(fx, score, "\n")
}

# Calculate p-values for alternative tests
for(fx in (names(exrates))) {
  y <- exrates[,fx]
  r <- log(y[2:length(y)]) - log(y[1:(length(y)-1)])
  r_convnet <- array_reshape(scale(r), c(1,M,1))
  pvalue_AVR <- AutoBoot.test(r, nboot = 500, wild = "Normal")$pval
  pvalue_AQ <- Auto.Q(r, lags = 10)$Pvalue
  cat(fx, pvalue_AVR, pvalue_AQ, "\n")
}

# ==================================================================================================
#    Federal reserve exchange rates
# ==================================================================================================

N <- 200000
M <- 1095
model_name = "FRNY_weekly.h5"
data_type = "monthly" # "weekly" or monthly"
if(data_type == "monthly") {
  M <- 251
  model_name = "FRNY_monthly.h5"
}
df <- data.frame(read_excel("data/FX_rates.xlsx", sheet = data_type, col_names = TRUE, na = "NA"))

H0 <- matrix(0, nrow = N, ncol = M)
H1 <- matrix(0, nrow = N, ncol = M)
for(i in 1:N) {
  H0[i,] <- AR(size = M, phi = 0, omega = runif(1, 0., 1E-3),
               alpha = runif(1, 0.85, 0.99), beta = 0, st = FALSE, inov = 1)
  H1[i,] <- AR(size = M, phi = runif(1, 0.01, 0.15), omega = runif(1, 0., 1E-3),
               alpha = runif(1, 0.85, 0.99), beta = 0, st = FALSE, inov = 1)
}
xT <- rbind(H0[1:N,], H1[1:N,])
yT <- c(rep(0, N), rep(1, N))
xT <- array_reshape(xT, c(dim(xT), 1))
m <- convnet(M)
history <- m %>% fit(xT, yT, batch_size = 32, epochs = 100)
save_model_hdf5(m, file = paste0("convnets/",model_name))

# Calculate distribution for null hypothesis
m <- load_model_hdf5(file = paste0("convnets/",model_name)) 
H0 <- matrix(0, nrow = N, ncol = M)
for(i in 1:N) {
  H0[i,] <- AR(size = M, phi = 0, omega = runif(1, 0.,1E-3),
               alpha = runif(1, 0.85, 0.99), beta = 0, st = FALSE, inov = 1)
}
scores_H0 <- m %>% predict(array_reshape(H0, c(dim(H0), 1)))

# Calculate critical values
crit_value_010 <- np_crit_value(scores_H0, 0.10, 0.1)
crit_value_005 <- np_crit_value(scores_H0, 0.05, 0.1)
crit_value_001 <- np_crit_value(scores_H0, 0.01, 0.1)
cat(crit_value_010, crit_value_005, crit_value_001)

# Calculate scores (test statistic)
for(fx in names(df[-1])) {
  y <- df[,fx]
  r <- log(y[2:length(y)]) - log(y[1:(length(y)-1)])
  r <- r[(length(r) - M + 1):length(r)]
  r_convnet <- array_reshape(scale(r), c(1,M,1))
  score <- m %>% predict(r_convnet, verbose = 0)
  cat(fx, score, "\n")
}

# Calculate p-values for alternative tests
for(fx in names(df[-1])) {
  y <- df[,fx]
  r <- log(y[2:length(y)]) - log(y[1:(length(y)-1)])
  r <- r[(length(r) - M + 1):length(r)]
  pvalue_AVR <- AutoBoot.test(r, nboot = 500, wild = "Normal")$pval
  pvalue_AQ <- Auto.Q(r, lags = 10)$Pvalue
  cat(fx, pvalue_AVR, pvalue_AQ, "\n")
}
