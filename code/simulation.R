remove(list = ls())
library(keras)
library(vrtest)
source("TS_gen.R")

# ==================================================================================================
#    Calculate power of automatic variance ratio test (Kim, 2009)
# ==================================================================================================

# inov = 1: GARCH
# inov = 2: Stochastic volatility

N = 1000
for(M in c(100,500,1000)) {  
  vr = rep(0, N)
  for(i in 1:N) {
    #y = AR(size = M, phi = 0.1, omega = 0.5, alpha = 0.1, beta = 0.75, st = FALSE, inov = 2)
    y = ARFIMA(size = M, omega = 0.5, alpha = 0.1, beta = 0.75, st = FALSE, inov = 2)
    vr[i] = AutoBoot.test(y, nboot = 500, wild = "Normal")$pval
  }
  print(100 * sum(vr < 0.05) / N)
}

# ==================================================================================================
#    Calculate power of automatic variance ratio test (Escanciano and Lobato, 2009)
# ==================================================================================================

N = 1000
for(M in c(100,500,1000)) {  
  vr = rep(0, N)
  for(i in 1:N) {
    #y = AR(size = M, phi = 0.1, omega = 0.5, alpha = 0.1, beta = 0.75, st = FALSE, inov = 2)
    y = ARFIMA(size = M, omega = 0.5, alpha = 0.1, beta = 0.75, st = FALSE, inov = 2)
    vr[i] = Auto.Q(y, lags = 10)$Pvalue
  }
  print(100 * sum(vr < 0.05) / N)
}

# ==================================================================================================
#    Calculate power of the convolutional network
# ==================================================================================================

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

N = 100000
power = NULL
for(M in c(100,500,1000)) {
  H0 = matrix(0, nrow = N, ncol = M)
  H1 = matrix(0, nrow = N, ncol = M)
  for(i in 1:N) {
    H0[i,] = AR(size = M, phi = 0, omega = 0.5, alpha = 0.1, beta = 0.75, st = FALSE, inov = 2)
    #H1[i,] = AR(size = M, phi = 0.1, omega = 0.5, alpha = 0.1, beta = 0.75, st = FALSE, inov = 2)
    H1[i,] = ARFIMA(size = M, omega = 0.5, alpha = 0.1, beta = 0.75, st = FALSE, inov = 2)
  }
  xT = rbind(H0[1:(N/2),], H1[1:(N/2),])
  xV = rbind(H0[(N/2+1):N,], H1[(N/2+1):N,])
  yT = c(rep(0, N/2), rep(1, N/2))
  yV = c(rep(0, N/2), rep(1, N/2))
  xT  = array_reshape(xT, c(dim(xT), 1))
  xV  = array_reshape(xV, c(dim(xV), 1))
  m = convnet(M)
  history = m %>% fit(xT, yT, batch_size = 32, epochs = 50, validation_data = list(xV, yV))
  scores  = m %>% predict(xV)
  scores_H0 = scores[1:(N/2)]
  scores_H1 = scores[(N/2+1):N]
  crit_value = quantile(scores_H0, probs = 0.95)
  print(100*sum((scores_H1 > crit_value) / length(scores_H1)))
  power = c(power, 100*sum((scores_H1 > crit_value) / length(scores_H1)))
}
power

