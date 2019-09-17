
library(MCox)
#parameters
p = 20
K = 2
n = 200
lambda = 0.01 # rate parameter in h0
rho = 1.0 #shape parameter in h0
rate_censored = 0.001 # rate parameter of the exponential distribution of C
#raw data
X = round(matrix(rnorm(n*p,0,1),n,p),1)
beta = (seq(p) -(p+1)/2)*2/p
beta[min(p,10):p] = 0
#weibull distribution baseline
#ref: https://stats.stackexchange.com/questions/135124/how-to-create-a-toy-survival-time-to-event-data-with-right-censoring
failure_time = (-log(runif(n)) / (lambda * exp(X %*% beta)))^(1/rho)
censored_time = rexp(n, rate = rate_censored)
time = pmin(failure_time, censored_time)
censored = as.numeric(failure_time > censored_time)
data <- data.frame(
    time = round(time,1),
    censored = censored,
    task = paste("task", sample.int(K,n, replace=TRUE),sep=""),
    X = X
)
time_index <- 1
censored_index <- 2
task_index <- 3

out = MCox(data, task_index, time_index, censored_index)

out$iError




# 
# obj = Preprocessing(data, task_index, time_index, censored_index)[[1]]
# obj
# 
# llk = LogLikelihood(obj, beta)
# llk$logLik
# llk = LogLikelihood(obj)
# llk$logLik
# 
# # 
# derivatives = PartialDerivatives(obj, 1)
# derivatives$gradient
# derivatives$hessian
# 
# out = GPG_Cycle_Backtracking(obj, 1e-1)
# out = GPG_Cycle(obj, 1e-1)
# out$beta
# out$gradient / out$hessian
# out = GPG_Descent(obj, 1e-1, 1, 0.8)
# abs(out$beta)>0.00001
# out$nCycles
# out$nUpdates
# out$beta
# plot(out$beta)
