#set.seed(1)
library(MCox)
#parameters
p = 50
K = 1
n = 10000
lambda = 0.01 # rate parameter in h0
rho = 1.0 #shape parameter in h0
rate_censored = 1e-3 # rate parameter of the exponential distribution of C
#raw data
X = round(matrix(rnorm(n*p,0,1),n,p),1)*2
beta = (seq(p) -(p+1)/2)*2/p
beta[min(p,5):p] = 0
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
par(mfrow=c(2,1), mar=c(2,2,2,2))

out = MCox(data, task_index, time_index, censored_index, maxIteration = 1000)

out$iError

matplot(t(out$betaNorm), type = "l", lty=1)
abline(h = abs(beta), col=1:6, lty=3)

out$nBeta
out$nCycles
out$nUpdates



# compare with glmnet
library(glmnet)
y=cbind(time=data$time+1,status=1-data$censored)
fit=glmnet(X,y,family="cox", lambda = out$lambda, weights = rep(1/n,n), 
           thresh = 1e-14, standardize=TRUE)

# overlay paths
matplot(abs(t(coef(fit))), add = T, type = "l",lty=2)
# max difference along the whole path
max(abs(abs(coef(fit)) - out$betaNorm))

matplot(abs(abs(t(coef(fit)))-t(out$betaNorm)), type = "l",lty=2)
# check if all entered at same point
out$nBeta-fit$df

