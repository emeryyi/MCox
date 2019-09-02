#parameters
p <- 7
K <- 4
n <- 30
#raw data
X = round(matrix(rnorm(n*p,0,1),n,p),1)
beta <- seq(p) -(p+1)/2
#eventually do weibull distribution or something else
y <- round(rexp(n,exp(X %*% beta)),1)
data <- data.frame(
    time = y,
    censored = rbinom(n,1,0.3),
    task = paste("task", sample.int(K,n, replace=TRUE),sep=""),
    X = X
)
time_index <- 1
censored_index <- 2
task_index <- 3

obj = Preprocessing(data, task_index, time_index, censored_index)[[1]]

# llk = LogLikelihood(obj)
# llk
# 
derivatives = PartialDerivatives(obj, 1)
derivatives$gradient
derivatives$hessian

out = GPG_Cycle(obj, 1e-1)


out$beta
out$gradient / out$hessian


out = GPG_Cycle(obj, 1e-5)
out
