library(lavaan)
library(psych)
library(rjags)
library(boot)
library(rsem)

data <- read.table('robustsem.ex2.txt')
colnames(data) <- c('V1','V2','V3','V4','income','grade')

data1 <- data[-1, 1:4]
data1[] <- lapply(data1, as.numeric)


#missing rate
proportion <- round(100*apply(data[-1,1:4], 2, function(x) mean(is.na(x))),3)

# descriptive statistics
s <- describe(data1)


data <- data1[,4]
skewness <- function(data, i) {
  d <- data[i]
  return(describe(d)$skew)
}

# bootstrap
bootstrap_result <- boot(data, skewness, R=1000)
bootstrap_se <- sd(bootstrap_result$t)  


# skewness and z
skew_value <- describe(data)$skew
z_value <- skew_value / bootstrap_se

# symmetrical distribution test
p_value <- 2 * (1 - pnorm(abs(z_value)))

print(paste("Z-value:", z_value))
round(z_value,2)
print(paste("P-value:", p_value))


#RMB GCM
Niter = 150000
burnIn = 75000

tau<-0.5
y <- data1
N <- nrow(y)
Time <- ncol(y)
Lambda <- matrix(c(rep(1,Time), (c(1:Time)-1)), nrow = Time)


#------------------------------------------
# GCM
dat <- list("N" = N, "y" = y, "tau" = tau, "Time" = ncol(y))

# Set initial values
initial=list(".RNG.name" = "base::Wichmann-Hill", ".RNG.seed" = 20231210)

jags.gmm.md3 <- jags.model( file = "gcm_median.txt", 
                            data=dat, inits= initial, n.chains=1, n.adapt=1000 )

params <- c("par")
samps.gmm.md3 <- coda.samples(jags.gmm.md3, params, n.iter = Niter)
smp.md3 <- window(samps.gmm.md3[[1]], start = burnIn)

#Geweke convergence test
geweke.md3 <- apply(smp.md3, 2, function(x) geweke.diag(x)$z)
geweke.md3

# parameters estimation for the RMB GCM
summary(smp.md3)$quantiles[1:5,3]

#traceplot
race_plot <- traceplot(samps.gmm.md3) 

# credible interval
# Convert MCMC samples to matrix format
samples_matrix <- as.matrix(samps.gmm.md3)[, 1:5]  # Extract first 5 parameters

# Compute 95% credible intervals for each parameter
credible_intervals <- apply(samples_matrix, 2, quantile, probs = c(0.025, 0.975))

# Print the results
print(credible_intervals)


#RMB GCM with the selection model
y <- data1
N <- nrow(y)
Time <- ncol(y)
Lambda <- matrix(c(rep(1,Time), (c(1:Time)-1)), nrow = Time)

#------------------------------------------
# GCM
dat <- list("N" = N, "y" = y, "tau" = tau, "Time" = ncol(y))

# Set initial values
initial=list(".RNG.name" = "base::Wichmann-Hill", ".RNG.seed" = 20231210)

jags.gmm.md4 <- jags.model( file = "gcm_median_MNAR.txt", 
                            data=dat, inits= initial, n.chains=1, n.adapt=1000 )
params <- c("par")
samps.gmm.md4 <- coda.samples(jags.gmm.md4, params, n.iter = Niter)

smp.md4 <- window(samps.gmm.md4[[1]], start = burnIn)
geweke.md4 <- apply(smp.md4, 2, function(x) geweke.diag(x)$z)

summary(smp.md4)
summary(smp.md4)$quantiles[1:5,3]

#Geweke convergence test
geweke.md3 <- apply(smp.md3, 2, function(x) geweke.diag(x)$z)
geweke.md3

# parameters estimation for the RMB GCM
summary(smp.md3)$quantiles[1:5,3]

#traceplot
race_plot <- traceplot(samps.gmm.md3) 

# credible interval
# Convert MCMC samples to matrix format
samples_matrix <- as.matrix(samps.gmm.md3)[, 1:5]  # Extract first 5 parameters

# Compute 95% credible intervals for each parameter
credible_intervals <- apply(samples_matrix, 2, quantile, probs = c(0.025, 0.975))

# Print the results
print(credible_intervals)


#FIML
gcmodel<-'i =~ 1*V1 + 1*V2 + 1*V3  + 1*V4
	s =~ 0*V1 + 1*V2 + 2*V3 + 3*V4'

res.FIML <- growth(gcmodel, data1, missing = "fiml")
parameterEstimates(res.FIML)

#TSRE
pat <- rsem.pattern(data1)
musig <- rsem.emmusig(pat, varphi=.1, max.it = 100000)
res2 <- growth(gcmodel, sample.cov=musig$sigma, sample.mean=musig$mu, sample.nobs=399,mimic='EQS')
ascov <- rsem.Ascov(pat, musig)
robust.se <- rsem.se(res2, ascov$Gamma)
res.est <- cbind(coef(res2), robust.se$se[[1]])[5:9,]
round(res.est,3)







