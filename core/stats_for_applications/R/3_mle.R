set.seed(1614)

# population parameters / sim data - X ~ Exponential(5)
lambda <- 5

n <- 1000
sample <- rexp(n, lambda)

# log likelihood / joint density function
# negative since optim is a minimizer
nll.exp <- function(lambda, x) {
    return(-1 * (log(lambda) - lambda * mean(x)))
}

# optimize the LL function
o <- optim(
    par = c(0.5), #initial param guess for lambda
    fn = nll.exp,
    x = sample,
    method = "BFGS"
)

(lambda_hat <- o$par)
#> 4.99

# which should equal 1 over the mean:
1 / mean(sample)
#> 4.99

# visually:
ll.exp <- \(p) -1 * nll.exp(p, x=sample)
params <- seq(1, 10, by = 0.5)
plot(params, sapply(params, ll.exp), ylab = "LogL(param | X)")
abline(v = lambda_hat, col = "red")


# confidence intervals:
## Recall that fisher information is equivalent to:
## -Expectation[second derivative of likelihood]
## The second derivative of the exponential's likelihood being:
##     -n / lambda^2
## Neither are random variables, so their expectation is themselves, leaving:
##     I = -(-n / lambda^2) = n / lambda^2
## Now, recall that the our estimator, "lambda_hat", converges in distribution to
## N(0, 1/I(theta*)), which means we know the distribution of our
## estimator is normal, and we know the variance, allowing us to calculate our
## classic confidence interval with the z-score and std. dev. (using
## I(lambda_hat) in place of I(theta*)), which simplifies for a 95% CI to:
##    lambda_hat +- 1.96 * (lambdd_hat / sqrt(n))
## Recalling that lambda_hat = 1 / X_bar (where X_bar is the sample mean),
##    1/X_bar +- 1.96 / (X_bar * sqrt(n))
sample_means <- replicate(1000, mean(rexp(n, lambda)))

# calculate CIs
data <- data.frame(
    fit = 1 / sample_means,
    lower = 1 / sample_means - 1.96 / (sample_means * sqrt(n)),
    upper = 1 / sample_means + 1.96 / (sample_means * sqrt(n))
)

# calculate proportion of CIs that include the true value
mean(data$lower < lambda & data$upper > lambda)
#> 0.95

# plot
plot(1:6, head(data$fit),
     ylim = range(c(data$lower, data$upper)),
     xlab = "Estimate Number", ylab = "Estimate and CI")
arrows(1:6, head(data$lower), 1:6, head(data$upper),
       length = 0.05, angle = 90, code = 3)
abline(a = lambda, b = 0, col = "red")


# variance of MLE:
## Recall that Var(theta_hat) = 1 / I(theta_hat)
mles <- replicate(1000, 1 / mean(rexp(n, lambda)))

# true variance
I <- n / (lambda^2) ## calculate Fisher Information
1 / I
# 0.025

# empirical variance
var(mles)
# 0.02337
