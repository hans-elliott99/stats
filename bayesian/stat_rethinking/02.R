#' Statistical Rethinking Chapter 2
#' 
#' What proportion of the globe is water?
#' 
#' In this chapter the toy example is a globe which we "toss" and randomly
#' select a point on, recording if the point is on Water or Land.
#' Through this sampling procedure, we can use Bayesian analysis to estimate the
#' true proportion of water on Earth.
#' 
#' We can think of this procedure as recording k successes (i.e., the number of
#' waters) out of n trials (i.e., tosses). This corresponds to a Binomially
#' distributed random variable with parameters n, k, and p, where n and k are
#' known to us due to the sampling procedure, but p is unknown.
#' p is the probability of observing a success (water). It is equivalent to the
#' proportion of water, so the true p is ~70%.
#' 
#' 
#' Stat. Rethinking presents the steps to Bayesian analysis as:
#' - Determine a story for the data and model the data generating process using
#'   DAG.
#'    - In this case, we can draw a simple DAG that shows directed edges pointing
#'      from n and p to "W" and "L", where W is the occurrence of a water and L
#'      the occurrence of a land. In other words, the proportion of water, p,
#'      affects the number of Ws and Ls that we see, as does the sample size
#'      automatically (for consistency with Binomial we might show directed
#'      edges from n and p point to k)
#' - Next we define a prior distribution for the parameter that needs to be
#'   estimated (potentially an uninformative prior).
#' - The final step is to use Bayes' rule to compute the posterior probability -
#'   the probability distribution over p conditional on the data we've observed.
#'   This is proportional to the likelihood times the prior (up to a normalizing
#'   constant)
#'    - In more complex cases, the posterior will be difficult or impossible to
#'      estimate analytically and so numerical solutions are needed to
#'      approximate it

#
# Grid Approximation of Posterior
#
# (simply evaluating the posterior distribution over a grid of possible values)
posterior_9C6 <- function(prior, ps) {
  # 6 waters in 9 globe tosses
  likelihood <- dbinom(6, size = 9, prob = ps)
  posterior <- likelihood * prior
  return(list(p = ps,
              prior = prior,
              posterior = posterior / sum(posterior)
              ))
}

len <- 100
p_grid <- seq(0, 1, length = len)

# uniform prior
prior <- rep(1, len)
mod <- posterior_9C6(rep(1, len), p_grid)
par(mfrow = c(1, 2))
plot(mod$p, mod$prior, type = "l", main = "prior")
plot(mod$p, mod$posterior, type = "l", main = "posterior")

# gate prior
prior <- ifelse(p_grid < 0.5, 0, 1)
mod <- posterior_9C6(prior, p_grid)
par(mfrow = c(1, 2))
plot(mod$p, mod$prior, type = "l", main = "prior")
plot(mod$p, mod$posterior, type = "l", main = "posterior")

# pointy prior
prior <- exp(-5 * abs(p_grid - 0.5))
mod <- posterior_9C6(prior, p_grid)
par(mfrow = c(1, 2))
plot(mod$p, mod$prior, type = "l", main = "prior")
plot(mod$p, mod$posterior, type = "l", main = "posterior")
par(mfrow = c(1,1))


#
# Quadratic Approximation & Maximum a Posteriori
#
# Grid approximation is too cumbersome in most cases.
# Quadratic/Gaussian approximation: posterior can be usefully approximated by a
# Gaussian.
# The approximation is usually equivalent to a Gaussian with mean equivalent to
# the MLE and std. dev. equivalent to the SE of the MLE. 
# 
# MAP: arg max of posterior distribution wrt p, which is proportional to:
#     arg max wrt p of (likelihood * prior)
# Likelihood (binomial):
#     choose(n, k) * p^k * (1 - p)^(n - k)
lik.binom <- function(n, k, p) {
  return( choose(n, k) * p^k * (1 - p)^(n - k) )
}
# If the prior is uninformative, this is equivalent to arg maxing the likelihood
# (i.e., MLE).
# Here I use the uniform prior Beta(1,1)
# (note, dbeta(ppoints(100), 1,1) == rep_len(1,100))


# data: sampled 6 waters in 9 globe tosses - estimate p (prob of water)
n <- 9
k <- 6
par(mfrow = c(1,3))

plot(ppoints(100), rep_len(1, 100), type = "l",
     main = "Beta(1,1)", xlab = "p", ylab = "Density")

plot(ppoints(100), lik.binom(n, k, ppoints(100)), type = "l",
     main = "likelihood", sub = "MLE, MLE +/- 2 Std Errors",
     xlab = "p", ylab = "Likelihood")

# maximum likelihood estimation
o <- optim(
  par = c(0.5),                    # initial paramater guess
  fn = lik.binom,                  # function to optimize
  n = n,                           # preset fn argument
  k = k,                           # preset fn argument
  control = list(fnscale = -1),    # scale fn by -1 since optim minimizes
  method = "Brent",                # optimization method
  lower = .Machine$double.eps,     # lower bound for par in Brent optim
  upper = 1 - .Machine$double.eps  # upper bound for par in Brent optim
)
p_mle <- o$par
abline(v = p_mle, col = "red")

# SE of the MLE:
# The SE of the MLE can be computed analytically using the Fisher information,
# https://math.stackexchange.com/questions/396982/fisher-information-of-a-binomial-distribution
# The variance of the MLE can be estimated by 1 / Information(p_hat).
# I(p) = n / p(1-p), Var(hat_p_mle) = p(1-p) / n
fisher.binom <- function(n, p) {
  return( n / (p * (1 - p)) )
}
var_mle <- 1 / fisher.binom(n = n, p = p_mle)
abline(v = p_mle - 2 * var_mle, col = "pink")
abline(v = p_mle + 2 * var_mle, col = "pink")


# Gaussian approximation
lik <- lik.binom(n, k, ppoints(100))
plot(ppoints(100), lik / sum(lik), type = "l",
     main = "posterior & quadratic approx",
     xlab = "p", ylab = "Density")
qa <- dnorm(ppoints(100), mean = p_mle, sd = sqrt(var_mle))
lines(ppoints(100), qa / sum(qa), col = "red")


# Repeat for different sample sizes...
#
qa_plot <- function(n, k) {
  p <- ppoints(100)
  lik <- lik.binom(n, k, p)
  # maximum likelihood estimation
  o <- optim(
    par = c(0.5),                    # initial paramater guess
    fn = lik.binom,                  # function to optimize
    n = n,                           # preset fn argument
    k = k,                           # preset fn argument
    control = list(fnscale = -1),    # scale fn by -1 since optim minimizes
    method = "Brent",                # optimization method
    lower = .Machine$double.eps,     # lower bound for par in Brent optim
    upper = 1 - .Machine$double.eps  # upper bound for par in Brent optim
  )
  p_mle <- o$par
  var_mle <- 1 / fisher.binom(n = n, p = p_mle)
  qa <- dnorm(p, mean = p_mle, sd = sqrt(var_mle)) # quadratic approx

  # plots
  par(mfrow = c(1,3))
  ## prior
  plot(p, rep_len(1, 100), type = "l",
       main = "Beta(1,1)", xlab = "p", ylab = "Density")
  ## likelihood
  plot(p, lik, type = "l",
       main = paste0("likelihood, n=",n, ", k=",k),
       sub = "MLE, MLE +/- 2 Std Errors",
       xlab = "p", ylab = "Likelihood")
  abline(v = p_mle, col = "red")
  abline(v = p_mle - 2 * var_mle, col = "pink")
  abline(v = p_mle + 2 * var_mle, col = "pink")
  # posterior + gaussian approximation
  plot(p, lik / sum(lik), type = "l",
       main = "posterior & quadratic approx",
       xlab = "p", ylab = "Density")
  lines(p, qa / sum(qa), col = "red", lty = 2)
  par(mfrow=c(1,1)) # reset
}


# approximation improves as sample size n increases
# (k could be held fixed, but we assume same fraction of water 0.7)
qa_plot(n =  9, k = floor(9 * 0.7))
qa_plot(n = 18, k = floor(18 * 0.7))
qa_plot(n = 36, k = floor(36 * 0.7))
qa_plot(n = 81, k = floor(81 * 0.7))
qa_plot(n = 729,k = floor(729 * 0.7))


# The posterior distribution for the binomial with an uninformative prior of
# Beta(1,1) is known to be a Beta distribution w/ parameters 1 + k and 1 + n - k,
# (where k is the number of waters in our sample)
# Thus, we can replicate the above posterior estimation:
curve(dbeta(x, k + 1, n - k + 1), from = 0, to = 1,
      main = "Beta Posterior & Quadratic Approx.",
      xlab = "p", ylab = "Density")
curve(dnorm(x, p_mle, sqrt(var_mle)), col = "red", lty = 2, add = TRUE)




