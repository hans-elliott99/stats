set.seed(1614)


# ===== NORMAL ==========
# population/true params - X ~ N(5, 2^2)
mu <- 5
sigma <- 2

# sample
n <- 1000
samples <- rnorm(n, mu, sigma)

# calculate estimate - 2 parameters, 2 sample moments
# first moment = E[X] ~= Mean(X^1)
# second moment = E[X^2] ~= Mean(X^2)
# mu = E(X) = first moment
# sigma^2 = E(X^2) - E(X)^2 = second moment - (first moment)^2
sample_means <- mean(samples)
sample_vars  <- (sum(samples^2) / n) - sum(samples / n)^2

sample_means
# 5.01
sample_vars
# 4.05


# ===== GAMMA ==========
# population parameters - X ~ Gamma(2, 2)
shape <- 2
scale <- 2

# sample
n <- 1000
sample <- rgamma(n, shape, rate)
# moment_1 = shape / rate
# moment_2 = (shape / rate^2) + (shape^2) / (rate^2)
# solving, find:
# shape = (moment_1^2) / (moment_2 - moment_1^2)
# scale = (moment_1) / (moment_2 - moment_1^2)

# sample moments
m1 <- mean(sample)
m2 <- mean(sample^2)

# parameter estimates
(shape <- m1^2 / (m2 - m1^2))
# 1.96
(scale <- m1 / (m2 - m1^2))
# 1.91
