#
# simulations from probability class
#

# PROBABILITY I ===============================================================

scrabble_sim <- function(tnum = c("E" = 5, "A" = 4, "N" = 3, "B" = 2),
                         C = c("E", "E", "A", "N"),
                         niter = 1000) {
    # A Scrabble bag contains various numbers of tiles of letters E,A,N,B
    #  (specified by `tnum`)
    #  (Original problem is 5E, 4A, 3N, 2B)
    # Let C be the event of drawing a specific group of tiles, for example, let
    #  C be the event of drawing 2 Es, 1 A, and 1 N.
    # What is P(C)?
    tiles <- c(
        rep("E", tnum["E"]),
        rep("A", tnum["A"]),
        rep("N", tnum["N"]),
        rep("B", tnum["B"])
    )
    true_p <- unname( (choose(tnum["E"], sum(C=="E")) *
                       choose(tnum["A"], sum(C=="A")) *
                       choose(tnum["N"], sum(C=="N")) *
                       choose(tnum["B"], sum(C=="B"))
                       )
                      / choose(sum(tnum), length(C)) )
    
    nC <- length(C)
    sC <- sort(C)
    res <- 0
    for (i in 1:niter) {
        x <- sample(tiles, size = nC, replace = FALSE)
        res <- res + all(sort(x) == sC) / niter
    }
    message("True probability: ", true_p)
    message("Estimated probability: ", res)
}

if (FALSE) {
scrabble_sim(niter = 10000)
scrabble_sim(C = c("E", "E", "A", "N", "B"), niter = 10000) # same true prob. as above
scrabble_sim(C = c("E", "E", "E", "E", "E"))
}


stick_break_sim <- function(niter = 1000, stick_len = 1000) {
    # We break a stick at a uniformly chosen random location.
    # Find the probability that the shorter piece is less than 1/5th of the original.
    stick <- seq(0, 1, length.out = stick_len)
    res <- 0
    for (i in 1:niter) {
      brk <- runif(1)
      lwr_len <- length(stick[stick < brk])
      upr_len <- length(stick[stick > brk])
      res <- res + ( lwr_len < stick_len/5 || upr_len < stick_len/5 ) / niter
    }
    message("True probability: ", 2/5)
    message("Estimated probability: ", res)
}
if (FALSE) {
  stick_break_sim(niter = 10000, stick_len = 1000)
}


uniform_dartboard_sim <- function(niter = 1000,
                                  board_r = 9,
                                  bullseye_r = 1/4) {
    # A disk-shaped dart board has radius 9, and its bullseye has radius 1/4.
    # Find the probability that a dart lands in the bullseye, if it is thrown
    # uniformly at random at the board.
    #
    # Note - uses "rejection sampling" so only about pi/4 of niter samples are
    # used!
    # (Because pi/4 is the ratio of the board to the square that contains it:
    #  let D = diameter of circle, then Square Area = D^2, Circle Area = pi*(D/2)^2,
    #  so Circle Area / Square Area = pi(D/2)^2 / D^2 = pi/4)
    stopifnot(bullseye_r < board_r)
    true_p <- (pi * (bullseye_r)^2) / (pi * board_r^2) # area of bullseye / area of board
    
    res <- 0
    nthrow <- 0
    nskip <- 0
    for (i in 1:niter) {
      x <- runif(1, min = -board_r, max = board_r)
      y <- runif(1, min = -board_r, max = board_r)
      if (x^2 + y^2 > board_r^2) {
        nskip <- nskip + 1
        next # not on board, ignore sample
      } 
      hit <- (x^2 + y^2) <= bullseye_r^2
      res <- res + as.numeric(hit) 
      nthrow <- nthrow + 1
    }
    message("True probability:      ", true_p)
    message("Estimated probability: ", res/nthrow)
    message("(Proportion of used samples: ", (niter-nskip)/niter, ")")
    return(c(true_p = true_p, est_p = res/nthrow))
}
if (FALSE) {
x <- uniform_dartboard_sim(niter = 10000, board_r = 9, bullseye_r = 1/4)
x <- uniform_dartboard_sim(niter = 10000, board_r = 9, bullseye_r = 1/4)

boardr <- 9
rads <- seq(0, boardr-1, by = 0.05)
ps <- vapply(rads,
             \(r) suppressMessages(uniform_dartboard_sim(1000, boardr, r)),
             numeric(2))

plot(rads, ps["true_p", ],
     type = "l",
     main = paste0("Probability of hitting bullseye on uniform dartboard (r = ",br, ")"),
     xlab = "Bullseye radius", ylab = "Probability of bullseye")
lines(rads, ps["est_p", ], col = "red", lty = 2)
legend("topright", c("True", "Estimated"),
       col = c("black", "red"), lty = c(1, 2))
}



die_until_k_sim <- function(win_num = 4, niter = 1000, report_k = c(1, 5, 10)) {
    # Roll a fair die repeatedly until we see the number 4 appear, then stop and
    # record the number of rolls it took ('k', the experiment's outcome).
    # Note: changing the 'winning number' from 4 won't change the probabilities
    # (still a 5/6 chance of losing and 1/6 chance of winning on each roll)
    # unless you add additional winning numbers!
    # Use report_k to change which 'k's have their probabilities reported (i.e.,
    # for which experiment outcomes do you want to know the probabilities).
    # Note:
    # This is actually the geometric distribution, which has PMF (1-p)^(k-1) * p,
    # where p is the probability of success and thus the number of winning numbers
    # out of 6.
    # (Technically this parameterization is the zero-truncated geometric
    #  distribution, and thus differs from R's which starts at k=0 instead of
    #  k=1. See `dgeom(0:1, p=1/6)`)
    stopifnot( all(win_num <= 6 & win_num > 0) )
    outcomes <- integer(100)
    for (i in 1:niter) {
        rolls <- sample(1:6, size = length(outcomes), replace = TRUE)
        first_k <- min(which(rolls %in% win_num))
        if (is.infinite(first_k)) next # drop sample if 4 not rolled in first 100
        outcomes[first_k] <- outcomes[first_k] + 1
    }
    p <- length(win_num) / 6 # probability of winning in each trial
    trues <- round(vapply(report_k, \(k) (1-p)^(k-1) * p, numeric(1)),
                   digits = 2)
    ests <- round(outcomes/sum(outcomes), digits = 3)[report_k]
    message("True probabilities:")
    print(setNames(trues, paste0("k = ", report_k)))
    message("Estimated probabilities:")
    print(setNames(ests, paste0("k = ", report_k)))
    return(outcomes)
}
if (FALSE) {
die_until_k_sim(win_num = 4, niter = 1000, report_k = c(1, 5, 10, 20))
dgeom(c(0,4,9,19), p=1/6) ## R's geometric dist.
# increase prob. of success in each trial by adding more winning numbers
die_until_k_sim(win_num = c(4,5), niter = 1000, report_k = c(1, 5, 10, 20))
die_until_k_sim(win_num = c(4,5,6), niter = 1000, report_k = c(1, 5, 10, 20))
}



alice_bob_conintoss <- function() {
    niter <- 10000
    ntoss_bob <- 10
    ntoss_alice <- ntoss_bob + 1
    o <- matrix(0, nrow = niter, ncol = 3) # losing, tied, winning
    y <- logical(niter)
    t <- matrix(0, nrow = niter, ncol = 2)
    for (i in seq(niter)) {
        ntails_bob <- sum(sample(c(0, 1), size = ntoss_bob, replace = TRUE))
        ntails_alice <- sum(sample(c(0, 1), size = ntoss_bob, replace = TRUE))
        if (ntails_alice < ntails_bob) {
            o[i, 1] <- 1 # already losing
        } else if (ntails_alice == ntails_bob) {
            o[i, 2] <- 1 # tied - can win or lose based on next toss
        } else {
            o[i, 3] <- 1 # already winning
        }
        ntails_alice <- ntails_alice + sample(c(0, 1), size = 1, replace = TRUE)
        y[i] <- ntails_alice > ntails_bob
        t[i, 1] <- ntails_bob
        t[i, 2] <- ntails_alice
    }
}




# probability of drawing Q, then K, then A from a standard deck
draw_QKA <- function() {
    cards <- rep(c(2:10, "J", "Q", "K", "A"), 4)
    
    res <- vapply(1:1e6, \(i) {
        all(sample(cards, size = 3, replace = FALSE) == c("Q", "K", "A"))
    }, logical(1))
    list(
        true = 8/16575,
        est = mean(res)
    )
}


res <- vapply(1:1e5, \(i) {
    n_rnds <- 4
    AnnBill <- data.frame(
        matrix(sample(c("R", "P", "S"), size = n_rnds * 2, replace = TRUE),
               ncol = 2, dimnames = list(NULL, c("Ann", "Bill")))
    )
    AnnBill <- within(AnnBill, {
        win <- ifelse(
            (Ann == "R" & Bill == "S") |
            (Ann == "P" & Bill == "R") |
            (Ann == "S" & Bill == "P"),
            TRUE, FALSE)
        ssum <- cumsum(win)
    })
    # if first win is on `n_rnds` round
    if (AnnBill[n_rnds, "ssum"] == 1 && AnnBill[n_rnds - 1, "ssum"] != 1)
        TRUE else FALSE
}, logical(1))
mean(res)
8/81



# PROBABILITY II ==============================================================


# draw 5 distinct numbers from set {1,2,...,40}
# what is the probability dist. over the sum of the 5 numbers?
sum_of_5_distinct <- function() {
    niter <- 1e9
    res <- vapply(1:niter, \(i) {
        sum(sample(1:40, size = 5, replace = FALSE))
    }, numeric(1))
    return(res) 
}

ex <- sum_of_5_distinct()
hist(ex,
     breaks = 100, col = "lightblue",
     main = "Sum of 5 distinct numbers from {1,2,...,40}")




# X ~ Geom(p). We know X is even. Y = X/2.
# What is the probability distribution of Y?
# i.e., what is P(Y = k | X even) ?
# Y ~ Geom(1 - q^2), where q = 1 - p
even_geom <- function(success_prob, niter = 1e6) {
    res <- vapply(1:niter, \(i) {
        x <- NA
        while (TRUE) {
            x <- rgeom(1, prob = success_prob)
            if (x %% 2 == 0) break
        }
        x
    }, numeric(1))
    return(res)
}

x <- even_geom(success_prob = 0.5, niter = 1e5)
hist(x/2,
     breaks = 100, col = "lightblue",
     main = "Y = X/2, X ~ Geom(0.5)",
     probability = TRUE)

x_true <- rgeom(1e5, prob = 1 - 0.5^2)
hist(x_true/2,
     breaks = 100, col = "lightblue",
     main = "X ~ Geom(1 - 0.5^2)",
     probability = TRUE)



# You flip a coin n times - but don't assume the coin is fair or that the flips
# are independent.
# Let X_n be the r.v. that counts the number of heads in the first n flips.
# Is the function that takes n -> P(X_n <= 5) always increasing or always decreasing?
dependent_coin_toss <- function(nsim, ntoss) {
    flipper <- function(last_flip) {
        if (last_flip == 1) {
            # if last flips is heads (1), we flip a fair coin
            sample(c(0,1), 1, prob = c(0.5, 0.5))
        } else {
            # if last flip is tails (0), more likely to get tails again
            sample(c(0,1), 1, prob = c(0.8, 0.2))
        }
    }
    experiment <- function(n) {
        vec = numeric(n)
        last = 1
        for (i in seq(n)) {
            last = flipper(last)
            vec[i] = last
        }
        return(vec)
    }
    mat <- matrix(NA, nrow = nsim, ncol = ntoss)
    for (sim in 1:nsim) {
        for (toss in 1:ntoss) {
            mat[sim, toss] <- sum(experiment(toss)) # X_n
        }
    }
    return(mat)
}
S <- dependent_coin_toss(1000, 30) # 10000 for much smoother curve (takes longer)
probs <- apply(S <= 5, 2, mean)
plot(probs, type = "l", col = "orange", lwd = 1,
     xlab = "n", ylab = "P(X_n <= 5)",
     main = "P(X_n <= 5) vs n")
head(probs, n = 10) # for n <= 5, P(X_n <= 5) = 1
# cool note:
# For the "biased" coin condition, if we make it more likely to flip tails after
# a tails, then P(Xn <= 5) declines relatively slowly in n compared to if we make
# it more likely to flip heads after flipping tails. In the first case, flipping T
# means we are likely to flip T again, so the number of heads is likely to be
# lower. In the second, flipping T makes it more likely you'll switch back to H,
# so the number of heads is likely to be higher, and P(Num heads <= 5) decreases
# quicker!
# E.g., set prob = c(0.2, 0.8) in flipper() for when last_flip == 0 and see a
# very steep decline in P(Xn <= 5) as n increases.



# X1, X2 are independent r.v.
# P(X1 = 1) = P(X1 = -1) = 1/2
# P(X2 = 1) = p, P(X2 = -1) = 1-p, for some 0 < p < 1
# Y = X1 * X2, are Y, X2 independent?
# Yes
y_x2_independence <- function() {
    niter <- 1e6
    out <- data.frame(
        x1 = sample(c(1, -1), niter, prob = c(0.5, 0.5), replace = TRUE),
        x2 = sample(c(1, -1), niter, prob = c(0.7, 0.3), replace = TRUE)
    )
    out$y <- out$x1 * out$x2
    return(out)
}
ex <- y_x2_independence()
cor(ex$x2, ex$y) # ~0
# jitter values
plot(jitter(ex$x2[1:1000]), jitter(ex$y[1:1000]), col = "blue",
     xlab = "X2", ylab = "Y = X1 * X2",
     main = "X1, X2 independent?")



# X1, X2, X3 are independent r.v., each Xi ~ Exp(lambda_i).
# Y = min(X1, X2, X3). What is the distribution of Y?
# Y ~ Exp(lambda_1 + lambda_2 + lambda_3)
indep_exp_min <- function(niter = 1e6) {
    out <- data.frame(
        x1 = rexp(niter, rate = 1),
        x2 = rexp(niter, rate = 2),
        x3 = rexp(niter, rate = 3)
    )
    out$y <- pmin(out$x1, out$x2, out$x3)
    return(out)
}
ex <- indep_exp_min()
true <- rexp(1e6, rate = 1 + 2 + 3)
mean(ex$y)
mean(true)
hist(ex$y, breaks = 100, col = "lightblue",
     main = "Y = min(X1, X2, X3), X_i ~ Exp(lambda_i)")
hist(true, breaks = 100, col = "lightblue",
     main = "Y ~ Exp(lambda_1 + lambda_2 + lambda_3)")


# X1, X2, X3 are iid Exp(lambda_i).
# I = argmin(X1, X2, X3). What is the distribution of I?
# P(I = 1) = P(X1 < X2, X1 < X3) = lambda_1 / (lambda_1 + lambda_2 + lambda_3)
indep_exp_argmin <- function(niter = 1e6,
                             rates = c(1, 2, 3)) {
    out <- matrix(NA, nrow = niter, ncol = length(rates))
    for (i in seq_along(rates)) {
        out[, i] <- rexp(niter, rate = rates[i])
    }
    out <- as.data.frame(out)
    out$I <- apply(out, 1, which.min)
    return(out)
}
ex <- indep_exp_argmin()
table(ex$I) / 1e6
c(1/6, 2/6, 3/6)


# X, Y are independent r.v.
# f_X(x) = 2x for 0 <= x <= 1, 0 otherwise
# Y ~ Unif(1, 2)
# Calculate P(Y - X >= 3/2). Calculated answer: 1/24
x_inv_cdf <- function(u) {
    ## F_x(x) = \int f_X(x) dx = x^2 for 0 <= x <= 1 [CDF]
    ## F_{x}^{-1}(u) = sqrt(u) for 0 <= u <= 1       [Inverse CDF]
    stopifnot(u >= 0 & u <= 1)
    return(sqrt(u))
}
X <- x_inv_cdf(runif(1e6))
Y <- runif(1e6, min = 1, max = 2)
mean(Y - X >= 3/2)
1/24
# We can sample a random point on the square given by (0,1),(1,0),(1,2) by
# sampling an x coord and Y coord from the PDFs of X and Y respectively.
# The upper-left triangle is the region where y - x >= 3/2.
# If X was also uniform, then the probability of being in the triangle would
# be 1/8 (from the area of the triangle / area of the square).
# But since X has PDF f(x) = 2x, there is more probability density on the right
# half of the square, and so a random point on the square is more likely to fall
# outside of the triangle (which is why the answer, 1/24, is less than 1/8).
x <- seq(0, 1, by = 0.01)
y <- seq(1, 2, by = 0.01)
z <- outer(x, y, function(x, y) y - x >= 3/2)
filled.contour(x, y, z,
               color.palette = function(n) c("black", "cyan"),
               plot.title = title(main = "Y - X >= 3/2"),
               xlab = "X", ylab = "Y",
               levels = c(0, 0.5, 1))





# Sample n numbers from {1, 2, ..., n} uniformly at random with replacement.
# Picks i, j (i < j) are a match if a_i = a_j.
# What is the expected total # of matches?
#
# X = # of matches = sum_{i=1}^{n-1} sum_{j=i+1}^{n} I(a_i = a_j)
# E[X] = sum_{i=1}^{n-1} sum_{j=i+1}^{n} E[I(a_i = a_j)] (linearity of E)
#      = sum_{i=1}^{n-1} sum_{j=i+1}^{n} P(a_i = a_j)
# P(a_i = a_j) = P(a_i = k, a_j = k) for k = 1, 2, ..., n
#              = sum_{k=1}^{n} P(a_i = k, a_j = k)
#              = sum_{k=1}^{n} P(a_i = k)P(a_j = k) (independence)
#              = n * (1/n)^2 = 1/n
# (n choose 2) pairs of i, j, so E[X] = (n choose 2) * 1/n = (n-1)/2

count_matches <- function(vec) {
    mtch <- 0
    for (i in seq(1, length(vec) - 1)) {
        
        for (j in seq(i + 1, length(vec))) {
            if (vec[i] == vec[j]) {
                mtch <- mtch + 1
            }
        }
    }
    return(mtch)
}

num_matches <- function(n, niter = 1e5) {
    mat <- matrix(sample.int(n, niter * n, replace = TRUE),
                  nrow = niter, ncol = n)
    matches <- apply(mat, 1, count_matches)
    return(matches)
}

ns <- seq(2, 10, by = 1)
vec <- numeric(length(ns))
for (i in seq_along(ns)) {
    vec[i] <- mean(num_matches(ns[i]))
}
plot(ns, vec, type = "b", xlab = "n", ylab = "Expected # of matches",
     main = "Expected # of matches in n draws from {1, 2, ..., n}")
lines(ns, (ns - 1)/2, col = "red") # theoretical result


mat <- matrix(
    sample(1:2, 2 * 1e5, replace = TRUE),
    nrow = 1e5
)
# P(X1 = 1, X2 = 1) = 1/4 = P(X1 = 1) * P(X2 = 1)
mean( apply(mat, 1, \(row) sum(row[1] == 1 && row[2] == 1)) )
mean( apply(mat, 1, \(row) sum(row == 1)) / 2 )^2
# P(X1 = 2, X2 = 2) = 1/4 = P(X1 = 2) * P(X2 = 2)
mean( apply(mat, 1, \(row) sum(row[1] == 2 && row[2] == 2)) )
mean( apply(mat, 1, \(row) sum(row == 2)) / 2 )^2 

# P(X1 = k, X2 = k) = P(X1 = 1, X2 = 1) + P(X1 = 2, X2 = 2) = 1/2
mean( apply(mat, 1, \(row) sum(row[1] == 1 && row[2] == 1)) ) +
    mean( apply(mat, 1, \(row) sum(row[1] == 2 && row[2] == 2)) )






# Accident occurs at point X on road of length L, uniformly at random.
# At time of accident, ambulance is at point Y on road, also uniformly at random.
# Find the expected distance between X and Y, assuming independence.
l <- 5
X <- runif(1e6, min = 0, max = l)
Y <- runif(1e6, min = 0, max = l)
Z <- abs(X - Y)
mean(Z)
l / 3
#
sim_fn <- function(l, niter) {
    if (niter %% 2 != 0) {
        stop("niter must be even")
    }
    xy <- runif(niter * 2, min = 0, max = l)
    z <- abs(
        xy[1:niter] - xy[(niter + 1):(2 * niter)]
    )
    return(mean(z))
}
# Check that the simulation matches the theoretical result.
ls <- seq(1, 10, by = 1)
vec <- numeric(length(ls))
for (i in seq_along(ls)) {
    vec[i] <- sim_fn(ls[i], 1e6)
}
plot(ls, vec, type = "b", xlab = "l", ylab = "E[|X - Y|]",
     main = "E[|X - Y|] for X, Y ~ Unif(0, l)")
lines(ls, ls / 3, col = "red")





# I have 4 different sweaters which I randomly sample each day.
# In a 5 day week, what is the expected number of unique sweaters worn?
#
# E[X] = 4 * (1 - (3/4)^5)
niter <- 1e6
nday <- 5
mat <- matrix(
    sample.int(4, size = nday * niter, replace = TRUE),
    nrow = niter,
    ncol = 5
)
X <- apply(mat, 1, \(row) length(unique(row)))
mean(X)




#
#
num_ducks <- 10
num_hunts <- 10
p <- 0.2 # prob. of hitting a duck per hunter

misses_count <- numeric(1e4)

for (i in seq_along(misses_count)) {
    hit_ducks <- numeric(num_ducks)
    # duck chosen by each hunter
    duck <- sample.int(num_ducks, num_hunts, replace = TRUE)
    # hit success of hunters 1:num_hunts
    hit <- sample(c(0, 1), num_hunts, prob = c(1 - p, p), replace = TRUE)
    # which ducks are hit
    hit_ducks[ duck[hit == 1] ] <- 1
    
    misses_count[i] <- num_ducks - sum(hit_ducks)
}

mean(misses_count)
num_ducks * (1 - p/num_ducks)^num_hunts

# niter <- 1e3
# hit <- sample(c(0, 1), niter * num_hunts, replace = TRUE, prob = c(1 - p, p))
# hit_ducks <- sample.int(num_ducks, niter * nh, replace = TRUE)
# hit_ducks[hit == 0] <- NA
# hit_ducks <- matrix(hit_ducks, nrow = niter)




# Generate independent Unif(0,1) r.v. until I see the first Un < 1/3.
# N is the number of draws I needed to do this.
# Y is the indicator r.v. that is 1 if Un < 1/9, 0 otherwise.
# Find joint pmf of (N,Y), marginals, and determine if N,Y independent.
#
# Joint PMF is (2/3)^{k-1}(1/9) for y = 1, (2/3)^{k-1}(2/9) for y = 0, and k = 1,2,...
# N ~ Geom(1/3) so p_N(k) = (2/3)^{k-1}(1/3).
# p_Y(y) = 1/3 if y = 1 and 2/3 if y = 0
# They are independent since the joint PMF factors into the product of the marginals.
sim <- data.frame(
    N = numeric(1e4),
    Y = numeric(1e4)
)
for (iter in 1:nrow(sim)) {
    ix <- numeric(0)
    while (length(ix) == 0) {
        Ui <- runif(50) ## generate 50 since highly unlikely that none < 1/3
        ix <- which(Ui < 1/3)
    }
    N <- min(ix)
    Y <- Ui[N] < 1/9
    sim[iter, ] <- c(N, Y)
}
mean(sim[, "Y"]) # expect 1/3
mean(sim[, "N"]) # N ~ Geom(1/3), so expect 1 / (1/3) = 3

mean(sim[sim$Y == 1, "N"]) # expect 3 due to independence
mean(sim[sim$Y == 0, "N"]) # expect 3 due to independence

# test the joint pmf
test <- expand.grid(
    n = seq(1, 10),
    y = c(0, 1)
)
test$truth <- (2/3)^(test$n - 1) * ifelse(test$y == 1, 1/9, 2/9)
test$sim <- NA_real_
for (i in 1:nrow(test)) {
    n <- test[i, "n"]
    y <- test[i, "y"]
    test[i, "sim"] <- nrow(sim[sim$N == n & sim$Y == y, ])/nrow(sim)
}
plot(test$truth, test$sim,
     xlab = "Truth", ylab = "Simulated",
     main = "Joint PMF Probabilities")
abline(a = 0, b = 1)




# You have n friends. Friend birth months are iid uniformly distributed.
# Find expected number of distinct birth months among n friends.
nfriend <- 13
res <- numeric(1e6)

# bdays <- sample.int(12, nfriend * 1e6, replace = TRUE)
# for (i in seq(1, nfriend * 1e6, by = 12)) {
#     res[i/12] <- length(unique(bdays[i:(i + 11)]))
# }
for (i in seq_along(res)) {
    bdays <- sample.int(12, nfriend, replace = TRUE)
    res[i] <- length(unique(bdays))
}
# expected number of distinct birth months among n friends
mean(res)
12 * (1 - (11/12)^nfriend)

# P(X = 12)?
mean(res == 12)
(factorial(nfriend)/factorial(2)) * 1/(12^12)









# Flip a fair coin 30 times. X is the number of heads in the first 20 flips,
#  Y is the number of heads in the last 20 flips. Find correlation coef.
# Note: X and Y are not independent, because they share 10 coin flips.
#   We find positive correlation (1/2) because if there are many heads in the
#   overlapping flips, then both X and Y will be higher.
#   The non-overlapping flips are independent of each other and don't contribute
#   to the covariance or correlation.
mat <- matrix(NA_integer_, nrow = 1e4, ncol = 2)
colnames(mat) <- c("X", "Y")
for (i in 1:nrow(mat)) {
    flips <- rbinom(30, 1, 0.5)
    mat[i, "X"] <- sum(flips[1:20])
    mat[i, "Y"] <- sum(flips[11:30])
}
cov(mat[, "X"], mat[, "Y"])
cor(mat[, "X"], mat[, "Y"])






# For r.v. X, Y with given expectations, variances, and covariance, find a value
#   for a s.t. X and X + aY are uncorrelated (ie, cov(X, Y) = 0)
S <- matrix(c(1, -1,
              -1, 9),
            nrow = 2, byrow = TRUE)
xy <- MASS::mvrnorm(n = 1e6,
                    mu = c(2, 1),
                    Sigma = S)
mean(xy[, 1])
var(xy[, 1])
mean(xy[, 2])
var(xy[, 2])
cov(xy[, 1], xy[, 2])
cor(xy[, 1], xy[, 2]) # expect -1/3 since corr = cov / sqrt(var1 * var2)
cov(xy[, 1], xy[, 1] + 1 * xy[, 2]) # expect 0


## another example
S <- matrix(c(2, -1,
              -1, 9),
            nrow = 2, byrow = TRUE)
xy <- MASS::mvrnorm(n = 1e6,
                    mu = c(2, 1),
                    Sigma = S)
cov(xy[, 1], xy[, 2])
cor(xy[, 1], xy[, 2]) # expect -1/sqrt(18)
# Found general formula: a = -cov(X, Y) / var(Y)
a <- -var(xy[, 1]) / cov(xy[, 1], xy[, 2])
cov(xy[, 1], xy[, 1] + c * xy[, 2]) # expect 0
true_a <- -S[1, 1] / S[1, 2]


# X,Y is a random point from triangle with vertices (0,0), (-1,1), (1,1).
# Find Cov(X, Y).
#
# Note, -1 < x < 1, 0 < y < 1, and y > |x|
x <- runif(1e6, -1, 1)
y <- runif(1e6, 0, 1)
drop_ix <- which(y < abs(x)) ## drop points outside triangle
dat <- data.frame(
    x = x[-drop_ix],
    y = y[-drop_ix]
)
xlin <- seq(-1, 1, length.out = 1000)
# plot triangle and sampled points
plot(dat[sample.int(nrow(dat), 10000), ], pch = ".")
lines(xlin, abs(xlin), col = "red")
lines(xlin, rep_len(1, length(xlin)), col = "red")
# E[X] = 0, E[Y] = 2/3, E[XY] = 0, Cov(X,Y) = 0
mean(dat$x) # expect 0
mean(dat$y) # expect 2/3
mean(dat$x * dat$y) # expect 0
cov(dat$x, dat$y) # expect 0
cor(dat$x, dat$y) # expect 0



# Let U1, U2, ..., be an iid seq of Unif(0, 1) r.v.
# Let N >= 2, be the first numbers s.t. U1 > U2 > ... U_{N-1} < U_N.
#  In other words, sample Unif(0, 1)'s until the sampled value is larger than
#  the previously sampled value. Then stop.
# Find E[N]
Ns <- numeric(1e5)
for (iter in seq_along(Ns)) {
    u <- runif(30) # very low prob that N > 30, safe to cutoff there
    N <- 1
    for (i in seq(2, length(u_order))) {
        N <- N + 1
        if (u[i] > u[i - 1])
            break
    }
    Ns[iter] <- N
}
mean(Ns)
exp(1)
hist(Ns)




mat <- matrix(NA_integer_, nrow = 1e6, ncol = 2)
colnames(mat) <- c("N1", "N2")
for (i in 1:nrow(mat)) {
    rolls <- sample.int(6, size = 60, replace = TRUE)
    mat[i, "N1"] <- which(rolls == 1)[1]
    mat[i, "N2"] <- which(rolls == 2)[1]
}
## throw away simulations where face 1 or 2 was never rolled (technically a very small bias)
mat <- mat[!is.na(mat[, 1]) & !is.na(mat[, 2]), ]
cov(mat[, "N1"], mat[, "N2"])







# Can we construct r.v. X ~ Ber(px) and Y ~ Ber(py) s.t. Corr(X, Y) = r,
#   for some numer -1 <= r <= 1?
# Yes, we sample X ~ Ber(px) and dependent on the outcome, sample Y with
#   a given probability, which can be calculated from the correlation formula.
# This is cool since Indicator rv's are Bernoulli rvs, which means we can exactly
# set the correlation between any 2 events if we have control over the condtional
# probabilities (i.e., P(Y = 1 | X = 1) and P(Y = 1 | X = 0)).
# (Well, except for some combinations of event probabilities that are not possible.)
#
CorrelatedEvents <- function(prob_x, prob_y, corr) {
    stopifnot(
        prob_x >= 0 && prob_x <= 1,
        prob_y >= 0 && prob_y <= 1,
        abs(corr) <= 1
    )
    # compute conditional probabilities needed for the given corr
    prob_y_given_x1 <- (
        (prob_x * prob_y) +
            corr * sqrt(prob_x * (1 - prob_x) * prob_y * (1 - prob_y))
    ) / prob_x
    if (prob_y_given_x1 < 0 || prob_y_given_x1 > 1)
        stop("Parameters lead to invalid probability P(Y = 1 |X = 1): ", prob_y_given_x1)
    
    prob_y_given_x0 <- (
        prob_y - prob_x * prob_y_given_x1 
    ) / (1 - prob_x)
    if (prob_y_given_x0 < 0 || prob_y_given_x0 > 1)
        stop("Parameters lead to invalid probability P(Y = 1 |X = 0): ", prob_y_given_x0)
    
    # fn to sample X ~ Ber(prob_x), Y ~ Ber(prob_y), s.t. cor(X, Y) = corr
    .sample_fn <- function(n = 1) {
        mat <- data.frame(
            X = rbinom(n, size = 1, prob = prob_x),
            # pre-compute the probabilities
            y_given_x1 = rbinom(n, size = 1, prob = prob_y_given_x1),
            y_given_x0 = rbinom(n, size = 1, prob = prob_y_given_x0)
        )
        mat[["Y"]] <- ifelse(mat[["X"]] == 1,
                             mat[["y_given_x1"]],
                             mat[["y_given_x0"]])
        
        mat <- mat[, c("X", "Y")]
        # # the more straightforward but slower way to do this:
        # for (i in nrow(mat)) {
        #     if (mat[i, "X"] == 1)
        #         mat[i, "Y"] <- rbinom(1, size = 1, prob = prob_y_given_x1)
        #     else
        #         mat[i, "Y"] <- rbinom(1, size = 1, prob = prob_y_given_x0)
        # }
        return(mat)
    }
    
    obj <- list(
        prob_x = prob_x,
        prob_y = prob_y,
        corr = corr,
        prob_y_given_x1 = prob_y_given_x1,
        prob_y_given_x0 = prob_y_given_x0,
        sample = .sample_fn
    )
    return(obj)
}

ce <- CorrelatedEvents(0.5, 0.5, 0.5)

ce$prob_y_given_x1
ce$prob_y_given_x0

samples <- ce$sample(1e6)
cor(samples$X, samples$Y)
mean(samples$X)
mean(samples$Y)
mean(samples[samples$X == 1, "Y"])
mean(samples[samples$X == 0, "Y"])





# X ~ N(0, 1)
# theta = {-1, 1} with equal probability
# Y = |X| * theta
# Is (X,Y) bivariate normal?
# No - linear combinations of X and Y need not be normal.
X <- rnorm(1e6)
theta <- sample(c(-1, 1), size = 1e6, replace = TRUE)
Y <- abs(X) * theta
mean(Y)
var(Y)
hist(X)
hist(Y)
a <- 10
b <- 10
hist( a*X + b*Y )



# parameters of the log-normal distribution
mu <- 1
sd <- 1
X <- rlnorm(1e6, meanlog = mu, sdlog = sd)

mean(X)
exp(mu + sd^2 / 2)

var(X)
(exp(sd^2) - 1) * exp(2 * mu + sd^2)




# X ~ Gamma(2, 1)
# Y|X ~ U(0, 1/X)
# Probability Y > 2?
# E[Y]?
x <- rgamma(1e6, 2, 1)
yx <- sapply(x, \(t) runif(1, 0, 1/t))

# P(Y > 2 | X = x) = 1 - 2x, for 0 < x < 1/2 (0 otherwise)
mean(yx[x > 0.9 & x < 1] > 2)     # 0
mean(yx[x > 0.09 & x < 0.1] > 2)  # 1 - 2*0.1 = 0.8
mean(yx[x > 0.19 & x < 0.2] > 2)  # 1 - 2*0.2 = 0.6
mean(yx) # 1/2






#
# N ~ Poi(lambda). Given N = n, generate N iid points uniformly distributed on
# the unit square. Given N >= 1, find the probability that all points lie below
# the diagonal y = x.
#
lambda <- 3

zeros <- 0
res <- logical(1e5)
for (i in seq_along(res)) {
    N <- rpois(1, lambda = lambda)
    if (N == 0) {
        zeros <- zeros + 1
        next
    }
    x <- runif(N)
    y <- runif(N)
    res[i] <- all(y <= x)
}
mean(res)
exp(-lambda/2) # probability if we allow N = 0... quite different for small lambda
exp(-lambda/2) - exp(-lambda)

# probability for different lambdas - generally becomes less likely as lambda
#   increases because larger lambda -> large avg number of trials -> more likely
#   to have a point above the diagonal the more you points you have.
ls <- seq(1, 10, length.out = 100)
plot(ls, exp(-ls/2), col = "red", type = "l") # prob. if not conditioning on N>=0
lines(ls, exp(-ls/2) - exp(-ls))





# There are N hunters, where N ~ Pois(lambda).
# A flock of ducks flies by and each hunter picks a duck uniformly at random, and
# independently of all other hunters.
# Each hunter hits its target with probability p.
# Say a flock of 20 ducks flies by...
#
# The number of hunters who choose duck i is X_i.
# Then given N = n, (X_1, ..., X_20) ~ Mult(n, 1/20, ..., 1/20).
# By Poisson thinning, the unconditional distributions of the X_i are
# Pois(lambda/20), and all X_i independent.
# Now given X_i = x, the number of hunters who hit duck i is Y_i ~ Bin(x, p).
# So by Poisson thinning, Y_i ~ Pois(lambda/20 * p). Again, all Y_i mutually
# independent.
#
# What is the probability all ducks escape unhurt?
# That occurs when Y_1, ..., Y_20 all equal 0.
# We can evaluate the joint PMF, which is the produce of the Y_i Poisson PMFs,
# to find exp(-lambda/20 * p)^20.
nduck <- 20
prob_hit <- 0.5
hunter_lambda <- 3
niter <- 1e4

ducks <- 1:nduck
duck_picks <- numeric(nduck)
duck_hits <- numeric(nduck)
allmiss <- 0

for (i in 1:niter) {
    nhunter <- rpois(1, lambda = hunter_lambda)
    duck_pick <- sample(ducks, nhunter, replace = TRUE)
    duck_hit <- rbinom(nhunter, size = 1, prob = prob_hit)
    if (sum(duck_hit) == 0) {
        allmiss <- allmiss + 1
        next
    }
    
    dpt <- table(duck_pick)
    dpi <- as.integer(names(dpt))
    duck_picks[dpi] <- duck_picks[dpi] + dpt
    
    dht <- table(duck_pick[duck_hit == 1])
    dhi <- as.integer(names(dht))
    duck_hits[dhi] <- duck_hits[dhi] + dht
}

duck_picks / niter
hunter_lambda / nduck

duck_hits / niter
hunter_lambda * prob_hit / nduck

allmiss / niter
exp(-hunter_lambda / nduck * prob_hit)^nduck




#
# Roll an n-sided die 5 times.
# Let X_1 count the number of times the first face is seen, X_2 the number of
# times the second face is seen, etc.
# What is the probability that X_1 == X_2?
#
library(data.table)
n <- 10

mat <- matrix(0, nrow = 1e6, ncol = n)
for (i in 1:nrow(mat)) {
    nums <- sample(1:n, size = 5, replace = TRUE)
    ntab <- table(nums)
    ix <- as.integer(names(ntab))
    mat[i, ix] <- mat[i, ix] + ntab
}
mat <- as.data.table(mat)
mat[V1 == V2, .N] / nrow(mat)
((n-2)/n)^5 +
    (20) * (1/n)^2 * ((n-2)/n)^3 +
    (30) * (1/n)^4 * ((n-2)/n)

