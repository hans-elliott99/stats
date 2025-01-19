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







