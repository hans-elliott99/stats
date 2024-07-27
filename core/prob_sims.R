#
# simulations from probability class
#


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

































