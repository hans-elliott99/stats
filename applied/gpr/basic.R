library(data.table)
library(ggplot2)


# https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(24)00476-8/fulltext
# suppl appendix 2, pg 15


sim_stgpr_data <- function(
        nregion = 3,
        ncountry_per_region = 4,
        nyear = 10,
        global_mean = NA,
        global_variance = NA,
        region_means = NA,
        region_variances = NA,
        country_means = NA,
        country_variances = NA,
        seed = sample.int(1e5, 1)
) {
    set.seed(seed)
    ncountry <- nregion * ncountry_per_region
    dat <- data.table(
        region_id = rep(1:nregion, each = ncountry_per_region * nyear),
        country = rep(1:ncountry_per_region, times = nregion * nyear)
    )
    dat[, country_id := .GRP, by = .(region_id, country)]
    dat[, country := NULL]
    dat[, year_id := rep(1:nyear), by = country_id]
    
    reg2country <- unique(dat[, .(region_id, country_id)])
    
    # global
    global_mean <- if (is.na(global_mean)) {
        rnorm(1, 3, 3)
    } else {
        global_mean
    }
    global_variance <- if (is.na(global_variance)) {
        rgamma(1, 6, 2)
    } else {
        global_variance
    }
    
    # region level
    region_means <- if (any(is.na(region_means))) {
        rnorm(nregion, global_mean, global_variance)
    } else {
        region_means
    }
    stopifnot(length(region_means) == nregion)
    # within region variance
    region_variances <- if (any(is.na(region_variances))) {
        rgamma(nregion, 3, 2)
    } else {
        region_variances
    }
    stopifnot(length(region_variances) == nregion)
    
    # country level
    country_means <- if (any(is.na(country_means))) {
        MASS::mvrnorm(
            n = 1,
            mu = region_means[reg2country$region_id],
            Sigma = diag(region_variances[reg2country$region_id])
        )
    } else {
        country_means
    }
    stopifnot(nrow(country_means) == ncountry)
    # within country variance
    country_variances <- if (any(is.na(country_variances))) {
        rgamma(ncountry, 1, 2)
    } else {
        country_variances
    }
    stopifnot(length(country_variances) == ncountry)
    
    # country values
    dat[, `:=` (
        global_mean = global_mean,
        global_variance = global_variance,
        region_mean = region_means[region_id],
        region_variance = region_variances[region_id],
        country_mean = country_means[country_id],
        country_variance = country_variances[country_id]
    )]
    dat[, value := rnorm(
        .N, country_mean[1], sqrt(country_variance[1])
    ), by = country_id]
    
    return(dat[])
}


assign_missing <- function(dat, n = NULL, prop = NULL) {
    if (!is.null(n) && !is.null(prop)) stop("specify 'n' or 'prop'")
    size <- if (!is.null(n)) n else ceiling(nrow(dat) * prop)
    ix <- sample.int(nrow(dat), size = size)
    dat[, missing := FALSE]
    dat[ix, missing := TRUE]
    return(dat[])
}



dat <- sim_stgpr_data(
    nregion = 3,
    ncountry_per_region = 4,
    nyear = 10,
    seed = 123
)
dat <- assign_missing(dat, prop = 0)


ggplot(dat, aes(x = year_id, y = value, color = as.factor(country_id))) +
    geom_point(aes(shape = missing)) +
    geom_hline(
        data = dat[, .(country_mean = country_mean[1]),
                   by = .(region_id, country_id)],
        aes(yintercept = country_mean, color = as.factor(country_id)),
        linetype = "dashed"
    ) +
    geom_hline(
        data = dat[, .(region_mean = region_mean[1]),
                   by = .(country_id)],
        aes(yintercept = region_mean),
        linetype = "solid", alpha = 0.5
    ) +
    facet_wrap(~country_id, scales = "free_y")





#
# stage 1: linear model prior
#
stage1 <- lme4::lmer(
    value ~ 1 + (1 | region_id / country_id),
    data = dat[missing == FALSE]
)
summary(stage1)

dat[, stage1_value := predict(stage1, newdata = dat)]

ggplot(dat, aes(x = value, y = stage1_value, color = as.factor(country_id))) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red")

ggplot(dat, aes(x = year_id)) +
    geom_point(aes(y = value)) +
    geom_line(aes(y = stage1_value, color = "predicted country mean"),
              linetype = "dashed") +
    geom_line(
        data = unique(dat[, .(year_id, country_id, country_mean)]),
        aes(y = country_mean, color = "true country mean"),
    ) +
    geom_line(
        data = unique(dat[, .(year_id, country_id, region_mean)]),
        aes(y = region_mean, color = "true region mean"),
    ) +
    geom_line(
        data = unique(dat[, .(year_id, country_id, global_mean)]),
        aes(y = global_mean, color = "true global mean"),
    ) +
    facet_wrap(~country_id, scales = "free_y") +
    scale_color_manual(
        values = c("predicted country mean" = "black",
                   "true country mean" = "red",
                   "true region mean" = "blue",
                   "true global mean" = "green")
    )




#
# stage 2: space-time smoothing
#
setorder(dat, year_id)
dat[, residual := value - stage1_value]


weights_time_smooth <- function(tvec, lambda) {
    #' tvec: time vector
    tops <- numeric(length(tvec))
    for (t in seq_along(tvec)) {
        tops[t] <- abs(tvec[t] - tvec[1])
    }
    w <- (1 - ( tops / (1 + which.max(tops)) ) ^ lambda)^3
    w <- c(rev(w), w[-1])
    return(w)
}

time_smooth_resid <- function(tvec, rvec, lambda = NULL, wvec = NULL) {
    if (is.null(wvec)) {
        if (is.null(lambda)) stop("specify 'lambda' or 'wvec'")
        wvec <- weights_time_smooth(tvec, lambda)
    }
    stopifnot(length(tvec) == length(rvec))
    stopifnot(length(wvec) == 2 * length(rvec) - 1)
    # find center of weights
    center <- which(wvec == max(wvec))
    # calculate weighted residuals
    out <- numeric(length(rvec))
    for (i in seq_along(tvec)) {
        # generate sliding window
        window_ix <- seq(center - tvec[i] + 1, length(wvec) - tvec[i] + 1)
        # normalize weights
        wp <- wvec[window_ix] / sum(wvec[window_ix])
        # weighted sum
        out[i] <- drop(wp %*% rvec) # sum(wp * rvec)
    }
    return(out)
}


### example
# higher lambda means more smoothing across time
wvec <- weights_time_smooth(seq(1, 10), lambda = 0.3)
plot(wvec)
rvec <- dat[country_id == 1, residual]
tvec <- dat[country_id == 1, year_id] 
pvec <- dat[country_id == 1, stage1_value]
out <- time_smooth_resid(tvec, rvec, wvec = wvec)

plot(tvec, pvec + rvec, type = "l")
points(tvec, pvec + out, col = "blue")
### 

time_weights <- weights_time_smooth(seq(1, max(dat$year_id)), lambda = 0.7)

setorder(dat, year_id)
dat[, resid_smooth := time_smooth_resid(year_id, residual, wvec = time_weights),
    by = country_id]

dat[country_id %in% 1:9,
    .(year_id, country_id, value, stage1_value,
      s2 = stage1_value + resid_smooth)] |>
    ggplot(aes(x = year_id)) +
    geom_line(aes(y = value, color = "raw data")) +
    geom_line(aes(y = stage1_value, color = "stage 1")) +
    geom_line(aes(y = s2, color = "stage 2")) +
    facet_wrap(~country_id, scales = "free_y")


# spatial smoothing...



# GPR






























