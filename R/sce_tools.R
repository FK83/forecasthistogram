#' @importFrom magrittr %>%
#' @importFrom dplyr if_else
#' @importFrom stats quantile integrate optim optimize pnorm rexp uniroot
#' @importFrom graphics lines matplot points
#' @importFrom utils head globalVariables

# define global variables (merely to avoid note when checking the package)
utils::globalVariables(c("x1", "x2", "y1", "y2"))

# handles NAs and transforms non-NAs to [0,1] interval
# (used for raw SCE data)
transform_p <- function(x){
  out <- x
  out[is.na(x)] <- 0
  (out/100)
}

#' Functions for the (Expected) Ranked Probability Score
#'
#' @param x object of class `forecasthistogram', or vector of probabilities (must be positive and sum to one)
#' @param kstar, realizing category (integer, min = 1, max = length(p))
#' @return \code{erps} returns the Expected Ranked Probability Score (ERPS). \code{rps} returns the Ranked Probability Score (RPS). \code{ebs} returns the expected Brier score (EBS).
#' @details Krüger and Pavlova (2020) propose to use the ERPS as a robust measure
#' of the uncertainty in a survey histogram. The RPS is due to Epstein (1969).
#' See also Section 4.2 in Krüger and Pavlova (2020).
#' @references
#' Epstein, E.S. (1969): `A Scoring System for Probability Forecasts of Ranked Categories', Journal of Applied Meteorology and Climatology 8, 985-987.
#'
#' Krüger, F., and L. Pavlova (2020): `Quantifying Subjective
#' Uncertainty in Survey Expectations'. Working Paper, \doi{10.18452/21401}.
#' @name rps_tools
#' @aliases rps
#' @aliases erps
#' @aliases ebs
NULL

#' @export
erps <- function(x){
  UseMethod("erps")
}
#' @export
ebs <- function(x){
  UseMethod("ebs")
}

#' @rdname rps_tools
#' @export
erps.numeric <- function(x){
  check_p(x)
  P <- cumsum(x)
  sum(P*(1-P))
}
#' @export
#' @rdname rps_tools
erps.forecasthistogram <- function(x){
  erps.numeric(x$p)
}

#' @rdname rps_tools
#' @export
ebs.numeric <- function(x){
  check_p(x)
  sum(x*(1-x))
}
#' @export
#' @rdname rps_tools
ebs.forecasthistogram <- function(x){
  ebs.numeric(x$p)
}




#' Helper function to convert numerical realization into bin indicator
#'
#' @export
#' @param y numerical value
#' @param ub vector of upper bounds of histogram bin ranges
num_to_bin <- function(y, ub){
  check_ub(ub, SCE = TRUE)
  min(which(y <= ub))
}

#' @export
rps <- function(x, kstar){
  UseMethod("rps")
}

#' @export
#' @rdname rps_tools
rps.numeric <- function(x, kstar){
  P <- cumsum(x)
  K <- length(x)
  if (kstar == 1){
    sum((1-P)^2)
  } else {
    ( sum(P[1:(kstar-1)]^2) + sum((1-P[kstar:K])^2) )
  }
}

#' @export
#' @rdname rps_tools
rps.forecasthistogram <- function(x, kstar){
  rps.numeric(x$p, kstar)
}

# Fast version of rps (here p is a matrix of probability forecasts)
rps_mat <- function(p, kstar){
  P <- apply(p, 1, cumsum) %>% t
  K <- ncol(p)
  if (kstar == 1){
    rowSums((1-P)^2)
  } else {
    ( rowSums(P[,1:(kstar-1),drop = FALSE]^2) + rowSums((1-P[,kstar:K,drop = FALSE])^2) )
  }
}

#' Function to visualize SCE histogram forecast using ggplot2
#'
#' @import ggplot2
#' @export
#' @param x object of class `forecasthistogram'
#' @param ... additional parameters (currently not in use)
#' @param outer_width support extension in case of infinite limits
#' @param ylim support limit for vertical axis (defaults to \code{NULL}, useful to make several comparable plots)
#' @return plot (ggplot2 object)
plot.forecasthistogram <- function(x, ylim = NULL, outer_width = 3, ...){
  ub <- x$ub
  p <- x$p
  # length of p vector
  lp <- length(p)
  # replace rightmost limit if it is infinite
  if (ub[lp] == Inf){
    ub[lp] <- ub[lp-1] + outer_width
  }
  # vector of bin widths
  bin_widths <- c(outer_width, diff(ub))
  # Transform probs (divide by bin width)
  p_t <- p/bin_widths
  # Auxiliary data frame
  d <- data.frame(x1 = c(ub[1] - outer_width, ub[-lp]),
                  x2 = c(ub[-lp], ub[lp]),
                  y1 = rep(0, lp),
                  y2 = p_t)
  # Make plot
  g <- ggplot(d, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2)) +
    geom_rect(col = "black", alpha = .5) + theme_minimal() +
    scale_x_continuous(breaks = ub[-lp]) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(size=.1, color="black"),
          text = element_text(size = 20))
  if (!is.null(ylim)){
    g <- g + ylim(ylim)
  }
  # If quantification is available for x: Add curve to plot
  if (!is.null(x$method)){
    # Add curve
    if (x$method == "gbeta"){
      ff <- function(z){
        # density of generalized beta distribution
        ret <- auxf(z, a = x$parameters[1], b = x$parameters[2],
                    l = x$parameters[3], r = x$parameters[4])
        # assign zero for NAs and values outside support
        zero_sel <- is.na(ret) | (z < x$parameters[3]) |
          (z > x$parameters[4])
        ret[zero_sel] <- 0
        # return vector
        ret
      }
    } else {
      ff <- function(z){
        dtri(z, a = x$parameters[1], b = x$parameters[2],
             c = x$parameters[3])
      }
    }
    g <- g +
      stat_function(fun = ff, n = 1000, geom = "line",
                    size = I(1.4)) + ylab("")
  }
  g
}

# Function to numerically recover quantile from CDF function
q_from_cdf <- function(cdf, alpha, lims){
  auxd <- function(x){
    (cdf(x)-alpha)^2
  }
  out_tmp <- optimize(auxd, interval = lims)
  if (out_tmp$objective > 1e-6){
    stop("Failed to find quantile")
  } else {
    out_tmp$minimum
  }
}

# Vectorized version of function q_from_cdf
get_q <- function(cdf, lims, quantile_levels = c(.25, .5, .75)){
  q <- rep(NA, length(quantile_levels))
  for (jj in seq_along(quantile_levels)){
    q[jj] <- q_from_cdf(cdf, alpha = quantile_levels[jj], lims = lims)
  }
  q
}

# Helper function to check whether non-zero histogram bins are next to each other
is_adjacent <- function(p){
  aux <- which(p != 0)
  out <- if_else(all(diff(aux) == 1), TRUE, FALSE)
  out
}

#' @export
quantify <- function(x, ...){
  UseMethod("quantify")
}

#' Fit parametric distribution to histogram probabilities
#'
#' @export
#' @param x object of class `forecasthistogram', or vector of probabilities (must be positive and sum to one)
#' @param ub vector of upper bounds
#' @param fit_support whether to choose support according to statistical criterium (default is \code{TRUE})
#' @param support_limit defaults to 38 (relevant only if \code{fit_support == FALSE})
#' @param ... other inputs (currently not in use)
#' @return object of class `forecasthistogram'. Use function \code{quantile.forecasthistogram} to compute quantiles, \code{plot.forecasthistogram} for plotting, and
#' \code{mean.forecasthistogram} for computing the implied mean.
#' @details The function fits a parametric distribution based on the histogram probabilities in \code{p}. The fitting procedure has been proposed by
#' Engelberg et al. (2009). Briefly, the method fits a generalized beta distribution if \code{p} features 3 or more bins,
#' and fits a simpler triangular distribution if \code{p} features at most two bins. See Krüger and Pavlova (2020, Appendix A) for details on the
#' present implementation.
#' @references
#' Engelberg, J., Manski, C.F., and J. Williams (2009): `Comparing the Point Predictions
#' and Subjective Probability Distributions of Professional Forecasters', Journal of Business & Economic Statistics 27,
#' 30-41, \doi{10.1198/jbes.2009.0003}
#'
#' Krüger, F., and L. Pavlova (2020): `Quantifying Subjective
#' Uncertainty in Survey Expectations'. Working Paper, \doi{10.18452/21401}.
#' @name quantify

#' @export
#' @rdname quantify
quantify.numeric <- function(x, ub, fit_support = TRUE, support_limit = 38, ...){
  nr_bins <- sum(x != 0)
  if (nr_bins < 3){
    #if (!is_adjacent(p)){
      #stop("In the two-bin case, EWM procedure requires adjacent bins")
    #}
    aux <- fit_hist_triangle2(ub = ub, p = x)
    mth <- "triangle"
  } else {
    aux <- fit_hist_gb(ub = ub, p = x, fit_support, support_limit)
    mth <- "gbeta"
  }
  aux$method <- mth
  aux$ub <- ub
  aux$p <- x
  aux$fit_support <- fit_support
  structure(aux, class = "forecasthistogram")
}

#' @export
#' @param x object of class `forecasthistogram'
#' @rdname quantify
quantify.forecasthistogram <- function(x, fit_support = TRUE, support_limit = 38, ...){
  quantify.numeric(x = x$p, ub = x$ub, fit_support, support_limit)
}

#' Compute quantiles of fitted histogram
#' @export
#' @param x object of class \code{forecasthistogram}
#' @param probs quantile levels
#' @param ... additional parameters, currently not in use
quantile.forecasthistogram <- function(x, probs = (1:9)/10, ...){
  # if x has not been quantified yet
  if (is.null(x$method)){
    x <- quantify.forecasthistogram(x)
    message("Histogram quantified for quantile computation")
  }
  # choose limits
  # try "narrow" limits
  q_lim1 <- max(c(x$support[1], x$moments[1]-5*x$moments[2]))
  q_lim2 <- min(c(x$support[2], x$moments[1]+5*x$moments[2]))
  check1 <- x$CDF(q_lim1) < min(probs)
  check2 <- x$CDF(q_lim2) > max(probs)
  # use "narrow" limits if wide enough, else use full support
  if (check1 & check2){
    q_lims = c(q_lim1, q_lim2)
  } else {
    q_lims = x$support
  }
  tmp <- get_q(cdf = x$CDF, lims = q_lims,
               quantile_levels = probs)
  names(tmp) <- paste0(probs*100, "%")
  tmp
}

#' Compute mean of fitted histogram
#' @export
#' @param x object of class \code{forecasthistogram}
#' @param ... additional parameters, currently not in use
mean.forecasthistogram <- function(x, ...){
  # if x has not been quantified yet
  if (is.null(x$method)){
    x <- quantify.forecasthistogram(x)
    message("Histogram quantified for mean computation")
  }
  x$moments[1]
}

fit_hist_gb <- function(ub, p, fit_support = TRUE, support_limit = 38){
  # largest finite value in ub
  largest_finite <- max(ub[ub < Inf])
  # Get support of histogram
  aux <- get_hist_support(ub, p)
  l <- aux$l
  r <- aux$r
  # Optional: Fix support
  # (kicks in only if support is unknown, i.e. null)
  if (!fit_support){
    if (is.null(l)) l <- -support_limit
    if (is.null(r)) r <- support_limit
  }
  # Estimate free parameters
  # Four scenarios, depending on whether left and/or right
  # limit of support is estimated
  if (is.null(l) & is.null(r)){
    # both limits unknown
    ff <- function(th){
      tmp <- pgbeta(ub, a = exp(th[1]) + 1,
                    b = exp(th[2]) + 1, l = th[3], r = th[4])
      mean((cumsum(p) - tmp)^2)
    }
    th0 <- c(0, 0, ub[1]-2, largest_finite + 2)
    cs <- 1
  } else if (is.null(l) & !is.null(r)){
    # left limit unknown, right limit known
    ff <- function(th){
      tmp <- pgbeta(ub, a = exp(th[1]) + 1,
                    b = exp(th[2]) + 1, l = th[3], r = r)
      mean((cumsum(p) - tmp)^2)
    }
    th0 <- c(0, 0, ub[1] - 2)
    cs <- 2
  } else if (!is.null(l) & is.null(r)){
    # left limit known, right limit unknown
    ff <- function(th){
      tmp <- pgbeta(ub, a = exp(th[1]) + 1,
                    b = exp(th[2]) + 1, l = l, r = th[3])
      mean((cumsum(p) - tmp)^2)
    }
    th0 <- c(0, 0, largest_finite + 2)
    cs <- 3
  } else if (!is.null(l) & !is.null(r)){
    # both limits known
    ff <- function(th){
      tmp <- pgbeta(ub, a = exp(th[1]) + 1,
                    b = exp(th[2]) + 1,
                    l = l, r = r)
      mean((cumsum(p) - tmp)^2)
    }
    th0 <- c(0, 0)
    cs <- 4
  }
  # Can be used for debugging (print parameter guesses)
  #ff2 <- function(par){
    #par %>% paste(collapse = "_") %>% print
    #return(ff(par))
  #}
  # Modify objective function
  # (Assign penalty when beta dist parameters get too large)
  # Following Armantier et al (2017, Appendix C): Set
  # maximal support limit to +/- 38
  # Also impose support limit greater than 12 since
  # the SCE bins end there
  if (cs == 1){
    # both limits unknown
    ff2 <- function(th){
      if_else(any(th[1:2] > 3.5) |
                (th[3] > ub[1]) |
                (th[4] < largest_finite) |
                any(abs(th[3:4]) > support_limit),
              10, ff(th))
    }
  } else if (cs == 2){
    # only left limit unknown
    ff2 <- function(th){
      if_else(any(th[1:2] > 3.5) |
                (th[3] > ub[1]) |
                (abs(th[3]) > support_limit),
              10, ff(th))
    }
  } else if (cs == 3){
    # only right limit unknown
    ff2 <- function(th){
      if_else(any(th[1:2] > 3.5) |
                (th[3] < largest_finite) |
                (th[3] > support_limit), 10, ff(th))
    }
  } else if (cs == 4){
    # both limits known
    ff2 <- function(th){
      if_else(any(th[1:2] > 3.5), 10, ff(th))
    }
  }

  # Run optimizer
  aux <- optim(th0, ff2)$par
  # Throw error if optimal parameters equal starting values
  if (any(aux == th0)){
    stop("Optimization seems to have failed")
  }
  # Compute transformed parameters and support limits
  if (cs == 1){
    pms <- c(exp(aux[1:2]) + 1, aux[3:4])
  } else if (cs == 2){
    pms <- c(exp(aux[1:2]) + 1, aux[3], r)
  } else if (cs == 3){
    pms <- c(exp(aux[1:2]) + 1, l, aux[3])
  } else if (cs == 4){
    pms <- c(exp(aux[1:2]) + 1, l, r)
  }
  # CDF
  CDF <- function(x){
    pgbeta(x, a = pms[1], b = pms[2],
           l = pms[3], r = pms[4])
  }
  # Support
  support <- c(pms[3], pms[4])
  # Moments
  moments <- moments_gb(pms[1], pms[2], pms[3], pms[4]) %>%
    unlist %>% unname
  # Return
  list(parameters = pms, CDF = CDF,
       support = support, moments = moments)
}

# Helper function for beta distribution
auxf <- function(z, a, b, l, r){
  d <- (beta(a,b)*((r-l)^(a+b-1)))
  out <- ((z-l)^(a-1))*((r-z)^(b-1))/d
  out
}

# CDF of generalized Beta distribution
# (Engelberg et al, Equation 1)
pgbeta <- function(x, a, b, l, r){
  out <- rep(0, length(x))
  ind1 <- which((x > l) & (x <= r))
  ind2 <- which(x > r)
  out[ind2] <- 1
  for (jj in ind1){
    tmp <- integrate(auxf, a = a, b = b,
                     l = l, r = r, lower = l,
                     upper = x[jj])$value
    out[jj] <- tmp
  }
  # Cut off at one
  out[out > 1] <- 1
  out
}

plot_gb <- function(a, b, l, r){
  rg <- seq(from = l, to = r, length.out = 1001)
  y <- auxf(rg, a, b, l, r)
  plot(rg, y, bty = "n", type = "l",
       xlab = "x", ylab = "f(x)")
}

moments_gb <- function(a, b, l, r){
  ff1 <- function(z) z*auxf(z, a = a, b = b, l = l, r = r)
  m <- integrate(ff1, lower = l, upper = r)$value
  ff2 <- function(z) (z^2)*auxf(z, a = a, b = b, l = l, r = r)
  v <- integrate(ff2, lower = l, upper = r)$value - m^2
  list(m = m, s = sqrt(v))
}

get_hist_support <- function(ub, p){
  # Lowest and highest bins used
  supp <- which(p != 0) %>% range

  if (supp[1] == 1){
    # left support limit unknown
    l <- NULL
  } else {
    # left support limit known
    l <- ub[supp[1] - 1]
  }
  if (supp[2] == length(ub)){
    # right support limit unknown
    r <- NULL
  } else {
    # right support limit known
    r <- ub[supp[2]]
  }

  list(l = l, r = r)

}

moments_mam <- function(mp, p){
  m <- sum(p*mp)
  s <- (sum(p*(mp^2)) - m^2) %>% sqrt
  list(m = m, s = s)
}

mp <- function(ub, support_limit = 38){
  cbind(c(-support_limit, ub[-length(ub)]),
        c(ub[-length(ub)], support_limit)) %>%
    rowMeans
}


triangle_get_b <- function(l, m, alpha){
  m + (m-l)*(1-alpha+sqrt(2*(1-alpha)))/(1+alpha)
}

triangle_get_a <- function(r, m, alpha){
  (m - (r-m)*(alpha + sqrt(2*alpha))/(2-alpha))
}

# Support limit follows Armantier et al (2017, Appendix C)
fit_hist_triangle <- function(ub, p, support_limit = 38){
  ub[ub == Inf] <- support_limit
  sel <- which(p != 0)
  psel <- p[sel]
  n_bins <- sum(p != 0)
  n_max <- length(p)
  if(n_bins == 1){
    if (sel > 1){
      a <- ub[sel - 1]
    } else {
      a <- -support_limit
    }
    b <- ub[sel]
  } else if (n_bins == 2){
    if (psel[1] < .5){
      b <- ub[sel[2]]
      a <- triangle_get_a(b, ub[sel[1]], psel[1])
    } else {
      if (sel[1] == 1){
        a <- -support_limit
      } else {
        a <- ub[sel[1]-1]
      }
      b <- triangle_get_b(a, ub[sel[1]], psel[1])
    }
  }
  c <- mean(c(a, b))
  parameters <- c(a, b, c)
  CDF <- function(x){
    ptri(x, a, b, c)
  }
  moments <- moments_tri(a, b, c)
  list(parameters = parameters,
       CDF = CDF,
       support = c(a, b),
       moments = moments)
}

# Alternative procedure for triangle fitting (based on optimization)
fit_hist_triangle2 <- function(ub, p, support_limit = 38){
  ub[ub == Inf] <- support_limit
  sel <- which(p != 0)
  psel <- p[sel]
  n_bins <- sum(p != 0)
  n_max <- length(p)
  # single-bin case: need to fix support limit
  if (n_bins == 1){
    if (sel > 1){
      a <- ub[sel - 1]
    } else {
      a <- -support_limit
    }
    b <- ub[sel]
  } else if (n_bins == 2){
    # two-bin case, adjacent bins: need to fix one limit
    # (otherwise solution is not unique)
    if (is_adjacent(p)){
      # find nonzero bins
      which_nz <- which(p != 0)
      # find width of bins !!! Assume that left endpoint of first interval
      # is at -Inf
      bin_width <- ub - c(-Inf, head(ub, -1))
      # fix endpoint of bin with higher density
      bin_dens <- (p/bin_width)
      # find bin with highest density
      which_max <- which.max(bin_dens)
      if (sum(bin_dens == max(bin_dens)) == 2){
        warning("Highest-density bin not unique. Arbitrarily chose leftmost bin")
      }
      if (which_max == which_nz[1]){
        # fix left endpoint
        a <- ub[which_nz[1]-1]
        # find right endpoint
        b <- optimize(f = ff_b, interval = c(a, 1e4),
                      ub = ub, p = p, a = a)$minimum
      } else if (which_max == which_nz[2]){
        # fix right endpoint
        b <- ub[which_nz[2]]
        # find left endpoint
        a <- optimize(f = ff_a, interval = c(-1e4, b),
                      ub = ub, p = p, b = b)$minimum
      }
    } else {
      # two-bin case, non-adjacent bins: optimize both endpoints
      # find starting values
      if (sel[1] == 1){
        par1 <- -support_limit
      } else {
        par1 <- ub[sel[1]-1]
      }
      par2 <- sqrt(ub[sel[2]-1]-par1)
      # run optimizer
      opt <- optim(par = c(par1, par2), fn = ff_ab,
                   ub = ub, p = p)
      stopifnot(opt$converence == 0)
      a <- opt$par[1]
      b <- a + opt$par[2]^2
    }
  }
  # c fixed (isosceles triangle)
  c <- mean(c(a, b))
  # summary information on fit
  parameters <- c(a, b, c)
  CDF <- function(x){
    ptri(x, a, b, c)
  }
  moments <- moments_tri(a, b, c)
  list(parameters = parameters,
       CDF = CDF,
       support = c(a, b),
       moments = moments)
}

# function for optimization of triangular distribution
# case of two parameters (left and right endpoint)
ff_ab <- function(theta, ub, p){
  a <- theta[1]
  b <- theta[1] + (theta[2]^2)
  c <- mean(c(a, b))
  d <- ptri(ub, a, b, c) - cumsum(p)
  sum(d^2)
}

# function for optimization of triangular distribution
# case of one parameter (left endpoint a, right one is fixed at b)
ff_a <- function(theta, ub, p, b){
  a <- theta
  c <- mean(c(a, b))
  d <- ptri(ub, a, b, c) - cumsum(p)
  sum(d^2)
}

# function for optimization of triangular distribution
# case of one parameter (right endpoint b, left one is fixed at a)
ff_b <- function(theta, ub, p, a){
  b <- theta
  c <- mean(c(a, b))
  d <- ptri(ub, a, b, c) - cumsum(p)
  sum(d^2)
}

#' Density and CDF of triangular distribution
#'
#' @param x evaluation point for density or CDF
#' @param a,b,c parameters of distribution (left limit, center, right limit)
#' @return \code{dtri} returns the density, \code{ptri} returns the CDF,
#' \code{moments_tri} returns the mean and standard deviation
#' @name triangular
NULL

#' @export
#' @rdname triangular
dtri <- function(x, a, b, c) {
  out <- rep(0, length(x))
  sel1 <- x >= a & x <= c
  out[sel1] <- 2*(x[sel1]-a)/((b-a)*(c-a))
  sel2 <- x > c & x <= b
  out[sel2] <- 2*(b-x[sel2])/((b-a)*(b-c))
  out
}

#' @export
#' @rdname triangular
ptri <- function(x, a, b, c) {
  out <- rep(0, length(x))
  sel1 <- x >= a & x <= c
  out[sel1] <- ((x[sel1]-a)^2)/((b-a)*(c-a))
  sel2 <- x > c & x <= b
  out[sel2] <- 1 - ((b-x[sel2])^2)/((b-a)*(b-c))
  out[x > b] <- 1
  out
}

#' @export
#' @rdname triangular
moments_tri <- function(a, b, c){
  m <- (a+b+c)/3
  v <- (a^2+b^2+c^2 -a*b-a*c-b*c)/18
  c(m, sqrt(v))
}

# Compute bounds on median, mean and mode, as implied by vector of probabilities p
compute_bounds <- function(p, lb, ub){
  # Median bound
  ind1 <- min(which(cumsum(p) >= .5))
  ind2 <- min(which(cumsum(p) > .5))
  med_bound <- c(lb[ind1], ub[ind2])
  # Mean bound
  p_sel <- p > 0
  mean_bound <- c(sum(p[p_sel]*lb[p_sel]), sum(p[p_sel]*ub[p_sel]))
  # Mode bound
  ind3 <- which.max(p)
  mode_bound <- c(lb[min(ind3)], ub[max(ind3)])
  # Return
  list(median = med_bound, mean = mean_bound, mode = mode_bound)
}

# check whether point forecast m is within bounds implied by (p, lb, ub)
within_bounds <- function(m, p, lb, ub){
  bound_list <- compute_bounds(p, lb, ub)
  sapply(bound_list, function(z) (m >= z[1] & m <= z[2])) %>% any
}

fit_moments <- function(ub, p, fit_support = TRUE){
  bins_used <- sum(p > 0)
  midp <- mp(ub)
  if (bins_used == 1){
    out <- c(midp[which(p > 0)], 1e-6)
  } else if (bins_used == 2){
    out <- moments_mam(midp, p) %>% unlist %>% unname
  } else {
    out <- fit_hist_gb(ub = ub, p = p, fit_support = fit_support)$moments
  }
  out
}

# function to check vector of probabilities
check_p <- function(p, SCE = FALSE){
  if (!is.vector(p)){
    stop("p must be a vector")
  }
  if (any(is.na(p))){
    stop("NAs in p cannot be handled")
  }
  if (any(p < 0) | any(p > 1)){
    stop("p contains values outside [0,1] interval")
  }
  if (abs(sum(p)-1) > 1e-10){
    stop("Values in p do not sum to one")
  }
  if (SCE){
    if (length(p) != 10){
      stop("p must contain 10 elements, in line with SCE question format")
    }
  }
}

# function to check vector of upper bounds
check_ub <- function(ub, SCE = FALSE){
  if (SCE){
    ub_sce <- c(-12, -8, -4, -2, 0, 2, 4, 8, 12, Inf)
    if (any(ub != ub_sce)){
      stop("Upper bounds do not conform to SCE format")
    }
  }
  if (any(diff(ub) <= 0)){
    stop("Upper bounds must be strictly increasing")
  }
}

#' Constructor function for forecasthistogram object
#'
#' @export
#' @param p vector of probabilities
#' @param ub vector of upper bounds
forecasthistogram <- function(p, ub){
  validate_forecasthistogram(new_forecasthistogram(p, ub))
}

new_forecasthistogram <- function(p, ub){
  x <- list(ub = ub, p = p)
  structure(x, class = "forecasthistogram")
}

validate_forecasthistogram <- function(x){
  # check vector of upper bounds
  check_ub(x$ub)
  # check vector of probabilities
  check_p(x$p)
  # return x
  x
}

#' Draw example histogram
#'
#' @param style Either "Gaussian" or "Random". "Gaussian" means that histogram has (approximate) bell shape.
#' @param n_bins Number of bins that receive strictly positive probability. Defaults to \code{NULL} (random number of bins). Relevant only if \code{style == "Random"}.
#' @param enforce_adjacent Whether to enforce adjacency of bins (defaults to \code{TRUE}). Relevant only if \code{style == "Random"}.
#' @details Bin definitions of random histogram are as in the Survey of Consumer Expectations data set.
#' @examples
#' # Draw random histogram
#' x <- forecasthistogram_example()
#' # Plot it
#' plot(x)
#' @export
forecasthistogram_example <- function(style = "Gaussian",
                                      n_bins = NULL,
                                      enforce_adjacent = TRUE){
  match.arg(style, c("Gaussian", "Random"))
  # Use SCE bin definitions (upper bounds)
  ub <- c(-12, -8, -4, -2, 0, 2, 4, 8, 12, Inf)
  n_total <- length(ub)
  if (style == "Gaussian"){
    ub2 <- ub[is.finite(ub)]
    n_total2 <- n_total - sum(is.infinite(ub))
    P <- pnorm(ub, mean = round(ub2[sample.int(n_total2, size = 1)]),
               sd = rexp(1, rate = 1))
    p <- c(P[1], diff(P))
  } else {
    # Determine number of bins with positive probability
    if (is.null(n_bins)){
      n_bins <- sample.int(n_total, size = 1)
    }
    # Choose positive bins
    pos_bins <- sample.int(n_total, size = n_bins, replace = FALSE) %>%
      sort
    if (enforce_adjacent){
      pos_bins <- pos_bins[1]:(pos_bins[1] + n_bins - 1)
    }
    # Draw probabilities
    p_aux <- rexp(n_bins) %>% (function(z) z/sum(z))
    p <- rep(0, n_total)
    p[pos_bins] <- p_aux
  }

  # Return forecasthistogram object
  forecasthistogram(p = p, ub = ub)
}

compare_fit <- function(ub, p, F, plot = FALSE){
  P <- cumsum(p)
  Fmatch <- F(ub)
  if (plot){
    matplot(ub, cbind(P, Fmatch), type = "n", bty = "n",
            ylab = "Cumulative probability")
    points(ub, P)
    lines(ub, Fmatch)
  }
  comp <- cbind(ub, P, Fmatch)
  objective <- sum((comp[,3]-comp[,2])^2)
  list(comp = comp, objective = objective)
}
