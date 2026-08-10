################################################################################
############# Wrapper function for simulation and checks #######################
################################################################################

# check if all necessary libraries are installed
# List of required packages
required_packages <- c("ggplot2", "dplyr")

# Check which are not installed
not_installed <- required_packages[
  !(required_packages %in% installed.packages()[, "Package"])
]

# Install missing packages
if (length(not_installed) > 0) {
  message(
    "Installing missing packages: ",
    paste(not_installed, collapse = ", ")
  )
  install.packages(not_installed, dependencies = TRUE)
}

# Load all packages
lapply(required_packages, library, character.only = TRUE)

message("All required packages are installed and loaded.")

# define restrict
restrict <- c("none", "weak", "moderate", "strong")

# set seed for reproducibility
set.seed(58294)

# Helper: variance with NA handling + finite guard
.var_non_na <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) < 2L) {
    return(NA_real_)
  }
  var(v)
}

# Helper: count unique finite values
.n_unique_finite <- function(v) {
  length(unique(v[is.finite(v)]))
}

# Validator function for dgp() output
validate_dgp_result <- function(
  res,
  beta1,
  restrict,
  n,
  formula,
  require_binary_for = NULL,
  min_random_n = 1000,
  min_var = 1e-8,
  min_unique_y_numeric_x = 4
) {
  # Global structure checks
  loc <- sprintf(
    "[beta1=%s, restrict=%s]",
    as.character(beta1),
    as.character(restrict)
  )
  if (!is.list(res)) {
    stop("dgp() must return a list ", loc)
  }
  expected <- c("convenience_x", "convenience_y")
  if (!setequal(names(res), expected) || length(res) != length(expected)) {
    stop(
      "dgp() must return a named list with exactly: ",
      paste(expected, collapse = ", "),
      " ",
      loc
    )
  }

  # Per-dataset structural checks
  for (nm in expected) {
    df <- res[[nm]]
    if (!is.data.frame(df)) {
      stop(sprintf("dgp()[[%s]] must be a data.frame %s", nm, loc))
    }

    f_exp <- formula(terms(as.formula(formula), data = df))
    response <- as.character(f_exp[[2]])
    predictors <- attr(terms(f_exp), "term.labels")
    rhs_vars <- all.vars(f_exp[[3]])
    missing_vars <- setdiff(rhs_vars, names(df))
    if (length(missing_vars) > 0) {
      stop(sprintf(
        "The formula references variable(s) not present in `data`: %s",
        paste(missing_vars, collapse = ", ")
      ))
    }

    if (!("y" %in% names(df))) {
      stop(sprintf("dgp()[[%s]] is missing column: y %s", nm, loc))
    }
    if (!("x" %in% names(df))) {
      stop(sprintf("dgp()[[%s]] is missing column: x %s", nm, loc))
    }
    if (!is.numeric(df$y)) {
      stop(sprintf("`y` must be numeric in dgp()[[%s]] %s", nm, loc))
    }
    if (all(is.na(df$y))) {
      stop("`y` contains only missing values.")
    }

    # check minimum variance in x and y
    vy <- .var_non_na(df$y)
    if (is.na(vy) || vy <= min_var) {
      stop(sprintf(
        "`y` variance too small (var=%.3g <= %.3g) in dgp()[[%s]] %s",
        vy,
        min_var,
        nm,
        loc
      ))
    }

    vx <- .var_non_na(df$x)
    if (is.na(vx) || vx <= min_var) {
      stop(sprintf(
        "`x` variance too small (var=%.3g <= %.3g) in dgp()[[%s]] %s",
        vx,
        min_var,
        nm,
        loc
      ))
    }

    # minimum unique values in y
    uy <- .n_unique_finite(df$y)
    if (uy < 4L) {
      stop(sprintf(
        "`y` must have at least 4 unique finite values (got %d) in dgp()[[%s]] %s",
        uy,
        nm,
        loc
      ))
    }

    # if x is numeric, require at least 4 unique finite values in x
    x_type <- classify_x(df$x)
    if (identical(x_type, "numeric")) {
      ux <- .n_unique_finite(df$x)
      if (ux < 4L) {
        stop(sprintf(
          "When `x` is numeric, `x` must have at least 4 unique finite values (got %d) in dgp()[[%s]] %s",
          ux,
          nm,
          loc
        ))
      }
    }
  }

  # equal-size check for convenience samples
  nx <- nrow(res$convenience_x)
  ny <- nrow(res$convenience_y)
  if (!identical(nx, ny)) {
    stop(sprintf(
      "`convenience_x` and `convenience_y` must have the same number of rows (got %d, vs %d) %s",
      nx,
      ny,
      loc
    ))
  }

  # minimum size check for random samples
  if (identical(restrict, "none")) {
    if (!isTRUE(all.equal(res$convenience_x, res$convenience_y))) {
      stop(
        "When restrict = 'none', `convenience_x` and `convenience_y` must be identical ",
        loc
      )
    }
  }

  invisible(TRUE)
}

# Helper: get unique non-NA values
.unique_non_na <- function(v) unique(v[!is.na(v)])

# Classify x: returns "numeric" or "binary" — otherwise errors
classify_x <- function(x) {
  # logical -> binary
  if (is.logical(x)) {
    return("binary")
  }
  # factor or ordered factor
  if (is.factor(x) || is.ordered(x)) {
    k <- nlevels(droplevels(x))
    if (k == 2) {
      return("binary")
    }
    stop(sprintf(
      "`x` is a factor with %d levels; only 2-level factors are allowed (binary).",
      k
    ))
  }
  # numeric (includes integer/double)
  if (is.numeric(x)) {
    vals <- .unique_non_na(x[is.finite(x)])
    if (length(vals) < 2) {
      stop("`x` must have at least 2 distinct non-missing values.")
    }
    if (length(vals) == 2) {
      return("binary")
    }
    return("numeric")
  }
  # everything else is not allowed
  stop(sprintf(
    "`x` must be numeric, logical, or a factor/ordered factor with 2 levels. Got class: %s",
    paste(class(x), collapse = "/")
  ))
}

# Check df has valid x (numeric or binary) and y (numeric), return x_type
check_x_y_types <- function(df) {
  x_type <- classify_x(df$x)
  # additional guard: if numeric, require >2 distinct values (not just 2)
  if (identical(x_type, "numeric")) {
    vals <- .unique_non_na(df$x[is.finite(df$x)])
    if (length(vals) <= 2) {
      stop(
        "`x` was classified as numeric but has ≤ 2 distinct values; this is effectively binary."
      )
    }
  }
  x_type
}

# MC over random samples from the DGP (parametric MC)
mc_parametric_random <- function(beta_true, n, n_match, formula, B, dgp_fun) {
  f <- as.formula(formula)
  est <- rep(NA_real_, B)

  for (b in seq_len(B)) {
    res <- dgp_fun(beta1 = beta_true, restrict = "none", n = n)
    df <- res$convenience_x
    # sample to match n
    idx <- if (nrow(df) >= n_match) {
      sample.int(nrow(df), n_match, replace = FALSE)
    } else {
      sample.int(nrow(df), n_match, replace = TRUE)
    }
    mf <- model.frame(f, data = df[idx, , drop = FALSE], na.action = na.omit)
    if (nrow(mf) < 2L) {
      next
    }
    fit <- try(lm(f, data = mf), silent = TRUE)
    if (inherits(fit, "try-error")) {
      next
    }
    cf <- coef(fit)
    # estimability check
    if (anyNA(cf)) {
      next
    }
    if ("x" %in% names(cf) && is.finite(cf["x"])) est[b] <- unname(cf["x"])
  }
  est[is.finite(est)]
}

# Check: is the random sample unbiased at this n
# Pass if beta_true lies in the empirical central 'level' interval of MC draws
check_random_unbiased_at_n_parametric <- function(
  beta_true,
  n,
  n_match,
  formula,
  B,
  dgp_fun,
  level = 0.9,
  abs_tol,
  rel_tol,
  min_ok = 100
) {
  mc <- mc_parametric_random(beta_true, n, n_match, formula, B, dgp_fun)
  if (length(mc) < min_ok) {
    stop(sprintf(
      "Too few MC estimates (got %d). Increase B or inspect DGP.",
      length(mc)
    ))
  }
  mc_mean <- mean(mc)
  ci <- quantile(mc, c((1 - level) / 2, 1 - (1 - level) / 2), names = FALSE)

  ok_ci <- (beta_true >= ci[1] && beta_true <= ci[2])

  if (!is.null(abs_tol) || !is.null(rel_tol)) {
    thr <- max(
      if (is.null(abs_tol)) 0 else abs_tol,
      if (is.null(rel_tol)) 0 else rel_tol * abs(beta_true)
    )
    ok <- ok_ci || (abs(mc_mean - beta_true) <= thr)
  } else {
    ok <- ok_ci
  }

  list(ok = ok, mc_mean = mc_mean, ci = ci, sd = sd(mc), B = length(mc))
}

# Overall checks function
checks <- function(
  niter,
  beta1,
  restrict,
  n,
  formula,
  dgp_fun,
  # MC preflight knobs
  B = 200L,
  level = 0.90,
  abs_tol = 0.02,
  rel_tol = 0.05,
  return_diag = FALSE
) {
  # input validation
  if (
    !is.numeric(niter) ||
      length(niter) != 1L ||
      !is.finite(niter) ||
      niter < 1 ||
      niter %% 1 != 0
  ) {
    stop("`niter` must be a single positive integer.")
  }
  if (
    !is.numeric(beta1) ||
      length(beta1) != 3 ||
      any(!is.finite(beta1)) ||
      any(beta1 == 0)
  ) {
    stop("`beta1` must be a numeric vector of length 3 with no zeros")
  }
  allowed_restrict <- c("none", "weak", "moderate", "strong")
  if (
    !is.character(restrict) ||
      length(restrict) != 4 ||
      !identical(restrict, allowed_restrict)
  ) {
    stop(
      "`restrict` must be exactly c('none','weak','moderate','strong') in this order."
    )
  }
  if (!inherits(formula, "formula")) {
    stop("`formula` must be an R formula (e.g., y ~ . or y ~ x + z).")
  }
  if (as.character(formula[[2]]) != "y") {
    stop("The response (left-hand side) of `formula` must be `y`.")
  }
  rhs_ok <- identical(as.character(formula[[3]]), ".") ||
    ("x" %in% all.vars(formula[[3]]))
  if (!rhs_ok) {
    stop(
      "The right-hand side of `formula` must be `.` or contain predictor `x` (e.g., y ~ x + z)."
    )
  }

  # helper to format warnings
  warnf <- function(fmt, ...) warning(sprintf(fmt, ...), call. = FALSE)

  # run structure + MC unbiasedness checks, collect diagnostics
  diag <- list()
  row_id <- 0L

  for (b1 in beta1) {
    for (rr in restrict) {
      # try to generate one dataset for this scenario (to get n_match).
      res <- tryCatch(
        dgp_fun(beta1 = b1, restrict = rr, n = n),
        error = function(e) e
      )

      if (inherits(res, "error")) {
        # record one row per sample type with the error, then continue
        for (nm in c("convenience_x", "convenience_y")) {
          row_id <- row_id + 1L
          diag[[row_id]] <- data.frame(
            beta1 = b1,
            restrict = rr,
            sample = nm,
            n_match = NA_integer_,
            mc_mean = NA_real_,
            ci_lo = NA_real_,
            ci_hi = NA_real_,
            sd = NA_real_,
            B = B,
            ok = FALSE,
            error_msg = sprintf("dgp() failed: %s", conditionMessage(res)),
            stringsAsFactors = FALSE
          )
        }
        warnf(
          "dgp() error: failed for beta1=%s, restrict=%s: %s",
          as.character(b1),
          rr,
          conditionMessage(res)
        )
        next
      }

      # validate_dgp_result might throw; soften to a warning and keep going if it does
      v_ok <- TRUE
      tryCatch(
        validate_dgp_result(
          res,
          beta1 = b1,
          restrict = rr,
          n = n,
          formula = formula
        ),
        error = function(e) {
          v_ok <<- FALSE
          warnf(
            "dgp() error: failed for beta1=%s, restrict=%s: %s",
            as.character(b1),
            rr,
            conditionMessage(e)
          )
        }
      )

      # if validation failed, record rows with the error for both samples
      for (nm in c("convenience_x", "convenience_y")) {
        n_match <- tryCatch(nrow(res[[nm]]), error = function(e) NA_integer_)

        f <- as.formula(formula)
        mf0 <- tryCatch(
          model.frame(f, data = res[[nm]], na.action = na.omit),
          error = function(e) NULL
        )

        if (is.null(mf0) || nrow(mf0) < 2L) {
          warnf(
            "Model frame too small/failed for beta1=%s, restrict=%s, sample=%s",
            as.character(b1),
            rr,
            nm
          )
        } else {
          fit0 <- try(lm(f, data = mf0), silent = TRUE)
          if (inherits(fit0, "try-error") || anyNA(coef(fit0))) {
            warnf(
              "Non-estimable model (error or NA coefficients) for beta1=%s, restrict=%s, sample=%s",
              as.character(b1),
              rr,
              nm
            )
          }
        }

        chk <- tryCatch(
          check_random_unbiased_at_n_parametric(
            beta_true = b1,
            n = n,
            n_match = n_match,
            formula = formula,
            B = B,
            dgp_fun = dgp_fun,
            level = level,
            abs_tol = abs_tol,
            rel_tol = rel_tol
          ),
          error = function(e) e
        )

        row_id <- row_id + 1L

        if (inherits(chk, "error")) {
          # record the error as a row and warn
          diag[[row_id]] <- data.frame(
            beta1 = b1,
            restrict = rr,
            sample = nm,
            n_match = n_match,
            mc_mean = NA_real_,
            ci_lo = NA_real_,
            ci_hi = NA_real_,
            sd = NA_real_,
            B = B,
            ok = FALSE,
            error_msg = sprintf(
              "check_random_unbiased_at_n_parametric() error: %s",
              conditionMessage(chk)
            ),
            stringsAsFactors = FALSE
          )
          warnf(
            "Check errored for beta1=%s, restrict=%s, sample=%s (n=%s): %s",
            as.character(b1),
            rr,
            nm,
            as.character(n_match),
            conditionMessage(chk)
          )
          next
        }

        # normal case: record results; if failed, warn instead of stop
        diag[[row_id]] <- data.frame(
          beta1 = b1,
          restrict = rr,
          sample = nm,
          n_match = n_match,
          mc_mean = chk$mc_mean,
          ci_lo = chk$ci[1],
          ci_hi = chk$ci[2],
          sd = chk$sd,
          B = chk$B,
          ok = chk$ok,
          error_msg = if (!chk$ok) {
            sprintf(
              paste0(
                "Biased estimate in unrestricted sample: MC mean=%.4g, %d%% MC=[%.4g, %.4g], ",
                "beta1=%.4g; appears biased at this n (by CI/tolerance)."
              ),
              chk$mc_mean,
              as.integer(100 * level),
              chk$ci[1],
              chk$ci[2],
              b1
            )
          } else {
            NA_character_
          },
          stringsAsFactors = FALSE
        )

        if (!chk$ok) {
          warnf(
            paste0(
              "Biased estimate in unrestricted sample: for beta1=%s, restrict=%s, sample=%s (n=%d): ",
              "MC mean=%.4g, %d%% MC=[%.4g, %.4g], beta1=%.4g. ",
              "Random baseline appears biased at this n (by CI and tolerance rules)."
            ),
            as.character(b1),
            rr,
            nm,
            n_match,
            chk$mc_mean,
            as.integer(100 * level),
            chk$ci[1],
            chk$ci[2],
            b1
          )
        }
      }
    }
  }

  diag_df <- do.call(rbind, diag)

  # success message if all checks passed
  all_ok <- all(diag_df$ok) && all(is.na(diag_df$error_msg))
  if (all_ok) {
    message(
      "Well done! All checks passed; you're good to go with the simulation."
    )
  }

  if (isTRUE(return_diag)) {
    return(diag_df)
  }
  invisible(diag_df)
}

# Compute Cohen's d and CI for two groups
cohens_d_ci <- function(x, y, conf.level = 0.95, na.rm = TRUE) {
  if (na.rm) {
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]
    y <- y[keep]
  }
  stopifnot(length(unique(x)) == 2)

  g <- split(y, x)
  y1 <- g[[1]]
  y2 <- g[[2]]

  n1 <- length(y1)
  n2 <- length(y2)

  # guard against tiny groups
  if (n1 < 2 || n2 < 2) {
    return(data.frame(
      d = NA_real_,
      se = NA_real_,
      ci.lower = NA_real_,
      ci.upper = NA_real_,
      n1 = n1,
      n2 = n2
    ))
  }

  # use doubles to avoid integer overflow
  n1d <- as.numeric(n1)
  n2d <- as.numeric(n2)

  m1 <- mean(y1)
  m2 <- mean(y2)
  s1 <- sd(y1)
  s2 <- sd(y2)

  sp2 <- ((n1d - 1) * s1^2 + (n2d - 1) * s2^2) / (n1d + n2d - 2)
  if (!is.finite(sp2) || sp2 <= 0) {
    return(data.frame(
      d = NA_real_,
      se = NA_real_,
      ci.lower = NA_real_,
      ci.upper = NA_real_,
      n1 = n1,
      n2 = n2
    ))
  }
  sp <- sqrt(sp2)

  # Cohen's d
  d <- (m2 - m1) / sp

  # SE and CI
  se_d <- sqrt((n1d + n2d) / (n1d * n2d) + (d^2) / (2 * (n1d + n2d - 2)))
  z <- qnorm(1 - (1 - conf.level) / 2)
  data.frame(
    d = d,
    se = se_d,
    ci.lower = d - z * se_d,
    ci.upper = d + z * se_d,
    n1 = n1,
    n2 = n2
  )
}

# Function to run one simulation
one_simulation_run <- function(beta1, restrict, n, formula, dgp_fun) {
  # helpers
  safe_scale <- function(v) {
    s <- sd(v, na.rm = TRUE)
    if (!is.finite(s) || s == 0) {
      return(rep(0, length(v)))
    }
    as.numeric(scale(v))
  }

  # when a sample can't be estimated, return a named NA vector
  make_na_params <- function(x_type = "numeric") {
    base <- c(
      coef = NA_real_,
      se = NA_real_,
      R2 = NA_real_,
      coef_std = NA_real_,
      se_std = NA_real_,
      R2_std = NA_real_
    )
    if (identical(x_type, "numeric")) {
      c(
        base,
        pearson = NA_real_,
        lower_CI_pearson = NA_real_,
        upper_CI_pearson = NA_real_
      )
    } else {
      c(base, cohen = NA_real_)
    }
  }

  # all combinations
  grid <- expand.grid(
    beta1 = beta1,
    restrict = restrict,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  # run dgp for all scenarios
  data <- lapply(seq_len(nrow(grid)), function(i) {
    b1 <- grid$beta1[i]
    rr <- grid$restrict[i]
    res <- tryCatch(
      dgp_fun(beta1 = b1, restrict = rr, n = n),
      error = function(e) e
    )

    return(res)
  })

  # name each result clearly
  names(data) <- with(grid, paste0("__beta_", beta1, "__restrict_", restrict))

  # function to run the linear regression model and extract the parameters of interest
  model_function <- function(data, formula) {
    # ensure formula object
    f_in <- as.formula(formula)
    # expand '.' using the columns in `data`
    f_exp <- formula(terms(f_in, data = data))

    # response and predictor names from the expanded formula
    response <- as.character(f_exp[[2]])
    predictors <- attr(terms(f_exp), "term.labels")

    # drop rows with NA in any variables used by the model
    used_vars <- unique(c(response, predictors))
    mf_idx <- complete.cases(data[, used_vars, drop = FALSE])
    df <- data[mf_idx, , drop = FALSE]

    # Fit unstandardized model
    model <- lm(f_exp, data = df)

    coef <- summary(model)$coefficients["x", "Estimate"]
    se <- summary(model)$coefficients["x", "Std. Error"]
    lower_CI_coef <- confint(model)["x", "2.5 %"]
    upper_CI_coef <- confint(model)["x", "97.5 %"]
    R2 <- ((summary(model)$coefficients["x", "Estimate"])^2 * var(df$x)) /
      var(df$y)

    var_x <- var(df$x, na.rm = TRUE)
    var_y <- var(df$y, na.rm = TRUE)

    # Standardized model (safe scaling; skip when y or x has zero variance)
    coef_std <- se_std <- lower_CI_coef_std <- upper_CI_coef_std <- NA_real_

    sd_y0 <- isTRUE(sd(df$y, na.rm = TRUE) == 0)
    sd_x0 <- isTRUE(sd(df$x, na.rm = TRUE) == 0)

    if (!sd_y0) {
      df$y_z <- safe_scale(df$y)
      # keep only predictors with non-zero variance for standardized fit
      keep_preds <- predictors[vapply(
        predictors,
        function(p) sd(df[[p]], na.rm = TRUE) > 0,
        logical(1)
      )]
      if (length(keep_preds) > 0) {
        for (p in keep_preds) {
          df[[paste0(p, "_z")]] <- safe_scale(df[[p]])
        }
        rhs_scaled <- paste(sprintf("%s_z", keep_preds), collapse = " + ")
        f_std <- as.formula(sprintf("y_z ~ %s", rhs_scaled))
        # only fit if resulting model frame has rows
        if (nrow(model.frame(f_std, data = df)) > 0L) {
          model_std <- lm(f_std, data = df)
          sm_std <- summary(model_std)

          # x may have been dropped if sd_x0
          if (!sd_x0 && "x" %in% keep_preds) {
            coef_std <- sm_std$coefficients["x_z", "Estimate"]
            se_std <- sm_std$coefficients["x_z", "Std. Error"]
            lower_CI_coef_std <- confint(model_std)["x_z", "2.5 %"]
            upper_CI_coef_std <- confint(model_std)["x_z", "97.5 %"]
          }
        }
      }
    }

    # Effect sizes / correlations with guards
    x_type <- check_x_y_types(df)
    if (x_type == "numeric") {
      ok <- complete.cases(df$y, df$x)
      yv <- df$y[ok]
      xv <- df$x[ok]
      pearson <- lower_CI_pearson <- upper_CI_pearson <- NA_real_
      # only if we actually have variation
      if (length(yv) > 1 && sd(yv) > 0 && sd(xv) > 0) {
        # derive r from beta and variances
        pearson <- coef * sqrt(var_x) / sqrt(var_y)
        lower_CI_pearson <- lower_CI_coef * sqrt(var_x) / sqrt(var_y)
        upper_CI_pearson <- upper_CI_coef * sqrt(var_x) / sqrt(var_y)
      }
      parameter <- c(
        coef = coef,
        lower_CI_coef = lower_CI_coef,
        upper_CI_coef = upper_CI_coef,
        R2 = R2,
        coef_std = coef_std,
        lower_CI_coef_std = lower_CI_coef_std,
        upper_CI_coef_std = upper_CI_coef_std,
        pearson = pearson,
        lower_CI_pearson = lower_CI_pearson,
        upper_CI_pearson = upper_CI_pearson
      )
    } else {
      # categorical x
      d <- cohens_d_ci(df$x, df$y)
      cohen <- d$d
      lower_CI_cohen <- d$ci.lower
      upper_CI_cohen <- d$ci.upper
      parameter <- c(
        coef = coef,
        lower_CI_coef = lower_CI_coef,
        upper_CI_coef = upper_CI_coef,
        R2 = R2,
        coef_std = coef_std,
        lower_CI_coef_std = lower_CI_coef_std,
        upper_CI_coef_std = upper_CI_coef_std,
        cohen = cohen,
        lower_CI_cohen = lower_CI_cohen,
        upper_CI_cohen = upper_CI_cohen
      )
    }

    return(parameter)
  }

  # compute estimates for all samples
  estimates <- lapply(data, function(samples) {
    stopifnot(is.list(samples))
    lapply(samples, function(df) {
      # robust to failures in model_function
      tryCatch(
        model_function(df, formula = formula),
        error = function(e) {
          # try to guess x type to size the NA vector
          xt <- try(check_x_y_types(df), silent = TRUE)
          make_na_params(if (inherits(xt, "try-error")) "numeric" else xt)
        }
      )
    })
  })

  # long formatting
  estimates_to_long <- function(estimates, grid) {
    stopifnot(length(estimates) == nrow(grid))
    long <- do.call(
      rbind,
      lapply(seq_along(estimates), function(i) {
        scen <- grid[i, , drop = FALSE]
        samples <- estimates[[i]]
        do.call(
          rbind,
          lapply(names(samples), function(sname) {
            v <- samples[[sname]]
            data.frame(
              beta1 = scen$beta1,
              restrict = as.character(scen$restrict),
              sample = sname,
              parameter = names(v),
              estimate = as.numeric(v),
              row.names = NULL,
              check.names = FALSE
            )
          })
        )
      })
    )
    rownames(long) <- NULL
    long
  }

  df_long <- estimates_to_long(estimates, grid)

  return(df_long)
}

# helper: find which model-matrix column(s) belong to the term "x"
# (robust to binary factor/logical where the column is not literally named "x")
get_x_cols <- function(f_exp, df) {
  mm <- model.matrix(delete.response(terms(f_exp)), data = df)
  asg <- attr(mm, "assign")
  tl <- attr(terms(f_exp), "term.labels")

  # which term-label index corresponds to the literal term "x"
  x_term_idx <- which(tl == "x")
  if (length(x_term_idx) != 1) {
    return(character(0))
  }

  x_cols <- colnames(mm)[asg == x_term_idx]
  setdiff(x_cols, "(Intercept)")
}

true_quantities <- function(
  beta1,
  n,
  dgp_fun,
  formula,
  B = 10000,
  seed = 58294
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  f_in <- as.formula(formula)
  restrict_levels <- c("none", "weak", "moderate", "strong")

  truth_rows <- vector("list", length(beta1) * B)
  var_rows <- vector("list", length(beta1) * length(restrict_levels) * B * 2L)

  it_truth <- 0L
  it_var <- 0L

  for (b1 in beta1) {
    for (iter in seq_len(B)) {
      # unrestricted draw for true parameter values
      res_none <- dgp_fun(beta1 = b1, restrict = "none", n = n)
      df <- res_none$convenience_x

      f_exp <- formula(terms(f_in, data = df))
      mm_full <- model.matrix(delete.response(terms(f_exp)), data = df)
      mm_noi <- mm_full[, colnames(mm_full) != "(Intercept)", drop = FALSE]

      y <- df$y
      ok <- complete.cases(cbind(y, mm_noi))
      y_cc <- y[ok]
      X_cc <- mm_noi[ok, , drop = FALSE]

      var_x <- var(df$x, na.rm = TRUE)
      var_y <- var(df$y, na.rm = TRUE)
      x_type <- classify_x(df$x)

      m0 <- m1 <- sd0 <- sd1 <- NA_real_
      n0 <- n1 <- NA_integer_

      if (identical(x_type, "binary")) {
        g <- droplevels(as.factor(df$x))
        lev <- levels(g)
        y0 <- df$y[g == lev[1]]
        y1 <- df$y[g == lev[2]]

        m0 <- mean(y0, na.rm = TRUE)
        m1 <- mean(y1, na.rm = TRUE)
        sd0 <- sd(y0, na.rm = TRUE)
        sd1 <- sd(y1, na.rm = TRUE)
        n0 <- sum(is.finite(y0))
        n1 <- sum(is.finite(y1))
      }

      x_cols <- get_x_cols(f_exp, df)
      usable <- (length(y_cc) >= 3 && ncol(X_cc) >= 1 && length(x_cols) == 1)

      Sxx <- sxy <- syy <- NA
      if (usable) {
        Sxx <- cov(X_cc)
        sxy <- cov(X_cc, y_cc)
        syy <- var(y_cc)
      }

      it_truth <- it_truth + 1L
      truth_rows[[it_truth]] <- list(
        beta1 = b1,
        x_type = x_type,
        var_x = var_x,
        var_y = var_y,
        m0 = m0,
        m1 = m1,
        sd0 = sd0,
        sd1 = sd1,
        n0 = n0,
        n1 = n1,
        Sxx = Sxx,
        sxy = sxy,
        syy = syy,
        Xcols = colnames(mm_noi),
        x_col = if (length(x_cols) == 1) x_cols else NA_character_
      )

      # same iteration: all restriction levels for variance reduction
      for (rr in restrict_levels) {
        res_rr <- if (rr == "none") {
          res_none
        } else {
          dgp_fun(beta1 = b1, restrict = rr, n = n)
        }

        it_var <- it_var + 1L
        var_rows[[it_var]] <- data.frame(
          beta1 = b1,
          restrict = rr,
          sample = "convenience_x",
          target_var = "x",
          variance = .var_non_na(res_rr$convenience_x$x),
          stringsAsFactors = FALSE
        )

        it_var <- it_var + 1L
        var_rows[[it_var]] <- data.frame(
          beta1 = b1,
          restrict = rr,
          sample = "convenience_y",
          target_var = "y",
          variance = .var_non_na(res_rr$convenience_y$y),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  scalars_df <- do.call(
    rbind,
    lapply(truth_rows, function(z) {
      data.frame(
        beta1 = z$beta1,
        x_type = z$x_type,
        var_x = z$var_x,
        var_y = z$var_y,
        m0 = z$m0,
        m1 = z$m1,
        sd0 = z$sd0,
        sd1 = z$sd1,
        n0 = z$n0,
        n1 = z$n1,
        stringsAsFactors = FALSE
      )
    })
  )

  summary_df <- scalars_df %>%
    group_by(beta1, x_type) %>%
    summarise(
      var_x_mean = mean(var_x, na.rm = TRUE),
      var_y_mean = mean(var_y, na.rm = TRUE),
      m0_mean = mean(m0, na.rm = TRUE),
      m1_mean = mean(m1, na.rm = TRUE),
      sd0_mean = mean(sd0, na.rm = TRUE),
      sd1_mean = mean(sd1, na.rm = TRUE),
      n0 = mean(n0, na.rm = TRUE),
      n1 = mean(n1, na.rm = TRUE),
      .groups = "drop"
    )

  true_tbl <- summary_df
  true_tbl$coef <- true_tbl$beta1
  true_tbl$R2 <- ((true_tbl$coef)^2 * true_tbl$var_x_mean) / true_tbl$var_y_mean
  true_tbl$pearson <- true_tbl$coef *
    sqrt(true_tbl$var_x_mean) /
    sqrt(true_tbl$var_y_mean)
  true_tbl$coef_std <- NA_real_
  true_tbl$cohen <- NA_real_

  for (b1 in unique(true_tbl$beta1)) {
    idx <- which(vapply(
      truth_rows,
      function(z) isTRUE(z$beta1 == b1),
      logical(1)
    ))

    good <- idx[vapply(
      truth_rows[idx],
      function(z) {
        is.matrix(z$Sxx) &&
          is.numeric(z$syy) &&
          isTRUE(is.finite(z$syy)) &&
          is.character(z$x_col) &&
          length(z$x_col) == 1 &&
          !is.na(z$x_col)
      },
      logical(1)
    )]

    if (length(good) < 5) {
      next
    }

    cols0 <- truth_rows[[good[1]]]$Xcols
    if (
      !all(vapply(
        good,
        function(ii) identical(truth_rows[[ii]]$Xcols, cols0),
        logical(1)
      ))
    ) {
      next
    }

    xcol0 <- truth_rows[[good[1]]]$x_col
    if (
      !all(vapply(
        good,
        function(ii) identical(truth_rows[[ii]]$x_col, xcol0),
        logical(1)
      ))
    ) {
      next
    }

    Sxx_bar <- Reduce(`+`, lapply(good, function(ii) truth_rows[[ii]]$Sxx)) /
      length(good)
    sxy_bar <- Reduce(`+`, lapply(good, function(ii) truth_rows[[ii]]$sxy)) /
      length(good)
    syy_bar <- mean(
      vapply(good, function(ii) truth_rows[[ii]]$syy, numeric(1)),
      na.rm = TRUE
    )

    sdX <- sqrt(diag(Sxx_bar))
    sdY <- sqrt(syy_bar)
    if (any(!is.finite(sdX)) || any(sdX == 0) || !is.finite(sdY) || sdY == 0) {
      next
    }

    Rxx_bar <- Sxx_bar / (sdX %o% sdX)
    rxy_bar <- sxy_bar / (sdX * sdY)

    beta_std_all <- tryCatch(solve(Rxx_bar, rxy_bar), error = function(e) {
      rep(NA_real_, length(rxy_bar))
    })
    names(beta_std_all) <- cols0

    true_tbl$coef_std[true_tbl$beta1 == b1] <- unname(beta_std_all[xcol0])
  }

  sp2 <- ((true_tbl$n1 - 1) *
    true_tbl$sd1_mean^2 +
    (true_tbl$n0 - 1) * true_tbl$sd0_mean^2) /
    (true_tbl$n0 + true_tbl$n1 - 2)
  true_tbl$cohen <- (true_tbl$m1_mean - true_tbl$m0_mean) / sqrt(sp2)

  raw_var <- bind_rows(var_rows)

  var_summary <- raw_var %>%
    group_by(beta1, restrict, sample, target_var) %>%
    summarise(
      mean_variance = mean(variance, na.rm = TRUE),
      sd_variance = sd(variance, na.rm = TRUE),
      .groups = "drop"
    )

  base_none <- var_summary %>%
    filter(restrict == "none") %>%
    select(beta1, sample, target_var, none_variance = mean_variance)

  var_tbl <- var_summary %>%
    left_join(base_none, by = c("beta1", "sample", "target_var")) %>%
    mutate(
      variance_reduction = 1 - (mean_variance / none_variance),
      variance_reduction_pct = 100 * variance_reduction
    )

  list(true_values = true_tbl, variance_reduction = var_tbl)
}


# Monte Carlo wrapper (uses your one_simulation_run function)
run_simulations <- function(niter, beta1, restrict, n, formula, dgp_fun) {
  seed <- 58294

  truth <- true_quantities(
    beta1 = beta1,
    n = n,
    dgp_fun = dgp_fun,
    formula = formula,
    B = 1000,
    seed = seed
  )

  data_true <- truth$true_values
  data_varred <- truth$variance_reduction

  data_list <- vector("list", niter)
  for (it in seq_len(niter)) {
    set.seed(seed + it)
    res <- one_simulation_run(
      beta1 = beta1,
      restrict = restrict,
      n = n,
      formula = formula,
      dgp_fun = dgp_fun
    )
    data_list[[it]] <- res %>% mutate(iter = it)
  }

  data <- as.data.frame(bind_rows(data_list))
  rownames(data) <- NULL

  data_joined <- data %>%
    left_join(data_true, by = "beta1") %>%
    left_join(
      data_varred %>%
        select(
          beta1,
          restrict,
          sample,
          target_var,
          true_variance = mean_variance,
          true_variance_none = none_variance,
          variance_reduction,
          variance_reduction_pct
        ),
      by = c("beta1", "restrict", "sample")
    )

  data_relbias <- data_joined %>%
    mutate(
      true_value = case_when(
        parameter == "coef" ~ coef,
        parameter == "coef_std" ~ coef_std,
        parameter == "R2" ~ R2,
        parameter == "pearson" ~ pearson,
        parameter == "cohen" ~ cohen,
        TRUE ~ NA_real_
      ),
      abs_bias = estimate - true_value,
      rel_bias = (estimate - true_value) / true_value,
      rel_bias_pct = rel_bias * 100,
      abs_bias_pct = abs_bias * 100
    ) %>%
    select(
      beta1,
      restrict,
      sample,
      parameter,
      estimate,
      iter,
      x_type,
      true_value,
      abs_bias,
      abs_bias_pct,
      rel_bias,
      rel_bias_pct,
      target_var,
      true_variance,
      true_variance_none,
      variance_reduction,
      variance_reduction_pct
    )

  return(data_relbias)
}

# Simulation output
# Plot mean relative bias with MC SE error bars
simulation_results <- function(sim, restrict) {
  if (unique(sim$x_type) == "binary") {
    keep_params <- c(
      "coef",
      "coef_std",
      "R2",
      "cohen"
    )
  } else {
    keep_params <- c(
      "coef",
      "coef_std",
      "R2",
      "pearson"
    )
  }

  sim <- sim %>% filter(parameter %in% keep_params)
  digits <- 2

  # compute summary (only convenience samples, drop 'full')
  summary_mc <- as.data.frame(sim) %>%
    group_by(parameter, sample, beta1, restrict) %>%
    summarise(
      n_iter = n(),
      mean_rel = mean(rel_bias, na.rm = TRUE),
      mean_abs = mean(abs_bias, na.rm = TRUE),
      sd_rel = sd(rel_bias, na.rm = TRUE),
      sd_abs = sd(abs_bias, na.rm = TRUE),
      mc_se = sd_rel / sqrt(n_iter),
      mc_se_abs = sd_abs / sqrt(n_iter),
      pct_lo = quantile(rel_bias, probs = 0.025, na.rm = TRUE),
      pct_hi = quantile(rel_bias, probs = 0.975, na.rm = TRUE),
      mean_rel_pct = 100 * mean_rel,
      mean_abs_pct = 100 * mean_abs,
      mc_se_pct = 100 * mc_se,
      mc_se_abs_pct = 100 * mc_se_abs,
      RMSE = sqrt(mean((abs_bias)^2, na.rm = TRUE)),
      RMSE_pct = 100 * RMSE,
      pct_NA = 100 * mean(is.na(rel_bias)),
      .groups = "drop"
    ) %>%
    arrange(sample, parameter, beta1, restrict)

  # prepare data for plotting
  plot_df <- summary_mc %>%
    mutate(
      mean = mean_rel_pct,
      se = mc_se_pct,
      restrict = factor(
        restrict,
        levels = c("none", "weak", "moderate", "strong")
      )
    ) %>%
    select(parameter, sample, beta1, restrict, mean, se, RMSE_pct)

  params_present <- sort(unique(plot_df$parameter))

  param_labels <- c(
    coef = "Unstd. slope",
    coef_std = "Std. slope",
    R2 = "R\u00B2",
    pearson = "Pearson r",
    cohen = "Cohen's d"
  )[params_present]

  lab_fx <- labeller(
    parameter = as_labeller(param_labels),
    restrict = label_value # keep your restrict labels as-is
  )

  # Plot for convenience samples restricted on X
  px_other <- plot_df %>%
    filter(sample == "convenience_x") %>%
    ggplot(aes(
      x = factor(beta1),
      y = mean,
      color = restrict,
      group = restrict
    )) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    geom_line(
      position = position_dodge(width = 0.6),
      linewidth = 0.5,
      show.legend = FALSE
    ) +
    geom_pointrange(
      aes(ymin = mean - se, ymax = mean + se),
      position = position_dodge(width = 0.6),
      size = 0.5
    ) +
    facet_grid(~parameter, labeller = lab_fx) +
    scale_color_manual(
      values = c(
        none = "black",
        weak = "#74a9cf",
        moderate = "#0570b0",
        strong = "#034e7b"
      )
    ) +
    scale_y_continuous(
      limits = c(-100, 100),
      breaks = seq(-100, 100, by = 20)
    ) +
    labs(
      x = "True effect size",
      y = "Relative bias (%)",
      color = "Restriction strength",
      title = "Mean relative bias ± Monte Carlo SE by sample size and restriction",
      subtitle = "Convenience sample restricted on X",
      caption = "Error bars show Monte Carlo SE; y-axis clipped to ±100%,"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      axis.text.x = element_text(angle = 20, hjust = 1),
      strip.placement = "outside",
      strip.text.y.left = element_text(face = "bold"),
      strip.text.x = element_text(face = "bold"),
      plot.title.position = "plot",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10, margin = margin(b = 6)),
      plot.caption = element_text(size = 9, color = "gray30")
    )

  print(px_other)
  ggsave(
    "Mean relative bias for convenience sample restricted on X.png",
    width = 10,
    height = 10
  )

  py_other <- plot_df %>%
    filter(sample == "convenience_y") %>%
    ggplot(aes(
      x = factor(beta1),
      y = mean,
      color = restrict,
      group = restrict
    )) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    geom_line(
      position = position_dodge(width = 0.6),
      linewidth = 0.5,
      show.legend = FALSE
    ) +
    geom_pointrange(
      aes(ymin = mean - se, ymax = mean + se),
      position = position_dodge(width = 0.6),
      size = 0.5
    ) +
    facet_grid(~parameter, labeller = lab_fx) +
    scale_color_manual(
      values = c(
        none = "black",
        weak = "#74a9cf",
        moderate = "#0570b0",
        strong = "#034e7b"
      )
    ) +
    scale_y_continuous(
      limits = c(-100, 100),
      breaks = seq(-100, 100, by = 20)
    ) +
    labs(
      x = "True effect size",
      y = "Relative bias (%)",
      color = "Restriction strength",
      title = "Mean relative bias ± Monte Carlo SE by sample size and restriction",
      subtitle = "Convenience sample restricted on Y",
      caption = "Error bars show Monte Carlo SE; y-axis clipped to ±100%"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      axis.text.x = element_text(angle = 20, hjust = 1),
      strip.placement = "outside",
      strip.text.y.left = element_text(face = "bold"),
      strip.text.x = element_text(face = "bold"),
      plot.title.position = "plot",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10, margin = margin(b = 6)),
      plot.caption = element_text(size = 9, color = "gray30")
    )

  print(py_other)
  ggsave(
    "Mean relative bias for convenience sample restricted on Y.png",
    width = 10,
    height = 10
  )
}

# Plot power by restriction level, effect size, and convenience-sample type
power <- function(sim) {
  if (unique(sim$x_type) == "binary") {
    families <- c("coef", "coef_std", "cohen")
  } else {
    families <- c("coef", "coef_std", "pearson")
  }

  keys <- c("beta1", "restrict", "sample", "iter")

  # helper to compute power for one parameter family
  get_power_one_family <- function(sim, fam, null = 0) {
    ln <- paste0("lower_CI_", fam)
    un <- paste0("upper_CI_", fam)

    # lower & upper rows for this family
    lower <- sim %>%
      filter(parameter == ln) %>%
      select(all_of(keys), estimate_lower = estimate)

    upper <- sim %>%
      filter(parameter == un) %>%
      select(all_of(keys), estimate_upper = estimate)

    wide <- left_join(lower, upper, by = keys)

    wide %>%
      mutate(
        reject_est = (estimate_upper < null) | (estimate_lower > null)
      ) %>%
      group_by(beta1, restrict, sample) %>%
      summarise(
        power_est = mean(reject_est, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        parameter_family = fam
      )
  }
  # compute power for all families and bind rows
  power_tbl <- bind_rows(lapply(families, function(f) {
    get_power_one_family(sim, f)
  }))

  none <- power_tbl %>%
    filter(restrict == "none", sample == "convenience_x") %>%
    mutate(sample = "none") # label this baseline however you like

  # repeat for each restriction level
  none_expanded <- bind_rows(
    none %>% mutate(restrict = "weak"),
    none %>% mutate(restrict = "moderate"),
    none %>% mutate(restrict = "strong")
  )

  # create power table for restricted samples only
  power_tbl <- power_tbl %>% filter(restrict != "none")

  # create power table for no restriction (truth line)
  power_tbl_with_truth <- bind_rows(power_tbl, none_expanded) %>%
    mutate(
      restrict = factor(
        restrict,
        levels = c("none", "weak", "moderate", "strong")
      )
    )

  # plotting
  param_labels <- c(
    coef = "Unstd. slope",
    coef_std = "Std. slope",
    R2 = "R\u00B2",
    pearson = "Pearson r",
    cohen = "Cohen's d"
  )

  pow_plot <- ggplot() +
    geom_line(
      data = power_tbl_with_truth,
      aes(x = restrict, y = power_est, color = sample, group = sample),
      position = position_dodge(width = 0.3)
    ) +
    geom_point(
      data = power_tbl_with_truth,
      aes(x = restrict, y = power_est, color = sample),
      position = position_dodge(width = 0.3),
      size = 2
    ) +
    facet_grid(
      beta1 ~ parameter_family,
      labeller = labeller(parameter_family = param_labels)
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    labs(
      x = "Restriction level",
      y = "Power",
      color = "Convenience sample",
      title = "Power by restriction level, effect size, and convenience-sample type",
      subtitle = "Dashed black line shows reference (\"true\") power"
    ) +
    scale_color_manual(
      values = c(
        none = "black",
        convenience_x = "#542788", # light blue
        convenience_y = "#b35806" # dark blue
      )
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  print(pow_plot)
  ggsave("Power Plot.png", width = 10, height = 10)
}
