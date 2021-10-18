## -------------------------------------------------------- #
## Author: Reto Buergin, rbuergin@gmx.ch
##
## Contents:
## R-Functions that extend the 'regsubset' function from 
## the R-package 'leaps'
## - regsubsets.AIC: Extracts the AIC of the best models
##      per number of predictors.
## - plot.regsubsets.selected: Visualize the best models
##      per number of predictors using the 'image' function.
##      Includes also the model without predictors.
## - plot.regsubsets.err: Visualize the estimated prediction
##      error (electively AIC, BIC, ...) from the best models
##      per number of predictors in a scatter plot.
## - regsubset.getbestmodel: Extract the best model from
##      a 'regsubset' object.
##
## Modifications:
## 2021-10-18: Add comments
## 2021-10-13: Create file
## -------------------------------------------------------- #

regsubsets.AIC <- function(object, data, yname, nullModel = FALSE) {
  require(leaps)
  stopifnot(inherits(object, "regsubsets"))
  rval <- rep(NA, nrow(summary(object)$which))
  for (i in 1:nrow(summary(object)$which)) {
    xnames <- object$xnames[summary(object)$which[i, ]]
    if ("(Intercept)" %in% xnames) {
      xnames <- setdiff(xnames, "(Intercept)")
    } else {
      xnames <- c("-1", xnames)
    }
    form <- as.formula(paste0(yname, " ~ ", paste0(xnames, collapse = " + ")))
    rval[i] <- as.numeric(AIC(lm(
      formula = form, 
      data = data)))
  }
  if (nullModel) {
    rval <- c(
      as.numeric(AIC(lm(
        formula = as.formula(paste0(yname, " ~ 1")), 
        data = data))),
      rval)
  }
  return(rval)
}

## -------------------------------------------------------- #

plot.regsubsets.selected <- function(
  x, col = c("white", "black"), ...) {
  require(leaps)
  stopifnot(inherits(x, "regsubsets"))
  which <- summary(x)$which
  if ("(Intercept)" %in% colnames(which)) {
    which <- rbind(rep(FALSE, ncol(which)), which)
    which[1, "(Intercept)"] <- TRUE
  }
  which <- t(which)
  colnames(which) <- 0:(ncol(which) - 1)
  image(
    z = which,
    xaxt = "n", yaxt = "n",
    x = 1:ncol(which),
    y = 1:nrow(which),
    xlab = "", ylab = "Number of Predictors",
    col = col)
  axis(side = 1, at = 1:nrow(which), labels = rownames(which), las = 2)
  axis(side = 2, at = 1:ncol(which), labels = colnames(which), las = 2)
  invisible()
}

## -------------------------------------------------------- #

plot.regsubsets.err <- function(
  x, data, yname, criteria = c("aic", "bic", "rss", "rsq", "adjr2"), ...) {
  require(leaps)
  stopifnot(inherits(x, "regsubsets"))
  criteria <- match.arg(criteria)
  err <- switch(
    criteria,
    aic = regsubsets.AIC(x, data, yname),
    bic = summary(x)$bic,
    rss = summary(x)$rss,
    rsq = summary(x)$rsq,
    adjr2 = summary(x)$adjr2)
  nullModel <- lm(
    formula = as.formula(paste0(yname, " ~ 1")), 
    data = data)
  err0 <- switch(
    criteria,
    aic = as.numeric(AIC(nullModel)),
    bic = as.numeric(BIC(nullModel)),
    rss = as.numeric(sum(resid(nullModel)^2)),
    rsq = as.numeric(summary(nullModel)$r.squared),
    adjr2 = as.numeric(summary(nullModel)$adj.r.squared))
  err <- c(err0, err)
  plot(
    x = 0:(length(err) - 1), 
    y = err,
    type = "b",
    xlab = "Number of Predictors",
    ylab = criteria)
  abline(v = which.min(err) - 1, col = "red", lty = 2)
  invisible()
}

## -------------------------------------------------------- #

regsubset.getbestmodel <- function(
  object, data, yname, criteria = c("aic", "bic", "rss", "rsq", "adjr2")) {
  require(leaps)
  stopifnot(inherits(object, "regsubsets"))
  criteria <- match.arg(criteria)
  err <- switch(
    criteria,
    aic = regsubsets.AIC(object, data, yname),
    bic = summary(object)$bic,
    rss = summary(object)$rss,
    rsq = -summary(object)$rsq,
    adjr2 = -summary(object)$adjr2)
  selected <- colnames(summary(object)$which)[
    summary(object)$which[which.min(err), ]]
  if ("(Intercept)" %in% selected) {
    selected <- setdiff(selected, "(Intercept)")
  } else {
    selected <- c("-1", selected)
  }
  rval <- lm(
    formula = as.formula(paste0(yname, " ~ ", paste(selected, collapse = " + "))),
    data = data)
  return(rval)
}

## EOF ---------------------------------------------------- #