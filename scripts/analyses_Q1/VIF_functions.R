
#### FUNCTIONS

# Variance Inflation Factor (VIF) functions to check for multi-collinearity among model predictors.  
# Citation: Mixed effects models and extensions in ecology with R. (2009).
#           Zuur, AF, Ieno, EN, Walker, N, Saveliev, AA, and Smith, GM.
#           Springer. http://www.highstat.com/book2.htm
corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1,dataz)
  lm_mod  <- lm(form,dataz)
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}

# Handy Variogram function to check for spatial residual autocorrelation.
# Citation: Mixed effects models and extensions in ecology with R. (2009).
#           Zuur, AF, Ieno, EN, Walker, N, Saveliev, AA, and Smith, GM.
#           Springer. http://www.highstat.com/book2.htm
MyVariogram <- function(x,y,z, MyDistance) {
  mydata      <- data.frame(z, x, y)
  coordinates(mydata)    <- c("x", "y")  
  Var <- variogram(z ~ 1, mydata, cutoff = MyDistance)
  data.frame(Var$np, Var$dist, Var$gamma)
}
