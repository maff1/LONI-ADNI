# GENERAL UTIL FUNCTIONS
# Last Observation to Carry Forward
na.locf <- function(x) {
  v <- !is.na(x)
  c(NA, x[v])[cumsum(v)+1]
}

LOCF <- function(x) {
  LOCF <- max(which(!is.na(x)))
  x[LOCF:length(x)] <- x[LOCF]
  return(x)
}

na.locf2 <- function(x) zoo::na.locf(x, na.rm = FALSE)
na.locf3 <- function(x) zoo::na.locf(x, na.rm = FALSE, fromLast=TRUE)

isnt_out_mad <- function(x, thres = 3, na.rm = TRUE) {
  abs(x - median(x, na.rm = na.rm)) <= thres * mad(x, na.rm = na.rm)
}

isnt_out_tukey <- function(x, k = 1.5, na.rm = TRUE) {
  quar <- quantile(x, probs = c(0.25, 0.75), na.rm = na.rm)
  iqr <- diff(quar)
  
  (quar[1] - k * iqr <= x) & (x <= quar[2] + k * iqr)
}

f_maha_outliers <- function(dat) {
  print(dim(dat))
  dat <- log2(dat + 1)
  md <- mahalanobis(dat, center = colMeans(dat), cov = cov(dat), tol = 1e-30)
  alpha <- .00001
  cutoff <- (qchisq(p = 1 - alpha, df = ncol(dat)))
  names_outliers_MH <- which(md > cutoff)
excluded_mh <- names_outliers_MH
data_clean_mh <- dat[-excluded_mh, ]
return(data_clean_mh)
}

`%nin%` <- Negate(`%in%`)

fold_cv=function(data,k){
  folds=cvTools::cvFolds(nrow(data),K=k)
  invisible(folds)
}

zscaleFunc <- function(x){(x-mean(x, na.rm = T))/sd(x, na.rm = T)}