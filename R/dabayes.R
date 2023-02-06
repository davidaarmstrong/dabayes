#' Convert JAGS output to Tibble
#'
#' Converts samples created by `rjags::coda.samples()` to a data frame similar to that
#' produced by the `tidybayes` package.
#'
#' @param coda_object An object of class `mcmc.list` produced by `rjags::coda.samples()`.
#' @param ... Other arguments to pass down, currently unimplemented.
#'
#' @importFrom dplyr tibble bind_rows
#' @export
coda_to_df <- function(coda_object, ...){
  mcpar <- attr(coda_object[[1]], "mcpar")
  if(length(mcpar) == 3){
    thin <- mcpar[3]
  }else{
    thiin <- 1
  }
  sq <- seq(mcpar[1], mcpar[2], by=thin)
  dats <- lapply(coda_object, function(x){
    nrow <- nrow(x)
    ncol <- ncol(x)
    tibble(.iteration = rep(sq, ncol),
           .draw = rep(1:nrow, ncol),
           param = rep(colnames(x), each=nrow),
           param_num = rep(1:ncol, each=nrow),
           val = c(x))
  })
  bind_rows(dats, .id=".chain")
}

#' Pairwise Comparisons of MCMC Output
#'
#' Generates pairwise comparisons of MCMC output that can then be either
#' summarised or plotted. Creates differences by calculating the difference of
#' lhs_name - rhs_name.
#'
#' @param samps A matrix of MCMC output or an `mcmc.list`.
#' @param hdi_prob The probability to be covered by the HDI, defaults to 0.95.
#' @param ... Other arguments to be passed down, currently unimplemented.
#'
#' @importFrom dplyr bind_rows `%>%` rowwise mutate filter
#' @importFrom runjags combine.mcmc
#' @importFrom coda HPDinterval as.mcmc
#' @importFrom rlang .data
#' @export
compare_params <- function(samps, hdi_prob = .95, ...){
  if(inherits(samps, "mcmc.list") | inherits(samps, "mcmc")){
    samps <- combine.mcmc(samps)
  }
  if(is.null(samps))stop("samps should be a matrix.\n")
  if(length(dim(samps)) != 2)stop("samps should be a matrix.\n")
  nms <- colnames(samps)
  if(is.null(nms)){
    nms <- paste0("x", 1:ncol(samps))
  }
  out <- expand.grid(rhs_num = 1:ncol(samps),
                     lhs_num = 1:ncol(samps)) %>%
    rowwise() %>%
    mutate(p_above_0 = mean(samps[,.data$lhs_num] - samps[, .data$rhs_num] > 0),
           diff_hdi_lwr = HPDinterval(as.mcmc(samps[, .data$lhs_num] - samps[, .data$rhs_num]), prob=hdi_prob)[1],
           diff_hdi_upr = HPDinterval(as.mcmc(samps[, .data$lhs_num] - samps[, .data$rhs_num]), prob=hdi_prob)[2],
           rhs_param = mean(samps[,.data$rhs_num]),
           lhs_param = mean(samps[,.data$lhs_num]),
    ) %>%
    mutate(rhs_name = nms[.data$rhs_num],
           lhs_name = nms[.data$lhs_num]) %>%
    filter(.data$rhs_num != .data$lhs_num)
  out
}


