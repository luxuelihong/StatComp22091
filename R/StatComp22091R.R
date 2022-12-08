#' @title  Calculate Delta using R
#' @description Calculate Delta from sample split paper for a single observation
#' @param data_i The i^th observation (Y_i, X_i, Z_i)
#' @param big_delta_symbolic Symbolic expression for Delta
#' @param params_first names of first stage parameters (theta) in a list
#' @param params_second names of second stage parameters (beta) in a list
#' @param  param_vals estimated values of first and second stage parameters
#' @return a matrix named big_delta of size row length(params_second) col length(params_first))
#' @export
delta_i = function(data_i, big_delta_symbolic, params_first,  params_second, param_vals 
){
  
  # Iterate through the symbolic expression for Delta and 
  # evaluate each element numerically.
  big_delta = matrix(NA, nrow = length(params_second), ncol = length(params_first))
  envir = as.list(c(data_i, param_vals))
  for (r in 1:nrow(big_delta)) {
    for (c in 1:ncol(big_delta)) {
      big_delta[r,c] = eval(
        big_delta_symbolic[r,c][[1]],
        envir = envir)
    }
  }
  
  return(big_delta)
}



