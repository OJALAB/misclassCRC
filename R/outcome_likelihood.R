#' Subtract Probabilities on the Log Scale
#' 
#' @noRd
log_sub_crc <- function(log_x, log_y) {

  difference <- log_y - log_x

  both_zero <- is.infinite(log_x) &
    log_x < 0 &
    is.infinite(log_y) &
    log_y < 0

  difference[both_zero] <- -Inf

  log_x + log(-expm1(difference))

}