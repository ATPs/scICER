#' Timestamped logger for scICER runtime messages
#'
#' @description
#' Prepends a wall-clock timestamp to messages for long-running workflow logs.
#'
#' @details
#' This wrapper delegates to \code{base::message()} while formatting all messages
#' as \code{[YYYY-mm-dd HH:MM:SS] <text>}. It is used throughout verbose paths to
#' make multi-stage runtime diagnostics easier to follow.
#'
#' @keywords internal
scice_message <- function(..., appendLF = TRUE, domain = NULL) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste0(..., collapse = "")
  base::message(sprintf("[%s] %s", ts, msg), appendLF = appendLF, domain = domain)
}
