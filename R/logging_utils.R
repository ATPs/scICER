#' Timestamped logger for scICER runtime messages
#' @keywords internal
scice_message <- function(..., appendLF = TRUE, domain = NULL) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste0(..., collapse = "")
  base::message(sprintf("[%s] %s", ts, msg), appendLF = appendLF, domain = domain)
}
