.kblock <- function(kf1i, kf0i, kf1j, kf0j, K, n) {
  1 / (K * (n - 1)) *
    ((K - 1) *
       (kf1i^2 - kf1i * kf0i + kf1j^2 - kf1j * kf0j) -
       K * (K - 1) * kf1i * kf1j)
}

.gate_labels <- function(tau, ngates, context = "score") {
  if (length(unique(tau)) == 1L) {
    warning(
      sprintf(
        paste0(
          "%s is constant. Gates are assigned at random, so the test ",
          "remains valid but has no power from the supplied ranking."
        ),
        context
      ),
      call. = FALSE
    )
  } else if (anyDuplicated(tau)) {
    warning(
      sprintf(
        paste0(
          "%s contains ties. Assignment among tied observations is ",
          "randomized to avoid dependence on row order."
        ),
        context
      ),
      call. = FALSE
    )
  }

  dplyr::ntile(rank(tau, ties.method = "random"), ngates)
}
