# Copied verbatim from
# buffy_methylation_architecture/scripts/cfdna_nucleosome_periodicity/_peak_detection_utils.r
# (already verified there against synthetic planted-period data and real
# mod_frags.txt rows) so mono-/di-nucleosome peak detection on the 5-base
# fragment-length histogram in this project uses the identical, already-
# tested method rather than a reimplementation that could silently drift.
# Copied rather than sourced across repos so this project stays self-
# contained and runnable without a checkout of buffy_methylation_architecture
# on the same machine.

# Centered moving-average smooth.
smooth_signal <- function(x, window) {
  n <- length(x)
  half <- window %/% 2L
  out <- numeric(n)
  for (i in seq_len(n)) {
    lo <- max(1L, i - half)
    hi <- min(n, i + half)
    out[i] <- mean(x[lo:hi])
  }
  out
}

# Strict > (not >=) against both neighbors - rules out flagging every bin
# in a flat plateau as its own peak without a separate collapsing step.
find_local_maxima <- function(x, min_prominence) {
  n <- length(x)
  is_peak <- logical(n)
  if (n < 3) return(is_peak)  # 2:(n-1) misbehaves (descends) if n < 2
  for (i in 2:(n - 1)) {
    is_peak[i] <- x[i] > x[i - 1] && x[i] > x[i + 1] && x[i] > min_prominence
  }
  is_peak
}

# Mean peak-to-peak spacing among the most prominent peaks. See the
# original in buffy_methylation_architecture for the full rationale on why
# selection must be magnitude-first-then-top-N-then-position, not
# position-first or top-N-without-a-magnitude-floor.
spacing_from_top_peaks <- function(peaks_dt, position_col, magnitude_col,
                                    min_spacing_prominence = -Inf, n_spacing_peaks = 5) {
  if (nrow(peaks_dt) < 2) return(list(spacings = numeric(0), mean_spacing = NA_real_, peaks_used = peaks_dt[0]))
  eligible <- peaks_dt[get(magnitude_col) > min_spacing_prominence]
  ordered_by_magnitude <- eligible[order(-get(magnitude_col))]
  top_peaks <- head(ordered_by_magnitude, n_spacing_peaks)
  top_peaks_by_position <- top_peaks[order(get(position_col))]
  if (nrow(top_peaks_by_position) < 2) {
    return(list(spacings = numeric(0), mean_spacing = NA_real_, peaks_used = top_peaks_by_position))
  }
  spacings <- diff(top_peaks_by_position[[position_col]])
  list(spacings = spacings, mean_spacing = mean(spacings), peaks_used = top_peaks_by_position)
}
