# got it from here: https://raw.githubusercontent.com/Arcadia-Science/prehgt/ba092b0e01a20688de56292883d7e35be0e2db08/bin/blastp_to_hgt_candidates_kingdom.R
gini <- function(x) {
  # calculate gini coefficient
  x_sorted <- sort(x)
  n <- length(x_sorted)
  return (1 - 2 * (sum((1:n) * x_sorted) / sum(x_sorted) - (n + 1) / 2) / n)
}

dmnd_cols <- c(
  "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send",
  "evalue", "bitscore", "qcovhsp", "qlen", "slen"
)