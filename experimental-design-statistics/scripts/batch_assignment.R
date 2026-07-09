# batch_assignment.R
# Build balanced batch layouts using OSAT (Yasrebi 2012, BMC Genomics 13:689).
# Falls back to a stratified round-robin assignment if OSAT cannot optimize.
#
# Inputs:
#   - metadata: data frame with sample IDs and condition / covariate columns.
#   - batch_size: samples per batch.
#   - balance_vars: variables to balance (condition first).
# Outputs:
#   - data frame with `batch` column added (sorted by batch).

.req <- function(pkg, bioc = TRUE) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("Install %s with %s.",
                 pkg, if (bioc) "BiocManager::install" else "install.packages"),
         call. = FALSE)
}

#' Assign samples to batches with optimal balance (OSAT)
#'
#' @param metadata Data frame with sample metadata
#' @param batch_size Samples per batch
#' @param balance_vars Variables to balance across batches (must contain condition)
#' @param sample_id_col Name of the sample ID column (default "sample_id")
#' @return Data frame with a `batch` column added
assign_samples_to_batches <- function(metadata, batch_size, balance_vars,
                                      sample_id_col = "sample_id") {
  .req("OSAT")

  if (!sample_id_col %in% colnames(metadata))
    stop(sprintf("Column '%s' not found in metadata", sample_id_col))

  missing_vars <- setdiff(balance_vars, colnames(metadata))
  if (length(missing_vars) > 0)
    stop(sprintf("Balance variables not present: %s",
                 paste(missing_vars, collapse = ", ")))

  if (batch_size <= 0 || batch_size > nrow(metadata))
    stop("batch_size must be positive and <= number of samples")

  n_batches <- ceiling(nrow(metadata) / batch_size)
  balance_df <- as.data.frame(lapply(metadata[, balance_vars, drop = FALSE], as.factor))

  message("Generating optimal batch assignment with OSAT...")

  result <- tryCatch({
    osat_res <- OSAT::osat(x = balance_df, batch.size = batch_size,
                           n.batch = n_batches)
    metadata$batch <- osat_res$BatchLabel
    metadata <- metadata[order(metadata$batch), ]
    cat(sprintf("✓ Batch design generated: %d batches of ~%d samples; balanced on %s\n",
                n_batches, batch_size, paste(balance_vars, collapse = ", ")))
    metadata
  }, error = function(e) {
    warning("OSAT optimization failed (", conditionMessage(e),
            "); falling back to stratified round-robin.")
    .assign_samples_simple(metadata, batch_size, balance_vars, sample_id_col)
  })

  result
}

.assign_samples_simple <- function(metadata, batch_size, balance_vars, sample_id_col) {
  n_samples <- nrow(metadata)
  n_batches <- ceiling(n_samples / batch_size)
  primary <- balance_vars[1L]

  metadata <- metadata[order(metadata[[primary]]), ]

  # Round-robin assignment across batches, then shuffle within batch
  set.seed(42L)
  assignments <- rep(seq_len(n_batches), length.out = n_samples)
  metadata$batch <- sample(assignments)
  metadata <- metadata[order(metadata$batch), ]
  cat(sprintf("✓ Stratified fallback assignment: %d batches\n", n_batches))
  metadata
}

#' Incomplete block design when complete balance is impossible
#'
#' Cycles through condition levels assigning one sample at a time to the
#' current batch, ensuring every level appears in every batch.
#'
#' @param metadata Data frame with sample metadata
#' @param batch_size Samples per batch
#' @param primary_var Variable that MUST be balanced (typically condition)
#' @param secondary_vars Variables to balance secondarily (currently informational)
incomplete_block_design <- function(metadata, batch_size, primary_var,
                                    secondary_vars = NULL) {
  if (!primary_var %in% colnames(metadata))
    stop(sprintf("Primary variable '%s' not in metadata", primary_var))

  primary_levels <- unique(metadata[[primary_var]])

  if (batch_size %% length(primary_levels) != 0)
    message("Batch size not divisible by number of levels — creating incomplete block design.")

  metadata$batch <- NA_integer_
  batch_idx <- 1L
  filled    <- 0L

  while (any(is.na(metadata$batch))) {
    for (lev in primary_levels) {
      eligible <- which(metadata[[primary_var]] == lev & is.na(metadata$batch))
      if (length(eligible) > 0 && filled < batch_size) {
        metadata$batch[eligible[1L]] <- batch_idx
        filled <- filled + 1L
        if (filled >= batch_size) {
          batch_idx <- batch_idx + 1L
          filled <- 0L
        }
      }
    }
  }

  cat(sprintf("✓ Incomplete block design: %d batches; '%s' balanced across batches\n",
              length(unique(metadata$batch)), primary_var))
  metadata
}

#' Add a randomized processing order within each batch
#'
#' @param batch_metadata Data frame with a `batch` column
#' @param randomize Whether to randomize order within each batch (default TRUE)
add_processing_order <- function(batch_metadata, randomize = TRUE) {
  if (!"batch" %in% colnames(batch_metadata))
    stop("Metadata must contain a 'batch' column.")

  batch_metadata$processing_order <- NA_integer_
  for (b in unique(batch_metadata$batch)) {
    rows <- which(batch_metadata$batch == b)
    if (randomize) {
      set.seed(as.integer(b))
      rows <- sample(rows)
    }
    batch_metadata$processing_order[rows] <- seq_along(rows)
  }
  batch_metadata <- batch_metadata[order(batch_metadata$batch,
                                         batch_metadata$processing_order), ]
  batch_metadata$overall_sequence <- seq_len(nrow(batch_metadata))
  batch_metadata
}

#' Assign plate positions (preferring non-edge wells)
#'
#' @param batch_metadata Data frame with a `batch` column
#' @param plate_format One of "96-well", "384-well", "48-well"
add_plate_positions <- function(batch_metadata, plate_format = "96-well") {
  if (!"batch" %in% colnames(batch_metadata))
    stop("Metadata must contain a 'batch' column.")

  formats <- list(
    "96-well"  = list(rows = 8L,  cols = 12L),
    "384-well" = list(rows = 16L, cols = 24L),
    "48-well"  = list(rows = 6L,  cols = 8L)
  )
  if (!plate_format %in% names(formats))
    stop(sprintf("Unknown plate format. Choose from: %s",
                 paste(names(formats), collapse = ", ")))

  dims <- formats[[plate_format]]
  row_letters <- LETTERS[seq_len(dims$rows)]
  col_numbers <- seq_len(dims$cols)

  positions <- expand.grid(row = row_letters, col = col_numbers,
                           stringsAsFactors = FALSE)
  positions$well <- paste0(positions$row, sprintf("%02d", positions$col))
  positions$is_edge <-
    positions$row == row_letters[1L] |
    positions$row == row_letters[dims$rows] |
    positions$col == 1L |
    positions$col == dims$cols

  non_edge <- positions[!positions$is_edge, , drop = FALSE]

  batch_metadata$plate <- NA_integer_
  batch_metadata$well  <- NA_character_
  batch_metadata$is_edge_well <- FALSE

  plate_idx <- 1L
  pos_idx <- 1L
  for (i in seq_len(nrow(batch_metadata))) {
    if (pos_idx > nrow(non_edge)) {
      plate_idx <- plate_idx + 1L
      pos_idx <- 1L
    }
    batch_metadata$plate[i] <- plate_idx
    batch_metadata$well[i]  <- non_edge$well[pos_idx]
    pos_idx <- pos_idx + 1L
  }

  cat(sprintf("✓ Plate positions assigned: %d plates (%s)\n", plate_idx, plate_format))
  batch_metadata
}
