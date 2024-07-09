
#' Load in data from Mac-seq experiments
#'
#' This functions is a derivative from Read10X function from the Seurat package.
#' Enables easy loading of sparse data matrices provided by 10X genomics.
#' @author Mark Li, \email{justdont}
#'
#' @param data_dir Directory containing the matrix.mtx, genes.tsv
#' (or features.tsv), and barcodes.tsv files provided by nextflow
#' pipeline using Starsolo.
#' Functionality for reading several data directories was removed
#' to allow reading in the metadata label matching to the Macseq
#' plate with minimum confusion.
#' Next version will consider integration of other modalities.
#'
#' @param label_dir Directory containing the matching metadata file in csv format.
#' @param gene_column Specify which column of genes.tsv or features.tsv
#' to use for gene names; default is 2
#' @param cell_column Specify which column of barcodes.tsv
#' to use for cell names; default is 1
#' @param unique.features Make feature names unique (default TRUE)
#' @param strip.suffix Remove trailing "-1" if present in all cell barcodes.
#'
#'
#' @return If features.csv indicates the data has multiple data types, a list
#'   containing a sparse matrix of the data from each type will be returned.
#'   Otherwise a sparse matrix containing the expression data will be returned.
#'
#' @importFrom Matrix readMM
#' @importFrom utils read.delim
#'
#' @export
#' @concept preprocessing
#'
#' @examples
#' \dontrun{
#' #
#' data_dir <- 'path/to/data/directory'
#' list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
#' expression_matrix <- Readmacseq(data_dir = data_dir)
#' macpie_object = CreateMacPie(counts = expression_matrix)
#' }
#'
readMacSeq <- function(
  data_dir,
  label_dir,
  gene_column = 2,
  cell_column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
) {

  full_data <- list()
  has_dt <- requireNamespace("data.table", quietly = TRUE)
    && requireNamespace("R.utils", quietly = TRUE)
  run <- data_dir
  metadata <- label_dir

  if (!dir.exists(paths = run)) {
    stop("Macseq data directory provided does not exist")
  }
  barcode_loc <- file.path(run, "barcodes.tsv")
  features_loc <- file.path(run, "features.tsv.gz")
  matrix_loc <- file.path(run, "matrix.mtx")

  if (!dir.exists(paths = metadata)) {
    stop("Metadata directory provided does not exist")
  }
  label_loc <- list.files(metadata, pattern = "\\.csv$", full.names = TRUE)

  if (!file.exists(barcode_loc)) {
    stop("Barcode file missing. Expecting ",
         basename(path = barcode_loc))
  }
  if (!file.exists(features_loc)) {
    stop("Features file missing. Expecting ",
         basename(path = features_loc))
  }
  if (!file.exists(matrix_loc)) {
    stop("Expression matrix file missing. Expecting ",
         basename(path = matrix_loc))
  }
  if (!file.exists(label_loc)) {
    stop("Metadata file missing. Expecting ", basename(path = label_loc))
  }
  warning("Make sure you only have 1 metadata file in the folder.")

  data <- readMM(file = matrix_loc)
  if (has_dt) {
    cell_barcodes <- as.data.frame(
      data.table::fread(barcode_loc, header = FALSE)
    )
  } else {
    cell_barcodes <- read.table(file = barcode_loc,
      header = FALSE, sep = "\t", row.names = NULL)
  }

  if (ncol(x = cell_barcodes) > 1) {
    cell_names <- cell_barcodes[, cell_column]
  } else {
    cell_names <- readLines(con = barcode_loc)
  }
  if (all(grepl(pattern = "\\-1$", x = cell_names)) && strip.suffix) {
    cell_names <- as.vector(x = as.character(x = sapply(
      X = cell_names,
      FUN = Seurat::ExtractField,
      field = 1,
      delim = "-"
    )))
  }
  if (is.null(x = names(x = data_dir))) {
    if (length(x = data_dir) < 2) {
      colnames(x = data) <- cell_names
    } else {
      colnames(x = data) <- paste0(1, "_", cell_names)
    }
  } else {
    colnames(x = data) <- paste0(names(x = data_dir)[1], "_", cell_names)
  }

  if (has_dt) {
    feature_names <- as.data.frame(data.table::fread(
      ifelse(pre_ver_3, gene_loc, features_loc), header = FALSE)
    )
  } else {
    feature_names <- read.delim(
      file = features_loc,
      header = FALSE,
      stringsAsFactors = FALSE
    )
  }

  if (any(is.na(x = feature_names[, gene_column]))) {
    warning(
      "Some features names are NA. Replacing NA names with
      ID from the opposite column requested.",
      call. = FALSE,
      immediate. = TRUE
    )
    na.features <- which(x = is.na(x = feature_names[, gene_column]))
    replacement_column <- ifelse(test = gene_column == 2, yes = 1, no = 2)
    feature_names[na.features, gene_column] <-
      feature_names[na.features, replacement_column]
  }
  if (unique.features) {
    fcols <- ncol(x = feature_names)
    if (fcols < gene_column) {
      stop(paste0("gene_column was set to ", gene_column,
        " but feature.tsv.gz (or genes.tsv) only has ", fcols, " columns.",
        " Try setting the gene_column argument to a value <= to ", fcols, "."))
    }
    rownames(x = data) <- make.unique(names = feature_names[, gene_column])
  }

  # In cell ranger 3.0, a third column specifying the type of data was added
  # and we will return each type of data as a separate matrix.
  # We decided to leave this in for now as "Gene Expression" column
  # still exists in the features.tsv file from StarSolo output.
  if (ncol(x = feature_names) > 2) {
    data_types <- factor(x = feature_names$V3)
    lvls <- levels(x = data_types)
    if (length(x = lvls) > 1 && length(x = full_data) == 0) {
      message("10X data contains more than one type
              and is being returned as a list containing matrices of each type.")
    }
    expr_name <- "Gene Expression"
    if (expr_name %in% lvls) { # Return Gene Expression first
      lvls <- c(expr_name, lvls[-which(x = lvls == expr_name)])
    }
    data <- lapply(
      X = lvls,
      FUN = function(l) {
        return(data[data_types == l, , drop = FALSE])
      }
    )
    names(x = data) <- lvls
  } else {
    data <- list(data)
  }

  full_data[[length(x = full_data) + 1]] <- data

  # Combine all the data from different directories
  # into one matrix. This assumes that all data
  # directories  have the same features files.

  list_of_data <- list()
  for (j in seq_along(x = full_data[[1]])) {
    list_of_data[[j]] <- do.call(cbind, lapply(X = full_data, FUN = `[[`, j))
  }
  names(x = list_of_data) <- names(x = full_data[[1]])
  # If multiple features, will return a list, otherwise
  # a matrix.
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  } else {
    return(list_of_data)
  }
}
