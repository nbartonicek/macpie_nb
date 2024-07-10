data_dir<-"inst/extdata/filtered_matrix/"
label_file<-"inst/extdata/metadata/metadata.csv"
gene_column = 2
cell_column = 1
unique_features = TRUE
strip_suffix = FALSE


#' data_dir <- 'path/to/data/directory'
#' list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
#' expression_matrix <- Readmacseq(data_dir = data_dir)
#' macpie_object = CreateMacPie(counts = expression_matrix)
#' }
#'
read_macseq <- function(
    data_dir,
    label_file,
    gene_column = 2,
    cell_column = 1,
    unique_features = TRUE,
    strip_suffix = FALSE
) {

  full_data <- list()
  has_dt <- requireNamespace("data.table", quietly = TRUE) &&
    requireNamespace("R.utils", quietly = TRUE)
  run <- data_dir
  metadata <- label_dir

  if (!dir.exists(paths = run)) {
    stop("Macseq data directory provided does not exist")
  }
  barcode_loc <- file.path(run, "barcodes.tsv.gz")
  features_loc <- file.path(run, "features.tsv.gz")
  matrix_loc <- file.path(run, "matrix.mtx.gz")

  if (!file.exists(label_file)) {
    stop("Metadata file provided does not exist")
  }
  label_file <- list.files(metadata, pattern = "\\.csv$", full.names = TRUE)

  if (!any(file.exists(c(barcode_loc,features_loc,matrix_loc)))) {
    stop(paste(
      "Missing one or more of count matrix files. Expecting",
      basename(path = barcode_loc),
      basename(path = features_loc),
      basename(path = matrix_loc))
    )
  }

  if (!file.exists(label_file)) {
    stop("Metadata file missing. Expecting ", basename(path = label_file))
  }
  warning("Make sure you only have 1 metadata file in the folder.")

  data <- readMM(file = matrix_loc)
  cell_barcodes <- as.data.frame(data.table::fread(barcode_loc, header = FALSE))


  if (ncol(x = cell_barcodes) > 1) {
    cell_names <- cell_barcodes[, cell_column]
  } else {
    cell_names <- readLines(con = barcode_loc)
  }
  if (all(grepl(pattern = "\\-1$", x = cell_names)) && strip_suffix) {
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
    feature_names <- as.data.frame(
      data.table::fread(ifelse(pre_ver_3, gene_loc, features_loc), header = FALSE)
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
    na_features <- which(is.na(feature_names[, gene_column]))
    replacement_column <- ifelse(gene_column == 2, yes = 1, no = 2)
    feature_names[na_features, gene_column] <-
      feature_names[na_features, replacement_column]
  }
  if (unique_features) {
    fcols <- ncol(x = feature_names)
    if (fcols < gene_column) {
      stop(
        paste0(
          "gene_column was set to ",
          gene_column,
          " but feature.tsv.gz (or genes.tsv) only has ",
          fcols,
          " columns. Try setting the gene_column argument to a value <= to ",
          fcols, ".")
      )
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
  for (j in seq_along(full_data[[1]])) {
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
