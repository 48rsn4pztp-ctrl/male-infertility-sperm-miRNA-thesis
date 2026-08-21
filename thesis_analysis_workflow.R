#!/usr/bin/env Rscript

###############################################################################
# Reproducible analysis workflow for the MSc sperm-miRNA thesis
#
# This public script reproduces the computational sequence described in the
# thesis without distributing participant-level counts or metadata. Authorised
# users must place the two required input files in a local data directory:
#
#   miRNA_expression.csv  first column = miRNA identifier; remaining columns = samples
#   demographics.txt      columns = Subject_ID, Age, BMI, Ethnicity, Status
#
# Status must contain Fertile, Borderline and Infertile. Only Fertile and
# Infertile samples are used for supervised development. Borderline samples are
# projected after the final model is fixed.
#
# Example:
#   THESIS_DATA_DIR=/path/to/private/data Rscript thesis_analysis_workflow.R
#
# Optional switches (TRUE/FALSE):
#   RUN_NESTED_COMPARISON       default TRUE
#   RUN_PARAMETER_OPTIMISATION  default TRUE
#   RUN_PANEL_REDUCTION         default TRUE
#   RUN_RANDOM_FOREST_CHECK     default TRUE
#   RUN_FUNCTIONAL_ANALYSIS     default FALSE (requires database access)
#
# This is a development workflow. Cross-validation estimates are internal and
# conditional on the candidate set selected before resampling; they are not an
# independent clinical validation.
###############################################################################

options(stringsAsFactors = FALSE, warn = 1)

###############################################################################
# 1. Configuration
###############################################################################

env_flag <- function(name, default = TRUE) {
  value <- Sys.getenv(name, unset = if (default) "TRUE" else "FALSE")
  tolower(value) %in% c("1", "true", "t", "yes", "y")
}

env_integer <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (is.na(value) || value < 1L) default else value
}

data_dir <- Sys.getenv("THESIS_DATA_DIR", unset = "data")
output_dir <- Sys.getenv("THESIS_OUTPUT_DIR", unset = "outputs")
counts_file <- Sys.getenv(
  "THESIS_COUNTS_FILE",
  unset = file.path(data_dir, "miRNA_expression.csv")
)
metadata_file <- Sys.getenv(
  "THESIS_METADATA_FILE",
  unset = file.path(data_dir, "demographics.txt")
)

run_nested_comparison <- env_flag("RUN_NESTED_COMPARISON", TRUE)
run_parameter_optimisation <- env_flag("RUN_PARAMETER_OPTIMISATION", TRUE)
run_panel_reduction <- env_flag("RUN_PANEL_REDUCTION", TRUE)
run_random_forest_check <- env_flag("RUN_RANDOM_FOREST_CHECK", TRUE)
run_functional_analysis <- env_flag("RUN_FUNCTIONAL_ANALYSIS", FALSE)

base_seed <- env_integer("THESIS_BASE_SEED", 20260611L)
outer_repeats <- env_integer("THESIS_OUTER_REPEATS", 15L)
panel_repeats <- env_integer("THESIS_PANEL_REPEATS", 50L)
outer_folds <- 5L
inner_folds <- 5L
positive_class <- "Infertile"
negative_class <- "Fertile"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "01_expression"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "02_model_comparison"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "03_elastic_net_tuning"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "04_panel_reduction"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "05_final_score"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "06_random_forest_check"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "07_functional_analysis"), recursive = TRUE, showWarnings = FALSE)

core_packages <- c(
  "edgeR", "limma", "caret", "glmnet", "randomForest", "rpart", "pROC"
)
missing_core <- core_packages[!vapply(core_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_core) > 0L) {
  stop(
    "Missing required R packages: ", paste(missing_core, collapse = ", "),
    "\nInstall CRAN packages with install.packages() and Bioconductor packages ",
    "with BiocManager::install()."
  )
}

###############################################################################
# 2. Reusable helpers
###############################################################################

write_table <- function(x, filename) {
  utils::write.csv(x, filename, row.names = FALSE, na = "")
}

rbind_fill <- function(rows) {
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(row) {
    missing <- setdiff(all_names, names(row))
    for (name in missing) row[[name]] <- NA
    row[, all_names, drop = FALSE]
  })
  do.call(rbind, rows)
}

normalise_ids <- function(x) gsub("\\.", "-", trimws(as.character(x)))

make_stratified_folds <- function(y, k = 5L, seed = 1L) {
  set.seed(seed)
  caret::createFolds(factor(y), k = k, returnTrain = FALSE)
}

make_foldid <- function(y, k = 5L, seed = 1L) {
  folds <- make_stratified_folds(y, k, seed)
  foldid <- integer(length(y))
  for (i in seq_along(folds)) foldid[folds[[i]]] <- i
  foldid
}

balanced_weights <- function(y) {
  y <- factor(y)
  counts <- table(y)
  weights <- length(y) / (length(counts) * counts)
  as.numeric(weights[as.character(y)])
}

auc_value <- function(truth, probability) {
  as.numeric(pROC::auc(pROC::roc(
    response = factor(truth, levels = c(negative_class, positive_class)),
    predictor = probability,
    levels = c(negative_class, positive_class),
    direction = "<",
    quiet = TRUE
  )))
}

classification_metrics <- function(truth, probability, threshold = 0.50) {
  truth <- factor(truth, levels = c(negative_class, positive_class))
  prediction <- factor(
    ifelse(probability >= threshold, positive_class, negative_class),
    levels = c(negative_class, positive_class)
  )
  tp <- sum(prediction == positive_class & truth == positive_class)
  tn <- sum(prediction == negative_class & truth == negative_class)
  fp <- sum(prediction == positive_class & truth == negative_class)
  fn <- sum(prediction == negative_class & truth == positive_class)
  sensitivity <- if ((tp + fn) == 0L) NA_real_ else tp / (tp + fn)
  specificity <- if ((tn + fp) == 0L) NA_real_ else tn / (tn + fp)
  data.frame(
    Threshold = threshold,
    Accuracy = (tp + tn) / length(truth),
    Sensitivity = sensitivity,
    Specificity = specificity,
    Balanced_Accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
    AUC = auc_value(truth, probability),
    TP = tp, TN = tn, FP = fp, FN = fn
  )
}

select_threshold <- function(truth, probability) {
  grid <- seq(0.05, 0.95, by = 0.01)
  table <- do.call(rbind, lapply(grid, function(threshold) {
    row <- classification_metrics(truth, probability, threshold)
    row$Youden_J <- row$Sensitivity + row$Specificity - 1
    row$Distance_from_0.50 <- abs(threshold - 0.50)
    row
  }))
  table <- table[order(
    -table$Youden_J,
    -table$Accuracy,
    table$Distance_from_0.50,
    table$Threshold
  ), , drop = FALSE]
  table[1, , drop = FALSE]
}

summarise_metrics <- function(table, group_columns) {
  split_key <- interaction(table[group_columns], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(table, split_key), function(data) {
    result <- data[1, group_columns, drop = FALSE]
    for (metric in c("AUC", "Balanced_Accuracy", "Accuracy", "Sensitivity", "Specificity")) {
      result[[paste0(metric, "_Mean")]] <- mean(data[[metric]], na.rm = TRUE)
      result[[paste0(metric, "_SD")]] <- stats::sd(data[[metric]], na.rm = TRUE)
    }
    result$N_Repeats <- nrow(data)
    result
  })
  combined <- do.call(rbind, rows)
  rownames(combined) <- NULL
  combined
}

best_tune_filter <- function(prediction_table, best_tune) {
  keep <- rep(TRUE, nrow(prediction_table))
  for (name in names(best_tune)) {
    if (name %in% names(prediction_table)) {
      if (is.numeric(best_tune[[name]])) {
        keep <- keep & abs(prediction_table[[name]] - best_tune[[name]]) < 1e-12
      } else {
        keep <- keep & prediction_table[[name]] == best_tune[[name]]
      }
    }
  }
  prediction_table[keep, , drop = FALSE]
}

fit_caret_model <- function(model_name, x_train, y_train, inner_index, seed) {
  inner_index_out <- lapply(
    inner_index,
    function(index) setdiff(seq_along(y_train), index)
  )
  control <- caret::trainControl(
    method = "cv",
    number = inner_folds,
    index = inner_index,
    indexOut = inner_index_out,
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = "final",
    sampling = "up",
    allowParallel = FALSE
  )
  set.seed(seed)
  if (model_name == "Elastic Net") {
    grid <- expand.grid(
      alpha = c(0.10, 0.50, 0.90, 1.00),
      lambda = exp(seq(log(1e-4), log(10), length.out = 15L))
    )
    caret::train(
      x = x_train, y = y_train, method = "glmnet", metric = "ROC",
      trControl = control, tuneGrid = grid, preProcess = c("center", "scale"),
      family = "binomial"
    )
  } else if (model_name == "Logistic Regression") {
    grid <- expand.grid(
      alpha = 0,
      lambda = exp(seq(log(1e-4), log(10), length.out = 20L))
    )
    caret::train(
      x = x_train, y = y_train, method = "glmnet", metric = "ROC",
      trControl = control, tuneGrid = grid, preProcess = c("center", "scale"),
      family = "binomial"
    )
  } else if (model_name == "Random Forest") {
    grid <- expand.grid(mtry = sort(unique(pmin(
      c(2L, 3L, 4L, 6L, 8L, ncol(x_train)), ncol(x_train)
    ))))
    caret::train(
      x = x_train, y = y_train, method = "rf", metric = "ROC",
      trControl = control, tuneGrid = grid, ntree = 500L,
      importance = FALSE
    )
  } else if (model_name == "Decision Tree") {
    grid <- expand.grid(cp = c(0, exp(seq(log(1e-4), log(0.316), length.out = 12L))))
    caret::train(
      x = x_train, y = y_train, method = "rpart", metric = "ROC",
      trControl = control, tuneGrid = grid
    )
  } else {
    stop("Unknown model: ", model_name)
  }
}

inner_oof_threshold <- function(fit) {
  predictions <- best_tune_filter(fit$pred, fit$bestTune)
  if (nrow(predictions) == 0L || !positive_class %in% names(predictions)) {
    stop("Could not recover inner out-of-fold probabilities from caret fit.")
  }
  select_threshold(predictions$obs, predictions[[positive_class]])
}

fit_cv_glmnet <- function(x, y, alpha, weighting, foldid, keep = TRUE) {
  y_numeric <- as.integer(y == positive_class)
  weights <- if (identical(weighting, "Balanced")) balanced_weights(y) else rep(1, length(y))
  glmnet::cv.glmnet(
    x = as.matrix(x), y = y_numeric, family = "binomial",
    alpha = alpha, weights = weights, foldid = foldid,
    standardize = TRUE, keep = keep, type.measure = "deviance"
  )
}

cv_oof_probability <- function(fit, lambda_rule) {
  selected_lambda <- fit[[lambda_rule]]
  lambda_index <- which.min(abs(fit$lambda - selected_lambda))
  probability <- as.numeric(fit$fit.preval[, lambda_index])
  if (any(probability < 0 | probability > 1, na.rm = TRUE)) probability <- stats::plogis(probability)
  probability
}

###############################################################################
# 3. Import, align, filter and test expression data
###############################################################################

if (!file.exists(counts_file)) stop("Counts file not found: ", counts_file)
if (!file.exists(metadata_file)) stop("Metadata file not found: ", metadata_file)

counts_data <- utils::read.csv(counts_file, check.names = FALSE)
if (ncol(counts_data) < 2L) stop("Counts file must contain a miRNA column and sample columns.")
mirna_ids <- make.unique(as.character(counts_data[[1]]))
counts <- as.matrix(counts_data[, -1, drop = FALSE])
storage.mode(counts) <- "numeric"
rownames(counts) <- mirna_ids
colnames(counts) <- normalise_ids(colnames(counts))

metadata <- utils::read.delim(metadata_file, check.names = FALSE)
required_columns <- c("Subject_ID", "Age", "BMI", "Ethnicity", "Status")
missing_columns <- setdiff(required_columns, names(metadata))
if (length(missing_columns) > 0L) {
  stop("Metadata is missing required columns: ", paste(missing_columns, collapse = ", "))
}
metadata$Subject_ID <- normalise_ids(metadata$Subject_ID)
metadata$Status <- trimws(as.character(metadata$Status))
metadata$Age <- as.numeric(metadata$Age)
metadata$BMI <- as.numeric(metadata$BMI)

common_ids <- intersect(metadata$Subject_ID, colnames(counts))
if (length(common_ids) == 0L) stop("No sample identifiers match between counts and metadata.")
metadata <- metadata[match(common_ids, metadata$Subject_ID), , drop = FALSE]
counts <- counts[, common_ids, drop = FALSE]
stopifnot(identical(metadata$Subject_ID, colnames(counts)))

labelled_index <- metadata$Status %in% c(negative_class, positive_class)
labelled_metadata <- metadata[labelled_index, , drop = FALSE]
labelled_counts <- counts[, labelled_index, drop = FALSE]
y <- factor(labelled_metadata$Status, levels = c(positive_class, negative_class))

if (sum(y == positive_class) < outer_folds) {
  stop("Too few Infertile samples for stratified five-fold analysis.")
}

dge <- edgeR::DGEList(counts = labelled_counts)
keep_expression <- rowSums(edgeR::cpm(dge) > 1) >= 15L
dge <- dge[keep_expression, , keep.lib.sizes = FALSE]
dge <- edgeR::calcNormFactors(dge, method = "TMM")

design_data <- data.frame(
  Status = stats::relevel(factor(labelled_metadata$Status), ref = negative_class),
  Age = labelled_metadata$Age,
  BMI = labelled_metadata$BMI
)
design <- stats::model.matrix(~ Age + BMI + Status, data = design_data)
voom_fit <- limma::voom(dge, design = design, plot = FALSE)
linear_fit <- limma::lmFit(voom_fit, design)
linear_fit <- limma::eBayes(linear_fit)
coefficient_name <- grep("^Status", colnames(design), value = TRUE)
if (length(coefficient_name) != 1L) stop("Could not identify the fertility-status coefficient.")
de_table <- limma::topTable(
  linear_fit, coef = coefficient_name, number = Inf, sort.by = "P", adjust.method = "BH"
)
de_table$miRNA <- rownames(de_table)
de_table <- de_table[, c("miRNA", setdiff(names(de_table), "miRNA")), drop = FALSE]
candidate_table <- de_table[!is.na(de_table$adj.P.Val) & de_table$adj.P.Val < 0.05, , drop = FALSE]
candidate_mirnas <- candidate_table$miRNA

if (length(candidate_mirnas) < 2L) stop("Fewer than two FDR-significant candidate miRNAs were found.")
x <- t(voom_fit$E[candidate_mirnas, , drop = FALSE])
storage.mode(x) <- "numeric"
rownames(x) <- labelled_metadata$Subject_ID

cohort_audit <- data.frame(
  Item = c(
    "Matched samples", "Fertile samples", "Borderline samples", "Infertile samples",
    "Annotated miRNA rows", "Expression-filtered miRNAs", "FDR-significant candidates"
  ),
  Value = c(
    nrow(metadata), sum(metadata$Status == negative_class), sum(metadata$Status == "Borderline"),
    sum(metadata$Status == positive_class), nrow(counts), nrow(dge), length(candidate_mirnas)
  )
)
write_table(cohort_audit, file.path(output_dir, "01_expression", "cohort_and_filter_audit.csv"))
write_table(de_table, file.path(output_dir, "01_expression", "differential_expression_all.csv"))
write_table(candidate_table, file.path(output_dir, "01_expression", "candidate_mirnas_fdr_lt_0.05.csv"))
write_table(
  data.frame(SampleID = rownames(x), Status = as.character(y), x, check.names = FALSE),
  file.path(output_dir, "01_expression", "candidate_expression_matrix_labelled.csv")
)

cat("Matched samples:", nrow(metadata), "\n")
cat("Labelled samples:", nrow(x), "(", sum(y == positive_class), "Infertile )\n")
cat("Expression-filtered miRNAs:", nrow(dge), "\n")
cat("FDR-significant candidates:", length(candidate_mirnas), "\n")

###############################################################################
# 4. Repeated nested five-fold comparison of four model classes
###############################################################################

model_names <- c("Elastic Net", "Random Forest", "Logistic Regression", "Decision Tree")
model_comparison_summary <- NULL

if (run_nested_comparison) {
  cat("\nRunning nested model comparison...\n")
  prediction_rows <- list()
  metric_rows <- list()
  tuning_rows <- list()
  prediction_index <- 1L
  metric_index <- 1L
  tuning_index <- 1L

  for (repeat_id in seq_len(outer_repeats)) {
    folds <- make_stratified_folds(y, outer_folds, base_seed + repeat_id)
    for (model_id in seq_along(model_names)) {
      model_name <- model_names[[model_id]]
      repeat_probability <- rep(NA_real_, length(y))
      repeat_prediction <- rep(NA_character_, length(y))
      repeat_fold <- rep(NA_integer_, length(y))
      repeat_threshold <- rep(NA_real_, length(y))

      for (fold_id in seq_along(folds)) {
        test_index <- folds[[fold_id]]
        train_index <- setdiff(seq_along(y), test_index)
        set.seed(base_seed + repeat_id * 100L + fold_id)
        inner_index <- caret::createFolds(
          droplevels(y[train_index]),
          k = inner_folds,
          returnTrain = TRUE
        )
        fit_seed <- base_seed + repeat_id * 10000L + fold_id * 100L + model_id
        fit <- fit_caret_model(
          model_name,
          x[train_index, , drop = FALSE],
          droplevels(y[train_index]),
          inner_index,
          fit_seed
        )
        threshold_row <- inner_oof_threshold(fit)
        probability <- stats::predict(
          fit, newdata = x[test_index, , drop = FALSE], type = "prob"
        )[[positive_class]]
        repeat_probability[test_index] <- probability
        repeat_prediction[test_index] <- ifelse(
          probability >= threshold_row$Threshold, positive_class, negative_class
        )
        repeat_fold[test_index] <- fold_id
        repeat_threshold[test_index] <- threshold_row$Threshold

        tuning_row <- data.frame(
          Repeat = repeat_id, Outer_Fold = fold_id, Model = model_name,
          Inner_Selected_Threshold = threshold_row$Threshold
        )
        for (name in names(fit$bestTune)) tuning_row[[name]] <- fit$bestTune[[name]]
        tuning_rows[[tuning_index]] <- tuning_row
        tuning_index <- tuning_index + 1L
      }

      prediction_rows[[prediction_index]] <- data.frame(
        Repeat = repeat_id,
        Model = model_name,
        SampleID = labelled_metadata$Subject_ID,
        Truth = as.character(y),
        Outer_Fold = repeat_fold,
        Inner_Selected_Threshold = repeat_threshold,
        Probability_Infertile = repeat_probability,
        Predicted_Class = repeat_prediction
      )
      metric_rows[[metric_index]] <- cbind(
        data.frame(Repeat = repeat_id, Model = model_name),
        classification_metrics(y, repeat_probability, threshold = 0.50)[, -1, drop = FALSE]
      )
      # Threshold-dependent metrics use the fold-specific inner thresholds.
      fold_specific <- factor(repeat_prediction, levels = c(negative_class, positive_class))
      truth <- factor(y, levels = c(negative_class, positive_class))
      tp <- sum(fold_specific == positive_class & truth == positive_class)
      tn <- sum(fold_specific == negative_class & truth == negative_class)
      fp <- sum(fold_specific == positive_class & truth == negative_class)
      fn <- sum(fold_specific == negative_class & truth == positive_class)
      metric_rows[[metric_index]]$Accuracy <- (tp + tn) / length(truth)
      metric_rows[[metric_index]]$Sensitivity <- tp / (tp + fn)
      metric_rows[[metric_index]]$Specificity <- tn / (tn + fp)
      metric_rows[[metric_index]]$Balanced_Accuracy <- mean(c(
        metric_rows[[metric_index]]$Sensitivity,
        metric_rows[[metric_index]]$Specificity
      ))
      metric_rows[[metric_index]]$TP <- tp
      metric_rows[[metric_index]]$TN <- tn
      metric_rows[[metric_index]]$FP <- fp
      metric_rows[[metric_index]]$FN <- fn
      prediction_index <- prediction_index + 1L
      metric_index <- metric_index + 1L
    }
    cat("  completed repeat", repeat_id, "of", outer_repeats, "\n")
  }

  model_predictions <- do.call(rbind, prediction_rows)
  model_metrics <- do.call(rbind, metric_rows)
  model_tuning <- rbind_fill(tuning_rows)
  model_comparison_summary <- summarise_metrics(model_metrics, "Model")
  model_comparison_summary <- model_comparison_summary[order(
    -model_comparison_summary$AUC_Mean,
    -model_comparison_summary$Balanced_Accuracy_Mean
  ), , drop = FALSE]

  write_table(model_predictions, file.path(
    output_dir, "02_model_comparison", "nested_outer_oof_predictions.csv"
  ))
  write_table(model_metrics, file.path(
    output_dir, "02_model_comparison", "nested_metrics_by_repeat.csv"
  ))
  write_table(model_tuning, file.path(
    output_dir, "02_model_comparison", "nested_selected_settings.csv"
  ))
  write_table(model_comparison_summary, file.path(
    output_dir, "02_model_comparison", "nested_model_comparison_summary.csv"
  ))
}

###############################################################################
# 5. Elastic Net parameter optimisation after model-class selection
###############################################################################

selected_alpha <- 0.25
selected_lambda_rule <- "lambda.min"
selected_weighting <- "Balanced"
parameter_summary <- NULL

if (run_parameter_optimisation) {
  cat("\nRunning Elastic Net parameter optimisation...\n")
  alpha_grid <- c(0.10, 0.25, 0.50, 0.75, 0.90)
  lambda_rule_grid <- c("lambda.min", "lambda.1se")
  weighting_grid <- c("Balanced", "None")
  tuning_rows <- list()
  tuning_index <- 1L

  for (repeat_id in seq_len(outer_repeats)) {
    foldid <- make_foldid(y, inner_folds, 50000L + base_seed + repeat_id - 1L)
    for (alpha in alpha_grid) {
      for (weighting in weighting_grid) {
        fit <- fit_cv_glmnet(x, y, alpha, weighting, foldid, keep = TRUE)
        for (lambda_rule in lambda_rule_grid) {
          probability <- cv_oof_probability(fit, lambda_rule)
          threshold <- select_threshold(y, probability)$Threshold
          metrics <- classification_metrics(y, probability, threshold)
          tuning_rows[[tuning_index]] <- cbind(
            data.frame(
              Repeat = repeat_id,
              Alpha = alpha,
              Lambda_Rule = lambda_rule,
              Class_Weighting = weighting
            ),
            metrics
          )
          tuning_index <- tuning_index + 1L
        }
      }
    }
  }

  parameter_metrics <- do.call(rbind, tuning_rows)
  parameter_summary <- summarise_metrics(
    parameter_metrics,
    c("Alpha", "Lambda_Rule", "Class_Weighting")
  )
  parameter_summary <- parameter_summary[order(
    -parameter_summary$AUC_Mean,
    parameter_summary$AUC_SD,
    abs(parameter_summary$Alpha - 0.50),
    parameter_summary$Alpha
  ), , drop = FALSE]
  selected_alpha <- parameter_summary$Alpha[[1]]
  selected_lambda_rule <- parameter_summary$Lambda_Rule[[1]]
  selected_weighting <- parameter_summary$Class_Weighting[[1]]
  parameter_summary$Selected <- FALSE
  parameter_summary$Selected[[1]] <- TRUE

  write_table(parameter_metrics, file.path(
    output_dir, "03_elastic_net_tuning", "parameter_metrics_by_repeat.csv"
  ))
  write_table(parameter_summary, file.path(
    output_dir, "03_elastic_net_tuning", "parameter_optimisation_summary.csv"
  ))
}

cat(
  "Selected Elastic Net setting: alpha =", selected_alpha,
  ",", selected_lambda_rule, ",", selected_weighting, "weights\n"
)

###############################################################################
# 6. Stability analysis, leave-one-out ablation and panel-size comparison
###############################################################################

fit_locked_enet <- function(x_train, y_train, seed, keep = FALSE) {
  foldid <- make_foldid(y_train, inner_folds, seed)
  fit_cv_glmnet(
    x_train, y_train, selected_alpha, selected_weighting, foldid, keep = keep
  )
}

evaluate_enet_panel <- function(panel, repeats = panel_repeats, threshold = 0.50) {
  metric_rows <- vector("list", repeats)
  prediction_rows <- vector("list", repeats)
  for (repeat_id in seq_len(repeats)) {
    seed <- base_seed + repeat_id - 1L
    folds <- make_stratified_folds(y, outer_folds, seed)
    probability <- rep(NA_real_, length(y))
    for (fold_id in seq_along(folds)) {
      test_index <- folds[[fold_id]]
      train_index <- setdiff(seq_along(y), test_index)
      fit <- fit_locked_enet(
        x[train_index, panel, drop = FALSE],
        droplevels(y[train_index]),
        seed + fold_id * 100L
      )
      probability[test_index] <- as.numeric(stats::predict(
        fit,
        newx = as.matrix(x[test_index, panel, drop = FALSE]),
        s = selected_lambda_rule,
        type = "response"
      ))
    }
    metric_rows[[repeat_id]] <- cbind(
      data.frame(Repeat = repeat_id, Seed = seed, Panel_Size = length(panel)),
      classification_metrics(y, probability, threshold)
    )
    prediction_rows[[repeat_id]] <- data.frame(
      Repeat = repeat_id,
      Seed = seed,
      SampleID = labelled_metadata$Subject_ID,
      Truth = as.character(y),
      Probability_Infertile = probability
    )
  }
  list(
    metrics = do.call(rbind, metric_rows),
    predictions = do.call(rbind, prediction_rows)
  )
}

final_panel <- candidate_mirnas

if (run_panel_reduction) {
  cat("\nRunning Elastic Net stability and panel reduction...\n")
  coefficient_rows <- list()
  coefficient_index <- 1L
  for (repeat_id in seq_len(panel_repeats)) {
    seed <- base_seed + repeat_id - 1L
    folds <- make_stratified_folds(y, outer_folds, seed)
    for (fold_id in seq_along(folds)) {
      test_index <- folds[[fold_id]]
      train_index <- setdiff(seq_along(y), test_index)
      fit <- fit_locked_enet(
        x[train_index, , drop = FALSE],
        droplevels(y[train_index]),
        seed + fold_id * 100L
      )
      coefficient_matrix <- as.matrix(stats::coef(fit, s = selected_lambda_rule))
      coefficient_rows[[coefficient_index]] <- data.frame(
        Repeat = repeat_id,
        Outer_Fold = fold_id,
        miRNA = rownames(coefficient_matrix)[-1],
        Coefficient = as.numeric(coefficient_matrix[-1, 1])
      )
      coefficient_index <- coefficient_index + 1L
    }
  }
  stability_coefficients <- do.call(rbind, coefficient_rows)
  stability_split <- split(stability_coefficients, stability_coefficients$miRNA)
  stability_summary <- do.call(rbind, lapply(stability_split, function(data) {
    data.frame(
      miRNA = data$miRNA[[1]],
      Fits = nrow(data),
      Selected_Fits = sum(data$Coefficient != 0),
      Selection_Frequency = mean(data$Coefficient != 0),
      Mean_Coefficient = mean(data$Coefficient),
      Mean_Absolute_Coefficient = mean(abs(data$Coefficient))
    )
  }))
  stability_summary <- stability_summary[order(
    -stability_summary$Selection_Frequency,
    -stability_summary$Mean_Absolute_Coefficient,
    stability_summary$miRNA
  ), , drop = FALSE]
  rownames(stability_summary) <- NULL

  # Leave-one-miRNA-out comparisons use the same repeated folds as the full set.
  ablation_rows <- list()
  panels <- c(list("Full candidate set" = candidate_mirnas), setNames(
    lapply(candidate_mirnas, function(mirna) setdiff(candidate_mirnas, mirna)),
    paste0("Without ", candidate_mirnas)
  ))
  for (panel_name in names(panels)) {
    result <- evaluate_enet_panel(panels[[panel_name]])
    summary <- summarise_metrics(result$metrics, "Panel_Size")
    summary$Panel <- panel_name
    ablation_rows[[panel_name]] <- summary
  }
  ablation_summary <- do.call(rbind, ablation_rows)
  rownames(ablation_summary) <- NULL

  # Cumulative panels are defined by stability rank; not every subset is tested.
  ranked_mirnas <- stability_summary$miRNA
  panel_size_rows <- list()
  panel_size_metrics <- list()
  for (size in 2:length(ranked_mirnas)) {
    panel <- ranked_mirnas[seq_len(size)]
    result <- evaluate_enet_panel(panel)
    summary <- summarise_metrics(result$metrics, "Panel_Size")
    summary$Panel <- paste(panel, collapse = "; ")
    panel_size_rows[[as.character(size)]] <- summary
    result$metrics$Panel <- paste(panel, collapse = "; ")
    panel_size_metrics[[as.character(size)]] <- result$metrics
  }
  panel_size_summary <- do.call(rbind, panel_size_rows)
  rownames(panel_size_summary) <- NULL
  panel_size_summary <- panel_size_summary[order(
    -panel_size_summary$AUC_Mean,
    panel_size_summary$Panel_Size
  ), , drop = FALSE]
  final_panel <- strsplit(panel_size_summary$Panel[[1]], "; ", fixed = TRUE)[[1]]

  write_table(stability_coefficients, file.path(
    output_dir, "04_panel_reduction", "stability_coefficients_250_fits.csv"
  ))
  write_table(stability_summary, file.path(
    output_dir, "04_panel_reduction", "stability_summary.csv"
  ))
  write_table(ablation_summary, file.path(
    output_dir, "04_panel_reduction", "leave_one_mirna_out_summary.csv"
  ))
  write_table(do.call(rbind, panel_size_metrics), file.path(
    output_dir, "04_panel_reduction", "panel_size_metrics_by_repeat.csv"
  ))
  write_table(panel_size_summary, file.path(
    output_dir, "04_panel_reduction", "panel_size_summary.csv"
  ))
  write_table(data.frame(miRNA = final_panel), file.path(
    output_dir, "04_panel_reduction", "selected_final_panel.csv"
  ))
}

cat("Selected panel:", paste(final_panel, collapse = ", "), "\n")

###############################################################################
# 7. Repeated OOF score and post-development Borderline projection
###############################################################################

# The score stage uses the explicitly recorded TMM-normalised log2-CPM scale
# (prior count 1) across the supplied libraries. Status labels are not used in
# this normalisation. A prospective assay would require a frozen reference.
dge_all <- edgeR::DGEList(counts = counts)
dge_all <- edgeR::calcNormFactors(dge_all, method = "TMM")
logcpm_all <- edgeR::cpm(dge_all, log = TRUE, prior.count = 1)
if (!all(final_panel %in% rownames(logcpm_all))) {
  stop("Final-panel miRNAs are absent from the full expression matrix.")
}
x_score_all <- t(logcpm_all[final_panel, , drop = FALSE])
x_score <- x_score_all[labelled_metadata$Subject_ID, , drop = FALSE]

score_prediction_rows <- vector("list", panel_repeats)
score_metric_rows <- vector("list", panel_repeats)
for (repeat_id in seq_len(panel_repeats)) {
  seed <- base_seed + repeat_id - 1L
  folds <- make_stratified_folds(y, outer_folds, seed)
  probability <- rep(NA_real_, length(y))
  fold_number <- rep(NA_integer_, length(y))
  for (fold_id in seq_along(folds)) {
    test_index <- folds[[fold_id]]
    train_index <- setdiff(seq_along(y), test_index)
    fit <- fit_locked_enet(
      x_score[train_index, , drop = FALSE],
      droplevels(y[train_index]),
      seed + fold_id * 100L
    )
    probability[test_index] <- as.numeric(stats::predict(
      fit,
      newx = as.matrix(x_score[test_index, , drop = FALSE]),
      s = selected_lambda_rule,
      type = "response"
    ))
    fold_number[test_index] <- fold_id
  }
  score_prediction_rows[[repeat_id]] <- data.frame(
    Repeat = repeat_id,
    Seed = seed,
    SampleID = labelled_metadata$Subject_ID,
    Truth = as.character(y),
    Fold = fold_number,
    Probability_Infertile = probability
  )
  score_metric_rows[[repeat_id]] <- cbind(
    data.frame(Repeat = repeat_id, Seed = seed),
    classification_metrics(y, probability, 0.50)
  )
}

score_oof <- do.call(rbind, score_prediction_rows)
score_metrics <- do.call(rbind, score_metric_rows)
mean_oof <- aggregate(
  Probability_Infertile ~ SampleID + Truth,
  data = score_oof,
  FUN = mean
)
mean_oof <- mean_oof[match(labelled_metadata$Subject_ID, mean_oof$SampleID), , drop = FALSE]
operating_point <- select_threshold(mean_oof$Truth, mean_oof$Probability_Infertile)
mean_oof_metrics <- classification_metrics(
  mean_oof$Truth,
  mean_oof$Probability_Infertile,
  operating_point$Threshold
)
mean_oof$Score_0_100 <- 100 * mean_oof$Probability_Infertile
mean_oof$Predicted_Class <- ifelse(
  mean_oof$Probability_Infertile >= operating_point$Threshold,
  positive_class,
  negative_class
)

final_foldid <- make_foldid(y, inner_folds, 99999L)
final_fit <- fit_locked_enet(x_score, y, 99999L)
final_coefficients <- as.matrix(stats::coef(final_fit, s = selected_lambda_rule))
coefficient_table <- data.frame(
  Term = rownames(final_coefficients),
  Coefficient = as.numeric(final_coefficients[, 1]),
  Nonzero = as.numeric(final_coefficients[, 1]) != 0
)
all_probability <- as.numeric(stats::predict(
  final_fit,
  newx = as.matrix(x_score_all),
  s = selected_lambda_rule,
  type = "response"
))
all_scores <- data.frame(
  SampleID = rownames(x_score_all),
  Status = metadata$Status[match(rownames(x_score_all), metadata$Subject_ID)],
  Probability_Infertile = all_probability,
  Score_0_100 = 100 * all_probability,
  Cutoff_Probability = operating_point$Threshold,
  Cutoff_Score = 100 * operating_point$Threshold
)

write_table(score_oof, file.path(output_dir, "05_final_score", "repeated_oof_predictions.csv"))
write_table(score_metrics, file.path(output_dir, "05_final_score", "repeated_oof_metrics.csv"))
write_table(mean_oof, file.path(output_dir, "05_final_score", "participant_mean_oof_scores.csv"))
write_table(mean_oof_metrics, file.path(output_dir, "05_final_score", "mean_oof_operating_metrics.csv"))
write_table(operating_point, file.path(output_dir, "05_final_score", "youden_operating_point.csv"))
write_table(coefficient_table, file.path(output_dir, "05_final_score", "final_model_coefficients.csv"))
write_table(all_scores, file.path(output_dir, "05_final_score", "all_sample_scores_local_only.csv"))
saveRDS(final_fit, file.path(output_dir, "05_final_score", "final_development_model.rds"))

###############################################################################
# 8. Random forest sensitivity check
###############################################################################

fit_predict_rf <- function(x_train, y_train, x_test, seed, mtry = NULL) {
  if (is.null(mtry)) mtry <- min(5L, ncol(x_train))
  class_table <- table(y_train)
  class_weights <- length(y_train) / (length(class_table) * class_table)
  set.seed(seed)
  fit <- randomForest::randomForest(
    x = as.data.frame(x_train),
    y = y_train,
    ntree = 500L,
    mtry = min(mtry, ncol(x_train)),
    nodesize = 1L,
    classwt = class_weights,
    importance = FALSE
  )
  as.numeric(
    stats::predict(fit, newdata = as.data.frame(x_test), type = "prob")[, positive_class]
  )
}

evaluate_rf_panel <- function(panel, repeats = panel_repeats) {
  rows <- vector("list", repeats)
  for (repeat_id in seq_len(repeats)) {
    seed <- 8000L + repeat_id
    folds <- make_stratified_folds(y, outer_folds, seed)
    probability <- rep(NA_real_, length(y))
    for (fold_id in seq_along(folds)) {
      test_index <- folds[[fold_id]]
      train_index <- setdiff(seq_along(y), test_index)
      probability[test_index] <- fit_predict_rf(
        x[train_index, panel, drop = FALSE],
        droplevels(y[train_index]),
        x[test_index, panel, drop = FALSE],
        seed + fold_id * 100L
      )
    }
    rows[[repeat_id]] <- cbind(
      data.frame(Repeat = repeat_id, Panel_Size = length(panel)),
      classification_metrics(y, probability, 0.50)
    )
  }
  do.call(rbind, rows)
}

if (run_random_forest_check) {
  cat("\nRunning random forest panel sensitivity check...\n")
  # Greedy backward deletion is model-specific and reuses the same development cohort.
  cache <- new.env(parent = emptyenv())
  cached_rf <- function(panel) {
    key <- paste(sort(panel), collapse = "|")
    if (!exists(key, envir = cache, inherits = FALSE)) {
      assign(key, evaluate_rf_panel(panel), envir = cache)
    }
    get(key, envir = cache, inherits = FALSE)
  }

  current_panel <- candidate_mirnas
  greedy_path <- list()
  step <- 1L
  initial_metrics <- cached_rf(current_panel)
  initial_summary <- summarise_metrics(initial_metrics, "Panel_Size")
  initial_summary$Step <- step
  initial_summary$Removed <- "None"
  initial_summary$Panel <- paste(current_panel, collapse = "; ")
  greedy_path[[step]] <- initial_summary

  while (length(current_panel) > 1L) {
    candidate_rows <- lapply(current_panel, function(removed) {
      panel <- setdiff(current_panel, removed)
      summary <- summarise_metrics(cached_rf(panel), "Panel_Size")
      summary$Removed <- removed
      summary$Panel <- paste(panel, collapse = "; ")
      summary
    })
    candidate_table <- do.call(rbind, candidate_rows)
    maximum_auc <- max(candidate_table$AUC_Mean)
    eligible <- candidate_table[candidate_table$AUC_Mean >= maximum_auc - 0.005, , drop = FALSE]
    selected <- eligible[order(
      -eligible$Balanced_Accuracy_Mean,
      -eligible$AUC_Mean,
      eligible$Panel_Size
    ), , drop = FALSE][1, , drop = FALSE]
    current_panel <- strsplit(selected$Panel, "; ", fixed = TRUE)[[1]]
    step <- step + 1L
    selected$Step <- step
    greedy_path[[step]] <- selected
  }
  rf_greedy_path <- do.call(rbind, greedy_path)
  rownames(rf_greedy_path) <- NULL
  rf_greedy_path <- rf_greedy_path[order(rf_greedy_path$Panel_Size), , drop = FALSE]
  best_rf <- rf_greedy_path[order(
    -rf_greedy_path$AUC_Mean,
    rf_greedy_path$Panel_Size
  ), , drop = FALSE][1, , drop = FALSE]

  rf_comparison_panels <- list(
    "All candidates" = candidate_mirnas,
    "Elastic Net selected panel" = final_panel,
    "Random forest selected panel" = strsplit(best_rf$Panel, "; ", fixed = TRUE)[[1]]
  )
  rf_rows <- lapply(names(rf_comparison_panels), function(name) {
    metrics <- cached_rf(rf_comparison_panels[[name]])
    summary <- summarise_metrics(metrics, "Panel_Size")
    summary$Panel_Name <- name
    summary$Members <- paste(rf_comparison_panels[[name]], collapse = "; ")
    summary
  })
  rf_comparison <- do.call(rbind, rf_rows)

  write_table(rf_greedy_path, file.path(
    output_dir, "06_random_forest_check", "random_forest_backward_ablation_path.csv"
  ))
  write_table(rf_comparison, file.path(
    output_dir, "06_random_forest_check", "random_forest_panel_comparison.csv"
  ))
}

###############################################################################
# 9. Optional multiMiR target retrieval and KEGG enrichment
###############################################################################

if (run_functional_analysis) {
  functional_packages <- c("multiMiR", "clusterProfiler")
  missing_functional <- functional_packages[
    !vapply(functional_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_functional) > 0L) {
    stop("Functional analysis requested but packages are missing: ", paste(missing_functional, collapse = ", "))
  }
  cat("\nRetrieving validated multiMiR targets and running KEGG enrichment...\n")
  target_result <- multiMiR::get_multimir(
    mirna = final_panel,
    table = "validated",
    summary = FALSE,
    legacy.out = FALSE
  )
  target_data <- target_result@data
  mirna_column <- intersect(
    c("mature_mirna_id", "mirna", "mature_mirna"), names(target_data)
  )[[1]]
  entrez_column <- intersect(
    c("target_entrez", "target_entrez_id", "entrez"), names(target_data)
  )[[1]]
  symbol_column <- intersect(
    c("target_symbol", "symbol", "target_gene"), names(target_data)
  )[[1]]
  targets <- target_data[, unique(c(mirna_column, symbol_column, entrez_column)), drop = FALSE]
  names(targets)[match(c(mirna_column, symbol_column, entrez_column), names(targets))] <-
    c("miRNA", "Target_Symbol", "Entrez_ID")
  targets$Entrez_ID <- as.character(targets$Entrez_ID)
  targets <- targets[
    !is.na(targets$Entrez_ID) & nzchar(targets$Entrez_ID), , drop = FALSE
  ]
  targets <- unique(targets)
  write_table(targets, file.path(
    output_dir, "07_functional_analysis", "validated_multimir_targets.csv"
  ))

  run_kegg <- function(entrez_ids, label) {
    result <- clusterProfiler::enrichKEGG(
      gene = unique(entrez_ids), organism = "hsa", pAdjustMethod = "BH",
      pvalueCutoff = 1, qvalueCutoff = 1
    )
    table <- as.data.frame(result)
    if (nrow(table) == 0L) return(table)
    table$Target_Set <- label
    table$Significant_FDR_0.05 <- !is.na(table$p.adjust) & table$p.adjust < 0.05
    broad_label <- grepl(
      "cancer|carcinoma|leukemia|infection|viral|bacterial|inflammatory disease",
      table$Description,
      ignore.case = TRUE
    )
    table$Retained_in_Focused_Display <- table$Significant_FDR_0.05 & !broad_label
    table
  }

  kegg_rows <- lapply(split(targets, targets$miRNA), function(data) {
    run_kegg(data$Entrez_ID, data$miRNA[[1]])
  })
  kegg_rows[["Pooled selected panel"]] <- run_kegg(
    targets$Entrez_ID, "Pooled selected panel"
  )
  kegg_table <- do.call(rbind, kegg_rows)
  write_table(kegg_table, file.path(
    output_dir, "07_functional_analysis", "kegg_complete_and_focused_results.csv"
  ))
}

###############################################################################
# 10. Reproducibility record
###############################################################################

configuration <- data.frame(
  Setting = c(
    "R version", "Base seed", "Outer folds", "Inner folds",
    "Nested comparison repeats", "Panel/score repeats", "Positive class",
    "Candidate threshold", "Elastic Net alpha", "Elastic Net lambda rule",
    "Class weighting", "Displayed operating-point rule", "Final panel"
  ),
  Value = c(
    R.version.string, base_seed, outer_folds, inner_folds,
    outer_repeats, panel_repeats, positive_class,
    "BH FDR < 0.05 after age/BMI-adjusted limma-voom model",
    selected_alpha, selected_lambda_rule, selected_weighting,
    "Youden index on participant-level mean repeated OOF probabilities",
    paste(final_panel, collapse = "; ")
  )
)
write_table(configuration, file.path(output_dir, "analysis_configuration.csv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

cat("\nAnalysis complete. Outputs written to:\n", normalizePath(output_dir), "\n")
