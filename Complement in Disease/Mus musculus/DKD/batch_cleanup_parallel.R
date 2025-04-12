# === batch_cleanup_parallel.R ===

library(DropletUtils)
library(SingleCellExperiment)
library(SoupX)
library(parallel)  # For parallel processing

# === SETTINGS ===
input_dir <- "/workspaces/Biochemistry-BS/Complement in Disease/Mus musculus/DKD/raw_data"
output_dir <- "/workspaces/Biochemistry-BS/Complement in Disease/Mus musculus/DKD/cleaned_data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Log file for tracking progress
log_file <- file.path(output_dir, "processing_log.txt")
cat("Starting batch processing...\n", file = log_file)

samples <- list(
  "GSM5594468" = "GSM5594468_E3019_raw_feature_bc_matrix.h5",
  "GSM5594469" = "GSM5594469_A3020_raw_feature_bc_matrix.h5",
  "GSM5594470" = "GSM5594470_F3021_raw_feature_bc_matrix.h5",
  "GSM5594471" = "GSM5594471_G3022_raw_feature_bc_matrix.h5",
  "GSM5594472" = "GSM5594472_H3023_raw_feature_bc_matrix.h5",
  "GSM5594473" = "GSM5594473_N3024_raw_feature_bc_matrix.h5",
  "GSM5594474" = "GSM5594474_B3025_raw_feature_bc_matrix.h5",
  "GSM5594475" = "GSM5594475_A3026_raw_feature_bc_matrix.h5",
  "GSM5594476" = "GSM5594476_B3027_raw_feature_bc_matrix.h5",
  "GSM5594477" = "GSM5594477_C3028_raw_feature_bc_matrix.h5",
  "GSM5594478" = "GSM5594478_A3029_raw_feature_bc_matrix.h5"
)

# === FUNCTION TO PROCESS A SINGLE SAMPLE ===
process_sample <- function(sample_id, file_name) {
  input_file <- file.path(input_dir, file_name)
  out_path <- file.path(output_dir, sample_id)

  cat("\n🔬 Processing", sample_id, "\n", file = log_file, append = TRUE)

  if (!file.exists(input_file)) {
    cat("⚠️  File not found:", input_file, "\n", file = log_file, append = TRUE)
    return(NULL)
  }

  tryCatch({
    # 1. Read 10x data
    cat("📖 Reading 10x data...\n", file = log_file, append = TRUE)
    sce <- read10xCounts(input_file)

    # 2. Run EmptyDrops
    cat("📉 Running EmptyDrops...\n", file = log_file, append = TRUE)
    ed <- emptyDrops(counts(sce))
    keep <- which(ed$FDR < 0.01)
    sce_filtered <- sce[, keep]

    # 3. Run SoupX
    cat("🧪 Running SoupX...\n", file = log_file, append = TRUE)
    sc <- SoupChannel(tod = counts(sce), toc = counts(sce_filtered))
    sc <- autoEstCont(sc)
    cleaned_counts <- adjustCounts(sc)

    # 4. Save
    cat("💾 Saving cleaned data...\n", file = log_file, append = TRUE)
    dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
    DropletUtils::write10xCounts(out_path, cleaned_counts, version = "3")
    cat("✅ Saved:", out_path, "\n", file = log_file, append = TRUE)

  }, error = function(e) {
    cat("❌ Error processing", sample_id, ":", e$message, "\n", file = log_file, append = TRUE)
  })
}

# === PARALLEL PROCESSING ===
num_cores <- detectCores() - 1  # Use 3 cores, reserve 1 for the system
cat("🖥️ Using", num_cores, "cores for parallel processing\n", file = log_file, append = TRUE)

mclapply(names(samples), function(sample_id) {
  process_sample(sample_id, samples[[sample_id]])
}, mc.cores = num_cores)

cat("✅ Batch processing completed!\n", file = log_file, append = TRUE)