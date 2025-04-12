# === batch_cleanup.R ===

library(DropletUtils)
library(SingleCellExperiment)
library(SoupX)

# === SETTINGS ===
input_dir <- "/workspaces/Biochemistry-BS/Complement in Disease/Mus musculus/DKD/raw_data"
output_dir <- "/workspaces/Biochemistry-BS/Complement in Disease/Mus musculus/DKD/cleaned_data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Log file for tracking progress
log_file <- file.path(output_dir, "processing_log.txt")
cat("Starting batch processing...\n", file = log_file)

# Test with just one sample
samples <- list(
  "GSM5594468" = "GSM5594468_E3019_raw_feature_bc_matrix.h5"
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
    cat("✅ Successfully read 10x data for", sample_id, "\n", file = log_file, append = TRUE)

    # 2. Run EmptyDrops
    cat("📉 Running EmptyDrops...\n", file = log_file, append = TRUE)
    ed <- emptyDrops(counts(sce))
    keep <- which(ed$FDR < 0.01)
    sce_filtered <- sce[, keep]
    cat("✅ Successfully filtered data for", sample_id, "\n", file = log_file, append = TRUE)

    # 3. Run SoupX
    cat("🧪 Running SoupX...\n", file = log_file, append = TRUE)
    sc <- SoupChannel(tod = counts(sce), toc = counts(sce_filtered))
    sc <- autoEstCont(sc)
    cleaned_counts <- adjustCounts(sc)
    cat("✅ Successfully ran SoupX for", sample_id, "\n", file = log_file, append = TRUE)

    # 4. Save
    cat("💾 Saving cleaned data...\n", file = log_file, append = TRUE)
    dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
    DropletUtils::write10xCounts(out_path, cleaned_counts, version = "3")
    cat("✅ Saved cleaned data for", sample_id, "to", out_path, "\n", file = log_file, append = TRUE)

  }, error = function(e) {
    cat("❌ Error processing", sample_id, ":", e$message, "\n", file = log_file, append = TRUE)
  })
}

# === SEQUENTIAL PROCESSING ===
cat("🖥️ Processing one sample sequentially\n", file = log_file, append = TRUE)

lapply(names(samples), function(sample_id) {
  process_sample(sample_id, samples[[sample_id]])
})

cat("✅ Batch processing completed!\n", file = log_file, append = TRUE)