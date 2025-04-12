# === batch_cleanup.R ===

library(DropletUtils)
library(SingleCellExperiment)
library(SoupX)

# === SETTINGS ===
input_dir <- "/workspaces/Biochemistry-BS/Complement in Disease/Mus musculus/DKD/raw_data"  # Path relative to repo root (inside Codespaces)
output_dir <- "/workspaces/Biochemistry-BS/Complement in Disease/Mus musculus/DKD/cleaned_data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

samples <- list(
  "GSM5594468" = "GSM5594468_E3019_raw_feature_bc_matrix.h5"
  # "GSM5594469" = "GSM5594469_A3020_raw_feature_bc_matrix.h5",
  # "GSM5594470" = "GSM5594470_F3021_raw_feature_bc_matrix.h5",
  # "GSM5594471" = "GSM5594471_G3022_raw_feature_bc_matrix.h5",
  # "GSM5594472" = "GSM5594472_H3023_raw_feature_bc_matrix.h5",
  # "GSM5594473" = "GSM5594473_N3024_raw_feature_bc_matrix.h5",
  # "GSM5594474" = "GSM5594474_B3025_raw_feature_bc_matrix.h5",
  # "GSM5594475" = "GSM5594475_A3026_raw_feature_bc_matrix.h5",
  # "GSM5594476" = "GSM5594476_B3027_raw_feature_bc_matrix.h5",
  # "GSM5594477" = "GSM5594477_C3028_raw_feature_bc_matrix.h5",
  # "GSM5594478" = "GSM5594478_A3029_raw_feature_bc_matrix.h5"
)

# === PROCESS SAMPLES ONE AT A TIME ===
for (sample_id in names(samples)) {
  input_file <- file.path(input_dir, samples[[sample_id]])
  out_path <- file.path(output_dir, sample_id)

  cat("\n🔬 Processing", sample_id, "\n")

  if (!file.exists(input_file)) {
    cat("⚠️  File not found:", input_file, "\n")
    next
  }

  tryCatch({
    sce <- read10xCounts(input_file)

    # 1. Run EmptyDrops
    cat("📉 Running EmptyDrops...\n")
    ed <- emptyDrops(counts(sce))
    keep <- which(ed$FDR < 0.01)
    sce_filtered <- sce[, keep]

    # 2. Run SoupX
    cat("🧪 Running SoupX...\n")
    sc <- SoupChannel(tod = counts(sce), toc = counts(sce_filtered))
    sc <- autoEstCont(sc)
    cleaned_counts <- adjustCounts(sc)

    # 3. Save
    dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
    DropletUtils::write10xCounts(out_path, cleaned_counts, version = "3")
    cat("✅ Saved:", out_path, "\n")

  }, error = function(e) {
    cat("❌ Error processing", sample_id, ":", e$message, "\n")
  })
}