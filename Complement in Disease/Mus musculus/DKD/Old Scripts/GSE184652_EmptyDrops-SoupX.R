library(DropletUtils)
library(SingleCellExperiment)
library(SoupX)

# === SETTINGS ===
input_dir <- '/Users/aumchampaneri/Databases/Mm DKD Two/GSE184652_RAW'
output_dir <- '/Users/aumchampaneri/Databases/Mm DKD Two/GSE184652_EmptyDrops-SoupX'

# Map of sample IDs to actual file names
sample_files <- c(
  GSM5594468 = "GSM5594468_E3019_raw_feature_bc_matrix.h5",
  GSM5594469 = "GSM5594469_A3020_raw_feature_bc_matrix.h5",
  GSM5594470 = "GSM5594470_F3021_raw_feature_bc_matrix.h5",
  GSM5594471 = "GSM5594471_G3022_raw_feature_bc_matrix.h5",
  GSM5594472 = "GSM5594472_H3023_raw_feature_bc_matrix.h5",
  GSM5594473 = "GSM5594473_N3024_raw_feature_bc_matrix.h5",
  GSM5594474 = "GSM5594474_B3025_raw_feature_bc_matrix.h5",
  GSM5594475 = "GSM5594475_A3026_raw_feature_bc_matrix.h5",
  GSM5594476 = "GSM5594476_B3027_raw_feature_bc_matrix.h5",
  GSM5594477 = "GSM5594477_C3028_raw_feature_bc_matrix.h5",
  GSM5594478 = "GSM5594478_A3029_raw_feature_bc_matrix.h5"
)

# === PROCESS EACH SAMPLE ===
for (sample_id in names(sample_files)) {
  cat("🔬 Processing", sample_id, "\n")

  input_path <- file.path(input_dir, sample_files[[sample_id]])
  sce <- read10xCounts(input_path)

  # 1. EmptyDrops
  cat("📉 Running EmptyDrops...\n")
  ed <- emptyDrops(counts(sce))
  keep <- which(ed$FDR < 0.01)
  sce_filtered <- sce[, keep]

  # 2. SoupX
  cat("🧪 Running SoupX ambient RNA removal...\n")
  sc <- SoupChannel(tod = counts(sce), toc = counts(sce_filtered))
  sc <- autoEstCont(sc)
  cleaned_counts <- adjustCounts(sc)

  # 3. Save cleaned matrix
  out_dir <- file.path(output_dir, sample_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  DropletUtils::write10xCounts(out_dir, cleaned_counts, version = "3")

  cat("✅ Saved cleaned matrix to", out_dir, "\n\n")
}