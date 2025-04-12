import os
import scanpy as sc
import anndata
import numpy as np
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector

# Enable R-Python data conversion
pandas2ri.activate()

# Import R packages
DropletUtils = importr('DropletUtils')
SoupX = importr('SoupX')

# === SETTINGS ===
input_dir = '/workspaces/Biochemistry-BS/Complement in Disease/Mus musculus/DKD/raw_data'
output_dir = '/workspaces/Biochemistry-BS/Complement in Disease/Mus musculus/DKD/cleaned_data'
os.makedirs(output_dir, exist_ok=True)

samples = {
    "GSM5594468": "GSM5594468_E3019_raw_feature_bc_matrix.h5",
    "GSM5594469": "GSM5594469_A3020_raw_feature_bc_matrix.h5",
    "GSM5594470": "GSM5594470_F3021_raw_feature_bc_matrix.h5",
    "GSM5594471": "GSM5594471_G3022_raw_feature_bc_matrix.h5",
    "GSM5594472": "GSM5594472_H3023_raw_feature_bc_matrix.h5",
    "GSM5594473": "GSM5594473_N3024_raw_feature_bc_matrix.h5",
    "GSM5594474": "GSM5594474_B3025_raw_feature_bc_matrix.h5",
    "GSM5594475": "GSM5594475_A3026_raw_feature_bc_matrix.h5",
    "GSM5594476": "GSM5594476_B3027_raw_feature_bc_matrix.h5",
    "GSM5594477": "GSM5594477_C3028_raw_feature_bc_matrix.h5",
    "GSM5594478": "GSM5594478_A3029_raw_feature_bc_matrix.h5"
}

# === PROCESS SAMPLES ONE AT A TIME ===
for sample_id, file_name in samples.items():
    input_file = os.path.join(input_dir, file_name)
    out_path = os.path.join(output_dir, sample_id)

    print(f"\n🔬 Processing {sample_id}")

    if not os.path.exists(input_file):
        print(f"⚠️  File not found: {input_file}")
        continue

    try:
        # 1. Read 10x data
        print("📖 Reading 10x data...")
        sce = DropletUtils.read10xCounts(input_file)

        # 2. Run EmptyDrops
        print("📉 Running EmptyDrops...")
        ed = DropletUtils.emptyDrops(sce.rx2('counts'))
        keep = np.where(np.array(ed.rx2('FDR')) < 0.01)[0]
        sce_filtered = sce.rx2('counts')[:, keep]

        # 3. Run SoupX
        print("🧪 Running SoupX...")
        sc = SoupX.SoupChannel(tod=sce.rx2('counts'), toc=sce_filtered)
        sc = SoupX.autoEstCont(sc)
        cleaned_counts = SoupX.adjustCounts(sc)

        # 4. Save
        print("💾 Saving cleaned data...")
        os.makedirs(out_path, exist_ok=True)
        DropletUtils.write10xCounts(out_path, cleaned_counts, version="3")
        print(f"✅ Saved: {out_path}")

    except Exception as e:
        print(f"❌ Error processing {sample_id}: {str(e)}")