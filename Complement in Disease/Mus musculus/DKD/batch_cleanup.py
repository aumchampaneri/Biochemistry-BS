import os
import subprocess
import time

# === SETTINGS ===
input_dir = "/Users/aumchampaneri/PycharmProjects/Biochemistry-BS/Complement in Disease/Mus musculus/DKD/raw_data"
output_dir = "/Users/aumchampaneri/PycharmProjects/Biochemistry-BS/Complement in Disease/Mus musculus/DKD/cleaned_data"
os.makedirs(output_dir, exist_ok=True)

# Define the sample to process
sample_id = "GSM5594468"
file_name = "GSM5594468_E3019_raw_feature_bc_matrix.h5"

# === FUNCTION TO RUN CELLBENDER ===
def run_cellbender(sample_id, file_name):
    input_file = os.path.join(input_dir, file_name)
    output_file = os.path.join(output_dir, f"{sample_id}_cleaned.h5")

    print(f"\n🔬 Processing {sample_id} with CellBender")

    if not os.path.exists(input_file):
        print(f"⚠️  File not found: {input_file}")
        return

    try:
        # Measure runtime
        start_time = time.time()

        # Run CellBender
        command = [
            "cellbender", "remove-background",
            "--input", input_file,
            "--output", output_file,
            "--cuda", "False",  # Set to "True" if you have a GPU
            "--expected-cells", "3000",  # Adjust based on your dataset
            "--total-droplets-included", "20000"  # Adjust based on your dataset
        ]
        print(f"📜 Running command: {' '.join(command)}")
        subprocess.run(command, check=True)

        # Measure runtime
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"✅ Successfully processed {sample_id} with CellBender in {elapsed_time:.2f} seconds")

    except subprocess.CalledProcessError as e:
        print(f"❌ Error processing {sample_id}: {str(e)}")

# === RUN THE SCRIPT FOR ONE SAMPLE ===
run_cellbender(sample_id, file_name)