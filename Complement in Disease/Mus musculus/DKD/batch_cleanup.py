import os
import subprocess
import time

# === SETTINGS ===
input_dir = '/Users/aumchampaneri/PycharmProjects/Biochemistry-BS/Complement in Disease/Mus musculus/DKD/raw_data'
output_dir = '/Users/aumchampaneri/PycharmProjects/Biochemistry-BS/Complement in Disease/Mus musculus/DKD/cleaned_data'
os.makedirs(output_dir, exist_ok=True)

# Define the sample to process
sample_id = "GSM5594468"
file_name = "GSM5594468_E3019_raw_feature_bc_matrix.h5"

# === FUNCTION TO RUN CELLBENDER ===
def run_cellbender(sample_id, file_name, use_cuda=False):
    input_file = os.path.join(input_dir, file_name)
    output_file = os.path.join(output_dir, f"{sample_id}_cleaned.h5")

    print(f"\n🔬 Processing {sample_id} with CellBender")

    if not os.path.exists(input_file):
        print(f"⚠️  File not found: {input_file}")
        return

    try:
        # Measure runtime
        start_time = time.time()

        # Build the CellBender command
        command = [
            "cellbender", "remove-background",
            "--input", input_file,
            "--output", output_file,
            # "--expected-cells", "3000",
            # "--total-droplets-included", "20000",
            "--zdim", "0"  # Disable checkpoint saving
        ]
        if use_cuda:
            command.append("--cuda")  # Add the --cuda flag only if GPU is enabled

        print(f"📜 Running command: {' '.join(command)}")
        subprocess.run(command, check=True)

        # Measure runtime
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"✅ Successfully processed {sample_id} with CellBender in {elapsed_time:.2f} seconds")

    except subprocess.CalledProcessError as e:
        print(f"❌ Error processing {sample_id}: {str(e)}")

# === RUN THE SCRIPT FOR ONE SAMPLE ===
run_cellbender(sample_id, file_name, use_cuda=False)  # Set use_cuda=True if you want to use a GPU