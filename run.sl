#!/bin/bash
#SBATCH --job-name=MCMC_E_L32_b6.0    # Job name: MCMC, Energy, L=32, beta=6.0
#SBATCH --output=outputs/%j/job_%j.out   # Standard output log
#SBATCH --error=outputs/%j/job_%j.err    # Standard error log
#SBATCH --nodes=1                      # Run on a single node
#SBATCH --ntasks=1                     # Run a single task (our script is serial)
#SBATCH --cpus-per-task=1              # Only needs one CPU
#SBATCH --partition=medium              # Or 'medium'/'long' - 'short' is too short
#SBATCH --mem=16G                       # 4GB of memory should be plenty
#SBATCH --time=24:00:00                # 24-hour run time

#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=rappapo@bc.edu     # Copied from your template

#======================================================================
# Part 1: Experiment Hyperparameters 🔬
#======================================================================
DESCRIPTION="Overnight L=32 plaquette model run (p=0.1, beta=6.0) to find energy plateaus."

# Define all your MCMC hyperparameters here
L_SIZE=32
PROB_P=0.1
BETA=6.0
T_MAX=1000000.0
SNAPSHOT_INTERVAL=100.0
NUM_RUNS=20

#======================================================================
# Part 2: Setup and Reproducibility Logging ✍️
#======================================================================

# Create the unique output directory
OUTPUT_DIR="outputs/${SLURM_JOB_ID}"
mkdir -p $OUTPUT_DIR

# 1. Copy the submission script itself into the output directory.
cp $0 $OUTPUT_DIR/submitted_script.sh

# 2. Create a log file for versioning information.
LOG_FILE="$OUTPUT_DIR/reproducibility.log"
echo "Reproducibility log for job $SLURM_JOB_ID" > $LOG_FILE
echo "Description: $DESCRIPTION" >> $LOG_FILE
echo "Submission time: $(date)" >> $LOG_FILE
echo "Running on node: $(hostname)" >> $LOG_FILE
echo "-------------------------------------------" >> $LOG_FILE

# 3. Log the Git commit hash (assuming you're in a git repo).
echo "Git Commit Hash:" >> $LOG_FILE
git rev-parse HEAD >> $LOG_FILE 2>&1
echo "" >> $LOG_FILE

# 4. Log any uncommitted changes.
echo "Uncommitted Changes (git status):" >> $LOG_FILE
git status --porcelain >> $LOG_FILE 2>&1
echo "" >> $LOG_FILE

# 5. Save the exact uncommitted changes as a diff patch.
git diff > $OUTPUT_DIR/uncommitted_changes.patch

# 6. Define the final output CSV path
#    This ensures the CSV is saved *inside* our new job folder
OUTPUT_FILE_PATH="$OUTPUT_DIR/energy_L${L_SIZE}_p${PROB_P}_b${BETA}.csv"

#======================================================================
# Part 3: Main Execution Block 🚀
#======================================================================

START_TIME=$SECONDS

# Navigate to your MCMC project directory
# (Based on your previous error message)
cd /home/rappapo/plaquetteMCMC

# Activate the Julia environment
module purge
module load julia  # Or your cluster's specific Julia module

echo "Starting Julia simulation..."

# Run the Julia script with all hyperparameters
julia plaquetteham.jl \
    --L $L_SIZE \
    --p $PROB_P \
    --beta $BETA \
    --T_max $T_MAX \
    --snapshot_interval $SNAPSHOT_INTERVAL \
    --num_runs $NUM_RUNS \
    --output_file $OUTPUT_FILE_PATH

END_TIME=$SECONDS
DURATION=$((END_TIME - START_TIME))

# Log the total runtime to your reproducibility file
echo "" >> $LOG_FILE
echo "-------------------------------------------" >> $LOG_FILE
echo "Total wall-clock runtime: ${DURATION} seconds" >> $LOG_FILE
echo "Which is $(($DURATION / 60)) minutes and $(($DURATION % 60)) seconds." >> $LOG_FILE

echo "Simulation finished. Results saved to $OUTPUT_DIR"