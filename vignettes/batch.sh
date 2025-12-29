#!/bin/bash
#SBATCH -A iicd
#SBATCH -J 1172
#SBATCH --time=12:00:00            
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=knd2127@columbia.edu
#SBATCH --array=1-102%20                # one task per R-file
#   SBATCH --output=logs/%x_%A_%a.out
#   SBATCH --error=logs/%x_%A_%a.err

set -e

ACTUAL_SAMPLE_INDEX=$((SLURM_ARRAY_TASK_ID + 0))

WORK_DIR="/burg-archive/iicd/users/knd2127/DriverSelectionSweep/vignettes"
cd "$WORK_DIR"
mkdir -p logs/Fabre

module purge
module load R/4.4.2
R_SCRIPT="inference_Fabre_parallel.r"

srun Rscript "$R_SCRIPT" "$ACTUAL_SAMPLE_INDEX"

echo ""
echo "=============================================================="
if [ $? -eq 0 ]; then
    echo "Successful run"
else
    echo "Error"
fi
echo "Task ID $SLURM_ARRAY_TASK_ID"
echo "Sample index: $ACTUAL_SAMPLE_INDEX"
echo "Time: $(date)"
echo "=============================================================="
