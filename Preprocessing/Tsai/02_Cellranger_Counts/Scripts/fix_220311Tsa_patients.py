#!/usr/bin/env python3
"""
Fix script for patients using the problematic 220311Tsa directory.

Creates symlink directories and regenerates Cell Ranger scripts.
"""

import os
import pandas as pd
from pathlib import Path

# Configuration
REPO_ROOT = "/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/ROSMAP-SingleNucleusRNAseq"
PIPELINE_DIR = f"{REPO_ROOT}/Preprocessing/Tsai/02_Cellranger_Counts"
INPUT_CSV = f"{REPO_ROOT}/Data/Tsai/All_ROSMAP_FASTQs.csv"

# Symlinks directory - create in permanent storage
SYMLINKS_DIR = f"{REPO_ROOT}/Data/Tsai/FASTQ_Symlinks"

# Problematic directory
PROBLEMATIC_DIR = "/nfs/picower001/lhtsailab/lab_shared/lpantano/data/220311Tsa"

# Cell Ranger settings
CELLRANGER_PATH = "/orcd/data/lhtsai/001/om2/mabdel03/apps/yard/cellranger-8.0.0"
TRANSCRIPTOME = "/orcd/data/lhtsai/001/om2/mabdel03/yard/references/human/refdata-gex-GRCh38-2020-A"
SCRATCH_ROOT = "/home/mabdel03/orcd/scratch"
CELLRANGER_OUTPUT = f"{SCRATCH_ROOT}/Tsai/Cellranger_Counts"
TRACKING_DIR = f"{PIPELINE_DIR}/Tracking"
LOGS_OUT = f"{PIPELINE_DIR}/Logs/Outs"
LOGS_ERR = f"{PIPELINE_DIR}/Logs/Errs"


def create_symlink_directory(projid, library_id, fastq_files):
    """Create a symlink directory for a patient's FASTQs."""
    symlink_dir = Path(SYMLINKS_DIR) / library_id
    symlink_dir.mkdir(parents=True, exist_ok=True)
    
    created = 0
    for fastq_path in fastq_files:
        fastq_path = Path(fastq_path)
        if fastq_path.exists():
            symlink_path = symlink_dir / fastq_path.name
            if symlink_path.exists():
                symlink_path.unlink()
            symlink_path.symlink_to(fastq_path)
            created += 1
    
    return str(symlink_dir), created


def generate_cellranger_script(projid, library_id, fastq_dir, batch_num):
    """Generate a Cell Ranger batch script."""
    output_dir = f"{CELLRANGER_OUTPUT}/{projid}"
    
    script = f"""#!/bin/bash
#SBATCH --job-name=cr_{projid}
#SBATCH --partition=mit_preemptable
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output={LOGS_OUT}/cellranger_{projid}_%j.out
#SBATCH --error={LOGS_ERR}/cellranger_{projid}_%j.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=FAIL

# Cell Ranger count for patient {projid} (Batch {batch_num})
# Library ID: {library_id}
# Using symlink directory (fixed)

set -e

# Add Cell Ranger to PATH
export PATH={CELLRANGER_PATH}:$PATH

# Create output directory
mkdir -p {output_dir}

# Run Cell Ranger count
cellranger count \\
    --create-bam=false \\
    --include-introns=true \\
    --nosecondary \\
    --r1-length=26 \\
    --id={library_id} \\
    --transcriptome={TRANSCRIPTOME} \\
    --sample={library_id} \\
    --fastqs={fastq_dir} \\
    --output-dir={output_dir}

# Mark completion
echo "{projid}" >> {TRACKING_DIR}/cellranger_completed.txt

echo "Cell Ranger completed for {projid}"
"""
    return script


def main():
    print("=" * 60)
    print("Fixing 220311Tsa patients")
    print("=" * 60)
    print()
    
    # Read the master CSV
    df = pd.read_csv(INPUT_CSV)
    print(f"CSV columns: {list(df.columns)}")
    print()
    
    # Find patients using the problematic directory (using correct column name)
    affected = df[df['source_dir'].str.contains('220311Tsa', na=False)]
    affected_projids = affected['projid'].unique()
    
    print(f"Found {len(affected_projids)} affected patients")
    print()
    
    # Read batch assignments
    batch_assignments = pd.read_csv(f"{TRACKING_DIR}/batch_assignments.csv")
    
    # Create symlinks directory
    Path(SYMLINKS_DIR).mkdir(parents=True, exist_ok=True)
    
    fixed_scripts = []
    
    for projid in affected_projids:
        patient_data = affected[affected['projid'] == projid]
        library_id = patient_data['Library_ID'].iloc[0]
        
        # Get all FASTQ files for this patient from 220311Tsa
        fastq_files = patient_data[patient_data['source_dir'].str.contains('220311Tsa', na=False)]['full_path'].tolist()
        
        # Filter to only R1 and R2 files (Cell Ranger needs these)
        r1_r2_files = [f for f in fastq_files if '_R1_' in f or '_R2_' in f]
        
        print(f"Patient {projid} ({library_id}):")
        print(f"  FASTQ files: {len(r1_r2_files)}")
        
        # Create symlink directory
        symlink_dir, created = create_symlink_directory(projid, library_id, r1_r2_files)
        print(f"  Symlinks created: {created} in {symlink_dir}")
        
        # Get batch number
        batch_info = batch_assignments[batch_assignments['projid'] == projid]
        if len(batch_info) > 0:
            batch_num = batch_info['batch'].iloc[0]
        else:
            batch_num = 1  # Default
        
        # Generate new script
        script_content = generate_cellranger_script(projid, library_id, symlink_dir, batch_num)
        
        # Write to batch scripts directory
        script_path = f"{PIPELINE_DIR}/Batch_Scripts/batch_{batch_num}/cellranger/{projid}_cellranger.sh"
        with open(script_path, 'w') as f:
            f.write(script_content)
        os.chmod(script_path, 0o755)
        
        print(f"  Script updated: {script_path}")
        fixed_scripts.append((projid, library_id, batch_num, script_path))
        print()
    
    print("=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"Created symlink directories for {len(affected_projids)} patients")
    print(f"Updated {len(fixed_scripts)} Cell Ranger scripts")
    print()
    print("Fixed patients by batch:")
    
    # Group by batch
    from collections import defaultdict
    by_batch = defaultdict(list)
    for projid, library_id, batch_num, script_path in fixed_scripts:
        by_batch[batch_num].append((projid, script_path))
    
    for batch_num in sorted(by_batch.keys()):
        patients = by_batch[batch_num]
        print(f"  Batch {batch_num}: {len(patients)} patients")
        for projid, script_path in patients:
            print(f"    - {projid}")
    
    # Write list of scripts to resubmit
    resubmit_file = f"{PIPELINE_DIR}/Scripts/resubmit_fixed.sh"
    with open(resubmit_file, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("# Resubmit fixed Cell Ranger jobs\n\n")
        for projid, library_id, batch_num, script_path in fixed_scripts:
            f.write(f"echo 'Submitting {projid}...'\n")
            f.write(f"sbatch {script_path}\n")
        f.write("\necho 'All fixed jobs submitted!'\n")
    os.chmod(resubmit_file, 0o755)
    
    print()
    print(f"Resubmit script created: {resubmit_file}")
    print()
    print("To resubmit failed jobs, run:")
    print(f"  {resubmit_file}")


if __name__ == "__main__":
    main()
