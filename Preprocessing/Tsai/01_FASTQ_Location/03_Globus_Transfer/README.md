# Globus FASTQ Transfer to Openmind

Transfer FASTQ files from Engaging cluster to Openmind using Globus CLI batch transfer.

## Overview

- **Source:** Engaging cluster (`/nfs/picower001/lhtsailab/...`)
- **Destination:** Openmind scratch (`/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs`)
- **Data:** ~5200 FASTQ files, ~9 TB, 480 patients

## Directory Structure

Files are organized for easy Cellranger processing:

```
/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/
├── <projid>/
│   └── <Library_ID>/
│       ├── [run_id/]           # Only for samples with multiple sequencing runs
│       │   └── *.fastq.gz
│       └── *.fastq.gz          # For samples with single sequencing run
```

For samples sequenced across multiple flow cells, a `run_id` subdirectory 
(e.g., `10x-4182G`) preserves all FASTQ files.

## Scripts

### 1. generate_globus_batch.py

Generates the Globus batch transfer file from the FASTQ CSV.

```bash
# Dry run - show summary without creating files
python Scripts/generate_globus_batch.py --dry-run

# Generate batch file
python Scripts/generate_globus_batch.py
```

Output: `Batch_Files/globus_batch.txt`

### 2. submit_globus_transfer.sh

Submits the batch transfer via Globus CLI.

```bash
# Check-only mode - verify login and endpoints
./Scripts/submit_globus_transfer.sh --check-only

# Submit the transfer
./Scripts/submit_globus_transfer.sh
```

### 3. verify_transfer.py

Verifies transfer completion on Openmind.

```bash
# Generate manifest of expected files
python Scripts/verify_transfer.py --generate-manifest

# Verify transfer (run on Openmind after transfer completes)
python Scripts/verify_transfer.py

# Verbose output
python Scripts/verify_transfer.py -v
```

## Usage Steps

1. **Generate batch file:**
   ```bash
   cd 03_Globus_Transfer
   python Scripts/generate_globus_batch.py
   ```

2. **Review the batch file:**
   ```bash
   head Batch_Files/globus_batch.txt
   wc -l Batch_Files/globus_batch.txt
   ```

3. **Submit transfer:**
   ```bash
   ./Scripts/submit_globus_transfer.sh
   ```

4. **Monitor transfer:**
   - Web UI: https://app.globus.org/activity
   - CLI: `globus task show <TASK_ID>`

5. **Verify completion (on Openmind):**
   ```bash
   python Scripts/verify_transfer.py
   ```

## Globus Endpoints

| Cluster   | Endpoint ID                          | Name |
|-----------|--------------------------------------|------|
| Engaging  | `ec54b570-cac5-47f7-b2a1-100c2078686f` | MIT ORCD Engaging Collection |
| Openmind  | `cbc6f8da-d37e-11eb-bde9-5111456017d9` | mithpc#openmind |

**Note:** The older `mithpc#engaging` endpoint (`c52fcff2-761c-11eb-8cfc-cd623f92e1c0`) may have connectivity issues. Use the MIT ORCD Engaging Collection instead.

## Cellranger Integration

After transfer, point Cellranger to the patient directory:

```bash
cellranger count \
    --fastqs=/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/<projid>/<Library_ID> \
    --sample=<Library_ID> \
    ...
```

For multi-run samples, Cellranger will automatically find all FASTQs in subdirectories.

## Notes

- The `/om/scratch/Mon/` directory may have retention policies - check cluster documentation
- Globus uses checksums to verify transfer integrity
- The conda environment for Globus CLI: `/home/mabdel03/conda_envs/globus_env`
