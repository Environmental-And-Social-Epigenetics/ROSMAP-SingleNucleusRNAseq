# Engaging-Openmind Globus Transfers

Scripts for transferring transcriptomics data between MIT HPC clusters using Globus.

## Globus Overview

Globus is the recommended method for large inter-cluster transfers. It handles
retries, checkpointing, and parallelism automatically.

**Two ways to use:**
1. **GUI**: https://app.globus.org/ (point-and-click, good for one-off transfers)
2. **CLI**: `globus transfer` command (scriptable, used by these scripts)

## Endpoint IDs

| Cluster   | Endpoint Name   | UUID                                   | Status |
|-----------|----------------|----------------------------------------|--------|
| Openmind  | mithpc#openmind | `cbc6f8da-d37e-11eb-bde9-5111456017d9` | Decommissioned March 2026 |
| Engaging  | MIT ORCD Engaging Collection | `ec54b570-cac5-47f7-b2a1-100c2078686f` | Active |

To find endpoint UUIDs: `globus endpoint search <NAME>`

Make sure UUIDs match the client email you authenticated with.

## Conda Environments

- **Openmind**: `/om2/user/mabdel03/conda_envs/globus_env`
- **Engaging**: `/home/mabdel03/conda_envs/globus_env`

## Usage

### Send data to another cluster

```bash
# Generate batch file listing all folders to transfer
bash Tsai/generate_batch.sh

# Submit the transfer
bash Tsai/send_globus.sh
```

### Receive data from another cluster

```bash
bash Tsai/receive_globus.sh
```

### Check transfer status

```bash
conda activate /om2/user/mabdel03/conda_envs/globus_env
globus task show <TASK_ID>
globus task wait <TASK_ID>
```

## Batch Transfers

For transferring many folders, use batch mode:

1. Generate a batch file (one line per folder):
   ```
   /source/path/folder1 /dest/path/folder1 --recursive
   /source/path/folder2 /dest/path/folder2 --recursive
   ```

2. Submit:
   ```bash
   globus transfer $SOURCE_ENDPOINT $DEST_ENDPOINT --batch batch_file.txt
   ```

The `generate_batch.sh` scripts automate step 1 for each cohort's data.
