# Config

This directory contains configuration files for the FASTQ Transfer Pipeline.

## Files

### `config.sh`
Central configuration file that defines all paths and parameters used by the pipeline scripts.

## Configuration Variables

### Conda Configuration
| Variable | Description |
|----------|-------------|
| `CONDA_INIT_SCRIPT` | Path to conda initialization script (miniforge3) |
| `CONDA_ENV` | Path to conda environment with pandas |

### Base Paths
| Variable | Description |
|----------|-------------|
| `PIPELINE_ROOT` | Root directory for the entire pipeline |
| `OLD_DIR` | Directory containing original CSVs and scripts |
| `NEW_DIR` | Directory for new outputs |

### Input CSVs
| Variable | Description |
|----------|-------------|
| `MASTER_CSV` | Path to `Tsai_To_Openmind.csv` (projid -> Library_ID -> path mapping) |
| `BACKUPS_CSV` | Path to `Reformatted_Backups.csv` (backup/fallback paths) |

### Output Directories
| Variable | Description |
|----------|-------------|
| `CSV_OUTPUT_DIR` | Where generated CSVs are saved |
| `FASTQ_OUTPUT_DIR` | Where organized FASTQ symlinks are created |
| `LOG_DIR` | Where SLURM logs are written |
| `TEMP_DIR` | Temporary files during processing |

### Output Files
| Variable | Description |
|----------|-------------|
| `ALL_FASTQS_CSV` | Master CSV with all FASTQ file paths |
| `VALIDATION_REPORT` | Validation report file |

### Parallelization Parameters
| Variable | Default | Description |
|----------|---------|-------------|
| `NUM_PARALLEL_JOBS` | 32 | Number of GNU parallel workers |
| `SLURM_TIME` | 10:00:00 | SLURM job time limit |
| `SLURM_MEM` | 64G | SLURM memory allocation |
| `SLURM_CPUS` | 32 | SLURM CPU allocation |
| `SLURM_MAIL_USER` | mabdel03@mit.edu | Email for SLURM notifications |

## Usage

The config file is sourced by other scripts:

```bash
source /path/to/Config/config.sh
```

All variables are exported and available to child processes.

## Helper Functions

### `ensure_dirs()`
Creates all necessary output directories if they don't exist.

```bash
source config.sh
ensure_dirs
```

## Customization

Edit `config.sh` to change:
- File paths (if moved to different location)
- SLURM parameters (time, memory, CPUs)
- Parallelization settings
- Email notifications

