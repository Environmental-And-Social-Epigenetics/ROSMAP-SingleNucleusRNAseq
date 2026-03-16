# Tsai Server (TsaiLabNAS) Transfers

SFTP-based upload/download scripts for transferring transcriptomics data to/from the
Tsai Lab NAS (Synology, DSM 7.3.2).

## Connection Details

- **Host**: `tsailabnas.mit.edu`
- **Port**: 22 (SFTP)
- **User**: `mabdel03`
- **Auth**: Password file at `~/.smb_tsailabnas`
- **NAS chroot**: `/LabShared` (SFTP root is `/`, which maps to `/LabShared` on the NAS)

## Usage

### Upload (single worker)

```bash
cd Data_Access/Transcriptomics/Tsai_Server/Tsai
bash upload_to_nas.sh upload
```

### Upload (parallel via SLURM, recommended for large datasets)

```bash
sbatch slurm_upload.sh
```

### Download

```bash
bash download_from_nas.sh
```

### Verify uploads

```bash
bash upload_to_nas.sh verify
```

## Configuration

Each cohort directory has a `.env` file controlling paths and SFTP settings.
Edit the `.env` to change data types, NAS paths, or SFTP tuning parameters.

## SFTP Tuning

- `SFTP_NUM_REQUESTS=64` - number of outstanding requests (higher = faster on high-latency links)
- `SFTP_CIPHER=aes128-gcm@openssh.com` - fast cipher
- `SFTP_BUFFER_SIZE=""` - leave empty for default 32KB (larger values cause "Outbound message too long" on this NAS)

## Known Issues

- **"Outbound message too long"**: Do NOT increase SFTP_BUFFER_SIZE above 32KB. The Synology NAS has a max packet size limit.
- **"Couldn't create directory: Failure"**: Expected for `mkdir` on existing directories. The `-` prefix in batch files makes this non-fatal.
- **WIN\Domain Users**: The NAS reports file ownership with spaces in the group name. The manifest parser handles this dynamically.
