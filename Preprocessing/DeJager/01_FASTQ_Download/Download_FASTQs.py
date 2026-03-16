"""
Download DeJager FASTQs from Synapse to the local filesystem.

Organizes FASTQs into directories by LibraryID.

Prerequisites:
    - Synapse credentials configured via `synapse login --rememberMe`
      or the SYNAPSE_AUTH_TOKEN environment variable.
    - CSV with columns: SynID, LibraryID
    - DEJAGER_FASTQS path set (via config/paths.sh or environment)

Usage:
    source config/paths.sh
    python Download_FASTQs.py
"""

import os
from pathlib import Path

import pandas as pd
import synapseclient


# Compute repo/workspace roots from script location (depth 3 from repo root)
_script_path = Path(__file__).resolve()
_repo_root = _script_path.parents[3]
_workspace_root = _repo_root.parent

# Resolve paths from environment (set by config/paths.sh) or use defaults
DATA_ROOT = os.environ.get(
    "DATA_ROOT",
    "__UNCONFIGURED__set_DATA_ROOT_in_paths_local_sh",
)
DEJAGER_FASTQS = os.environ.get(
    "DEJAGER_FASTQS",
    str(_workspace_root / "DeJager_Data" / "FASTQs"),
)

SYNAPSE_CSV = os.path.join(
    DATA_ROOT, "Data/DeJager/FASTQs_Download/FASTQ_Download_CSVs/Synapse_FASTQ_IDs.csv"
)

df = pd.read_csv(SYNAPSE_CSV)
root = DEJAGER_FASTQS

# Create library directories
for libID in set(df["LibraryID"]):
    os.makedirs(os.path.join(root, libID), exist_ok=True)

# Login via stored credentials or SYNAPSE_AUTH_TOKEN env var
# NEVER hardcode auth tokens in source code
syn = synapseclient.Synapse()
syn.login()

for row in df.index:
    entity = syn.get(
        df["SynID"][row], downloadLocation=os.path.join(root, df["LibraryID"][row])
    )
