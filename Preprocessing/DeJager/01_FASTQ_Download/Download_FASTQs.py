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

import pandas as pd
import synapseclient


# Resolve paths from environment (set by config/paths.sh) or use defaults
DATA_ROOT = os.environ.get(
    "DATA_ROOT",
    "/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis",
)
DEJAGER_FASTQS = os.environ.get(
    "DEJAGER_FASTQS",
    os.path.join(os.environ.get("SCRATCH_ROOT", "/om/scratch/Mon/mabdel03"), "FASTQs"),
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
