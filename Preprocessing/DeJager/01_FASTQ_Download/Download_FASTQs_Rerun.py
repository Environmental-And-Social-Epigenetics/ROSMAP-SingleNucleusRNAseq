"""
Re-download failed DeJager FASTQs from Synapse.

Checks each file against what already exists on disk and only downloads
files that are missing.  Uses stored Synapse credentials (never hardcode
auth tokens in source code).

Usage:
    source config/paths.sh
    python Download_FASTQs_Rerun.py
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

# Login via stored credentials or SYNAPSE_AUTH_TOKEN env var
syn = synapseclient.Synapse()
syn.login()

# Iterate over the rows of the DataFrame
for row in df.index:
    library_id = df["LibraryID"][row]
    syn_id = df["SynID"][row]
    file_name = df["FileName"][row]
    download_dir = os.path.join(root, library_id)

    # Ensure the download directory exists
    os.makedirs(download_dir, exist_ok=True)

    # Check if the specific file is already present
    downloaded_files = os.listdir(download_dir)
    if file_name not in downloaded_files:
        try:
            entity = syn.get(syn_id, downloadLocation=download_dir)
            print(f"Downloaded {file_name} (SynID {syn_id}) for Library ID {library_id}")
        except Exception as e:
            print(f"Exception downloading {file_name} (SynID {syn_id}) for Library ID {library_id}: {e}")
    else:
        print(f"{file_name} already downloaded for Library ID {library_id}")
