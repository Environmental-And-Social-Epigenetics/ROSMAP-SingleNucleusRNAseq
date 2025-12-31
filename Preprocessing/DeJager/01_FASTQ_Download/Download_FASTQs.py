"""
Script to download all of the DeJager FASTQs to Openmind
Organizes FASTQs into directories by LibraryID
NOT PARTITIONED --> DOWNLOAD ALL FASTQs 
"""

import synapseclient
import os
import pandas as pd
from synapseclient import Project, File, Folder


df = pd.read_csv('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/DeJager/FASTQs_Download/FASTQ_Download_CSVs/Synapse_FASTQ_IDs.csv')
root = '/om/scratch/Mon/mabdel03/FASTQs'

for libID in set(df['LibraryID']):
    os.mkdir(os.path.join(root, libID))


syn = synapseclient.Synapse()
syn.login(authToken="eyJ0eXAiOiJKV1QiLCJraWQiOiJXN05OOldMSlQ6SjVSSzpMN1RMOlQ3TDc6M1ZYNjpKRU9VOjY0NFI6VTNJWDo1S1oyOjdaQ0s6RlBUSCIsImFsZyI6IlJTMjU2In0.eyJhY2Nlc3MiOnsic2NvcGUiOlsidmlldyIsImRvd25sb2FkIiwibW9kaWZ5Il0sIm9pZGNfY2xhaW1zIjp7fX0sInRva2VuX3R5cGUiOiJQRVJTT05BTF9BQ0NFU1NfVE9LRU4iLCJpc3MiOiJodHRwczovL3JlcG8tcHJvZC5wcm9kLnNhZ2ViYXNlLm9yZy9hdXRoL3YxIiwiYXVkIjoiMCIsIm5iZiI6MTc1MDUzMjQ2MSwiaWF0IjoxNzUwNTMyNDYxLCJqdGkiOiIyMTkyNSIsInN1YiI6IjM0OTYxNTMifQ.qkR8uCe4QnlS8GZRo89C2535wcdVEFwjCNre5q7LGWiMHT8wUAD4Fr_mIwLmeCg2HmSxl79CE9CwIs1RcZLRKG95ZB0GW4r62X9NUACZClfPPJYgwnrQiXtx2nHbfB9DtIFcyYUlpDgxBcstd92eg8ew7MJ_YrX2DnD8vlHuyZIgfWUoGbcd0L3pwPwK0X2ePxrVwiK6lY6EUJ5QU5qOfiRAKLNi8e9hwN1lXtD7EGK3MWdkmIa3NRrAn3aIZQnPAY9zevH9VkqF_Ui6jNCE2V6CDiBN6r0mE1jFcVmOcfRtnnYR7KIXpIZJpgVjYMyYqEmS1G2wrOe12-XqRj2dJg")


for row in df.index:
    entity = syn.get(df['SynID'][row], downloadLocation=os.path.join(root, df['LibraryID'][row]))
