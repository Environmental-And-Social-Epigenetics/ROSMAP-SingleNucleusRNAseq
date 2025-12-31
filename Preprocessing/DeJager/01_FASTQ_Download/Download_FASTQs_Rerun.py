"""
Script to download all of the DeJager FASTQs to Openmind
Organizes FASTQs into directories by LibraryID
NOT PARTITIONED --> DOWNLOAD ALL FASTQs 
"""

import synapseclient
# import os
import pandas as pd
from synapseclient import Project, File, Folder


# df = pd.read_csv('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/DeJager/FASTQs_Download/FASTQ_Download_CSVs/Synapse_FASTQ_IDs.csv')
# root = '/om/scratch/Sun/mabdel03/ROSMAP_SC/DeJager/Preprocessing/FASTQs'

# syn = synapseclient.login(authToken="eyJ0eXAiOiJKV1QiLCJraWQiOiJXN05OOldMSlQ6SjVSSzpMN1RMOlQ3TDc6M1ZYNjpKRU9VOjY0NFI6VTNJWDo1S1oyOjdaQ0s6RlBUSCIsImFsZyI6IlJTMjU2In0.eyJhY2Nlc3MiOnsic2NvcGUiOlsidmlldyIsImRvd25sb2FkIiwibW9kaWZ5Il0sIm9pZGNfY2xhaW1zIjp7fX0sInRva2VuX3R5cGUiOiJQRVJTT05BTF9BQ0NFU1NfVE9LRU4iLCJpc3MiOiJodHRwczovL3JlcG8tcHJvZC5wcm9kLnNhZ2ViYXNlLm9yZy9hdXRoL3YxIiwiYXVkIjoiMCIsIm5iZiI6MTcyMzc1MTgyNiwiaWF0IjoxNzIzNzUxODI2LCJqdGkiOiIxMTEyNCIsInN1YiI6IjM0NzM0MDkifQ.GLH2upv_u-0tgXKJUB7y1zilPH8YOwQENqvzhZbrU6tuRmrVsmLBfPRCOxnoc34-C7tltg3LG3Yh8SOnGVNbrwnnPfc9uUW4kNhyUdp9ZExp4qyPaxS1otB1FgWL9745vQnTiQqNQJwSgSiCb25xCMnFhKwg29TVlMRuYfSS3lNeQWNWiwfy_C9KP2DA26JjCg8d7VbQbp_4rVZv90YK7LJv00hIzGjAVO_jC0ufDTteE1CXZYEdiPLu-wyV1im7qXywSHwmOLsFZa7jGgr1srnBfCZo4ZNX-iklOZajW0P-Gfl_HsySpi1E_wPqwFqlWKhLzQ2dXCdRSFJs8aaqJw")

# for row in df.index:
# 	if not os.listdir(os.path.join(root, df['LibraryID'][row])):
# 		try:
# 			entity = syn.get(df['SynID'][row], downloadLocation=os.path.join(root, df['LibraryID'][row]))
# 		except:
# 			print(f'exception downloading synID {df['SynID'][row]} for library ID {df['LibraryID'][row]}')
# 	else:
# 		print(f'{df['LibraryID'][row]} already downloaded')


import os
from synapseclient import Synapse


df = pd.read_csv('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/DeJager/FASTQs_Download/FASTQ_Download_CSVs/Synapse_FASTQ_IDs.csv')
root = '/om/scratch/Mon/mabdel03/FASTQs'


syn = synapseclient.Synapse()
syn.login(authToken="eyJ0eXAiOiJKV1QiLCJraWQiOiJXN05OOldMSlQ6SjVSSzpMN1RMOlQ3TDc6M1ZYNjpKRU9VOjY0NFI6VTNJWDo1S1oyOjdaQ0s6RlBUSCIsImFsZyI6IlJTMjU2In0.eyJhY2Nlc3MiOnsic2NvcGUiOlsidmlldyIsImRvd25sb2FkIiwibW9kaWZ5Il0sIm9pZGNfY2xhaW1zIjp7fX0sInRva2VuX3R5cGUiOiJQRVJTT05BTF9BQ0NFU1NfVE9LRU4iLCJpc3MiOiJodHRwczovL3JlcG8tcHJvZC5wcm9kLnNhZ2ViYXNlLm9yZy9hdXRoL3YxIiwiYXVkIjoiMCIsIm5iZiI6MTc1Mzc3MTc0OSwiaWF0IjoxNzUzNzcxNzQ5LCJqdGkiOiIyMzY0OSIsInN1YiI6IjM0OTYxNTMifQ.fmkXrw-eB2dQh7l8ZuwRfLSWZKaYePiaBjKPS_ZRJNv9YUGV3ApwVO43gKIhV_l7cCB9l8-GSgOEz_l-YerrlYknaJje4M1FddZxOxxXO-OeHPtEMyFurqX5eHxs9pqkJKREZ108fHKekY4RVtdirQargfkKqC8ktyo2pHhTyzSiq7R1gRq9qQZir38ITTaFShxrjo3w3u_-68bXoUFFNqayEraXOjPqoTx9xTNVDCE9XjDfs-Zvi6g1AL5WkVosE35HqUjSlWvA7OpLJ4kxNXbb6_ZBIY8624XG5K6_fRXe9MlyC0wze-J3s7vJ39G4LVSTmZrdUp9KI9hLG3f69w")

# Iterate over the rows of the DataFrame
for row in df.index:
    library_id = df['LibraryID'][row]
    syn_id = df['SynID'][row]
    file_name = df['FileName'][row]  # Get the file name from the DataFrame
    download_dir = os.path.join(root, library_id)

    # Ensure the download directory exists
    if not os.path.exists(download_dir):
        os.makedirs(download_dir)

    # List all files currently in the download directory
    downloaded_files = os.listdir(download_dir)

    # Check if the specific file is present
    if file_name not in downloaded_files:
        try:
            # Download the file from Synapse
            entity = syn.get(syn_id, downloadLocation=download_dir)
            print(f'Downloaded {file_name} (SynID {syn_id}) for Library ID {library_id}')
        except Exception as e:
            print(f'Exception downloading {file_name} (SynID {syn_id}) for Library ID {library_id}: {e}')
    else:
        print(f'{file_name} already downloaded for Library ID {library_id}')
