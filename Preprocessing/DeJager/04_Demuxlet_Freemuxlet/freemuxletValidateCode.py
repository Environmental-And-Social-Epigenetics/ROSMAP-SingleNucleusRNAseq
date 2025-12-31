import pandas as pd
import gzip
import shutil

with gzip.open("output/freemux2.clust1.samples.gz", 'rb') as f_in:
    with open("output/freemux2.clust1.samples", 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

freemuxlet_df = pd.read_csv('output/freemux2.clust1.samples',sep = '\t')

annotations_df = pd.read_csv('cell-annotation.csv')

freemuxlet_df.rename(columns={'BARCODE': 'barcode'}, inplace=True)

#subset cell annotation for ths particular library (190409-B5-B)

annotations_df = annotations_df[annotations_df['barcode'].str.startswith("190409-B5-B_")]
annotations_df['barcode'] = annotations_df['barcode'].str.lstrip("190409-B5-B_")

# Assuming both dataframes have a common column 'cell_id'
merged_df = pd.merge(freemuxlet_df, annotations_df, on='barcode')

#compare barcodes in both and their respective clusterings
#aka how similar is the clusteirng of these two

from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, confusion_matrix
import seaborn as sns
import matplotlib.pyplot as plt

ari = adjusted_rand_score(merged_df['BEST.GUESS'], merged_df['individualID'])
nmi = normalized_mutual_info_score(merged_df['BEST.GUESS'], merged_df['individualID'])

print(f'Adjusted Rand Index: {ari}')
print(f'Normalized Mutual Information: {nmi}')

differently_classified_count = (merged_df['BEST.GUESS'] != merged_df['individualID']).sum()
print(f"Number of differently classified members: {differently_classified_count}")


# Plot confusion matrix
conf_matrix = confusion_matrix(merged_df['individualID'], merged_df['BEST.GUESS'])
sns.heatmap(conf_matrix, annot=True, fmt='d')
plt.xlabel('FreeMUXLET Cluster')
plt.ylabel('Annotation Cluster')
# plt.show()
plt.savefig("heatmapOfFreemuxOverlap070124.jpg")



freemuxClusters = {}
ROSMAPClusters = {}

for index, row in merged_df.iterrows():
	if row['BEST.GUESS'] in freemuxClusters:
		freemuxClusters[row['BEST.GUESS']].append(row['barcode'])
	else:
		freemuxClusters[row['BEST.GUESS']]=[row['barcode']]
	if row['individualID'] in ROSMAPClusters:
		ROSMAPClusters[row['individualID']].append(row['barcode'])
	else:
		ROSMAPClusters[row['individualID']]=[row['barcode']]
import numpy as np
matrix = np.zeros((len(freemuxClusters),len(ROSMAPClusters)), dtype=int)

# Set the element at index (5, 7) to 8

numFree=0
for keyFreemux in freemuxClusters:
	numROS=0
	for keyROSMAP in ROSMAPClusters:
		overlap = [element for element in freemuxClusters[keyFreemux] if element in ROSMAPClusters[keyROSMAP]]
		matrix[numFree, numROS] = len(overlap)
		numROS=numROS+1
	numFree=numFree+1

print(list(freemuxClusters.keys()))
print(list(ROSMAPClusters.keys()))
