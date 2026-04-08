# DEPRECATED: Legacy data preparation script from Openmind era.
# The current Processing pipeline (Processing/*/Pipeline/) supersedes this.
# Paths may still reference decommissioned Openmind directories.
#

import anndata as ad
import scanpy as sc
import os
os.environ["HDF5_USE_FILE_LOCKING"] = "TRUE"
import numpy as np

import tempfile
import zarr

from functools import reduce
import scanpy.external as sce

import seaborn as sns
import torch
from rich import print
from pathlib import Path
import harmonypy as hm
import matplotlib.pyplot as plt

import pandas as pd


libraries = ['482428', '2518573', '2899847', '3283241', '3713990', '3889845', '6107196', '6804844', '7253015', '7265221', '7311370', '8109170', '8132197', '9650662', '9841821', '14184286', '14498577', '15218541', '16068769', '17929065', '18414513', '18659212', '18920002', '20907534', '21362537', '22396591', '22776575', '22789958', '22868024', '23004922', '23690880', '24039289', '24680888', '24747976', '25300551', '26569730', '26631069', '27586957', '29933130', '31509843', '31726180', '32383679', '32705437', '33411712', '33501827', '34542628', '34726040', '34962204', '35941263', '36492755', '36830117', '37030589', '37065652', '37178462', '37436329', '37527863', '38264019', '39484737', '39989287', '41285665', '41773404', '42988567', '43074402', '43485807', '44019405', '44299049', '44671043', '45115248', '45566083', '46000440', '46246604', '46251007', '48023497', '49676048', '50100806', '50100932', '50101523', '50102114', '50102376', '50103129', '50105301', '50105699', '50106280', '50106442', '50106730', '50107583', '50107619', '50107745', '50107907', '50108886', '50108912', '50109477', '50109639', '50109927', '50111971', '50197261', '50300084', '50300408', '50300822', '50301099', '50301125', '50301387', '50301675', '50302266', '50302554', '50400259', '50400385', '50400835', '50401002', '50402567', '50402693', '50402729', '50403446', '50404299', '50405330', '50406057', '50409406', '50410319', '50500136', '50500550', '51520126', '51624179', '51815338', '51826377', '52311825', '52940446', '53772202', '54324848', '57978756', '60725338', '60747316', '60848460', '61029627', '61142759', '61344957', '61827429', '62404688', '62574447', '63188799', '64291829', '64336939', '64493027', '65001607', '65002324', '65292866', '65499271', '65652206', '65736039', '66406040', '66754397', '66924745', '67185070', '67429065', '68015667', '68525196', '68539908', '68745332', '69982533', '70153803', '70625336', '70816595', '70883578', '70917649', '71626373', '71648351', '72650337', '72777797', '74284255', '74753465', '75169847', '75990666', '76733461', '76867532', '77143621', '77239958', '77891596', '78353027', '78452313', '80790863', '81086436', '82317494', '82684406', '83034844', '83173570', '83984043', '84417209', '84642424', '84653463', '85065193', '85171938', '86934089', '87410220', '87645445', '87779516', '89001223', '89164957', '89546375', '89903942', '90267190', '90780976', '90942860', '91018909', '91444029', '91707643', '92393245', '92629514', '93815598', '94092977', '94144536', '94430339', '94974890', '95491648', '96095092', '98953007', '86408244', '50102088']
batch_deidentified = pd.read_csv("batch_mapping_deidentified.tsv", sep = '\t')
batchiden = pd.read_csv('full_batch_mapping.tsv', sep='\t')

# ast = ad.read_h5ad("ast.h5ad")
# oli = ad.read_h5ad("oli.h5ad")
mic = ad.read_h5ad("mic.h5ad")

# exc1 = ad.read_h5ad("exc1.h5ad")
# exc2 = ad.read_h5ad("exc2.h5ad")
# exc3 = ad.read_h5ad("exc3.h5ad")
# exc1.layers["raw"] = exc1.raw.X
# exc2.layers["raw"] = exc2.raw.X
# exc3.layers["raw"] = exc3.raw.X
# exc = ad.concat([exc1,exc2,exc3], merge="same")

# exc.raw = ad.AnnData(X=exc.layers["raw"],var=exc.var,obs=exc.obs)
# subjectsAst = ast.obs['subject'].unique()
# subjectsOli = oli.obs['subject'].unique()
subjectsMic = mic.obs['subject'].unique()
# subjectsExc = exc.obs['subject'].unique()

variables=['mic']
#'ast','oli',

# if True:
# 	subjectMatches = {}
# 	for i in subjectsAst:
# 		result = batch_deidentified[batch_deidentified["subject"] == i]
# 		subjectMatches[i]=[]
# 		print(type(result))
# 		for idx, row in result.iterrows():
# 			subjectMatches[i].append(row["dataset"])

# 	subjectMatchesSecond = {}
# 	meta = pd.read_csv('scPFC_432_withACEandSIvariables.csv')
# 	for k in subjectMatches:
# 		subjectMatchesSecond[k]=[]
# 		for l in subjectMatches[k]:
# 			result = batchiden[batchiden["dataset"] == l]
# 			for idx, row in result.iterrows():
# 				subjectMatchesSecond[k].append(row["projid"])

# 	subjectMatchesFinal={}
# 	for k in subjectMatchesSecond:
# 		subjectMatchesFinal[k]=subjectMatchesSecond[k][0]

# 	ast.obs["patient_id"] = ast.obs["subject"].map(subjectMatchesFinal)
# 	ast2 = ast.raw.to_adata()
# 	ast2.var_names = ast2.var['_index']
# 	ast.var_names= ast2.var_names
# 	ast2.var = ast2.var.drop(columns="_index")
# 	ast.raw=ast2
# 	ast.write("ast2.h5ad")

# if True:
# 	subjectMatches = {}
# 	for i in subjectsOli:
# 		result = batch_deidentified[batch_deidentified["subject"] == i]
# 		subjectMatches[i]=[]
# 		print(type(result))
# 		for idx, row in result.iterrows():
# 			subjectMatches[i].append(row["dataset"])

# 	subjectMatchesSecond = {}
# 	meta = pd.read_csv('scPFC_432_withACEandSIvariables.csv')
# 	for k in subjectMatches:
# 		subjectMatchesSecond[k]=[]
# 		for l in subjectMatches[k]:
# 			result = batchiden[batchiden["dataset"] == l]
# 			for idx, row in result.iterrows():
# 				subjectMatchesSecond[k].append(row["projid"])

# 	subjectMatchesFinal={}
# 	for k in subjectMatchesSecond:
# 		subjectMatchesFinal[k]=subjectMatchesSecond[k][0]

# 	oli.obs["patient_id"] = oli.obs["subject"].map(subjectMatchesFinal)
# 	oli2 = oli.raw.to_adata()
# 	oli2.var_names = oli2.var['_index']
# 	oli.var_names= oli2.var_names
# 	oli2.var = oli2.var.drop(columns="_index")
# 	oli.raw=oli2
# 	oli.write("oli2.h5ad")

# if True:
# 	subjectMatches = {}
# 	for i in subjectsMic:
# 		result = batch_deidentified[batch_deidentified["subject"] == i]
# 		subjectMatches[i]=[]
# 		print(type(result))
# 		for idx, row in result.iterrows():
# 			subjectMatches[i].append(row["dataset"])

# 	subjectMatchesSecond = {}
# 	meta = pd.read_csv('scPFC_432_withACEandSIvariables.csv')
# 	for k in subjectMatches:
# 		subjectMatchesSecond[k]=[]
# 		for l in subjectMatches[k]:
# 			result = batchiden[batchiden["dataset"] == l]
# 			for idx, row in result.iterrows():
# 				subjectMatchesSecond[k].append(row["projid"])

# 	subjectMatchesFinal={}
# 	for k in subjectMatchesSecond:
# 		subjectMatchesFinal[k]=subjectMatchesSecond[k][0]

# 	mic.obs["patient_id"] = exc.obs["subject"].map(subjectMatchesFinal)
# 	exc.write("exc2.h5ad")

if True:
	subjectMatches = {}
	for i in subjectsMic:
		result = batch_deidentified[batch_deidentified["subject"] == i]
		subjectMatches[i]=[]
		print(type(result))
		for idx, row in result.iterrows():
			subjectMatches[i].append(row["dataset"])

	subjectMatchesSecond = {}
	meta = pd.read_csv('scPFC_432_withACEandSIvariables.csv')
	for k in subjectMatches:
		subjectMatchesSecond[k]=[]
		for l in subjectMatches[k]:
			result = batchiden[batchiden["dataset"] == l]
			for idx, row in result.iterrows():
				subjectMatchesSecond[k].append(row["projid"])

	subjectMatchesFinal={}
	for k in subjectMatchesSecond:
		subjectMatchesFinal[k]=subjectMatchesSecond[k][0]

	mic.obs["patient_id"] = mic.obs["subject"].map(subjectMatchesFinal)
	mic2 = mic.raw.to_adata()
	mic2.var_names = mic2.var['_index']
	mic.var_names= mic2.var_names
	mic2.var = mic2.var.drop(columns="_index")
	mic.raw=mic2
	mic.write("mic2.h5ad")
