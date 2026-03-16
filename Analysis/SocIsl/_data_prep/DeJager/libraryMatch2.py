import csv
import os
result_table = []

with open("cell-annotation.csv", 'r') as csvfile:
	csvreader = csv.reader(csvfile)
	next(csvreader)
	for row in csvreader:
		patient = row[1]
		key = row[0]
		# print(key)
		keyOverall = key.split('_')
		library = keyOverall[0]
		barcode = keyOverall[1]
		result_table.append({"Cell Barcode": barcode, "Assigned Patient": patient, "Library": library})


with open("ROSMAP_clinical.csv", 'r') as csvfile:
	csvreader = csv.reader(csvfile)
	for row in csvreader:
		value = row[0]
		checkHere = row[len(row)-3]
		for entry in result_table:
			if entry["Assigned Patient"] == checkHere:
				entry["Assigned Patient"] = value

print(result_table[:5])

output_csv="newCellAnno.csv"

file_exists = os.path.exists(output_csv)
with open(output_csv, mode="a", newline="") as file:
	fieldnames = ["Cell Barcode", "Assigned Patient", "Library"]
	writer = csv.DictWriter(file, fieldnames=fieldnames)
	if not file_exists:
		writer.writeheader()
	for row in result_table: 
		writer.writerow(row)
# # cellAnno=pd.read_csv("cell-annotation.csv")
# result_dict = {}

# # Load the CSV file
# fastqs=[]

# with open("Fastq_paths_432_PFC_HM_updated_edited.csv", 'r') as csvfile:
# 	csvreader = csv.reader(csvfile)
# 	for row in csvreader:
# 		if len(row)>1:
# 			fastqs.append(row[1])

# print(fastqs)

# filtered_list=[]

# with open("scPFC_432_withACEandSIvariables.csv", 'r') as csvfile:
# 	csvreader = csv.reader(csvfile)
# 	for row in csvreader:
# 		filtered_list.append(row[5])

# print(len(list(set(filtered_list).intersection(libraries))))
# with open("dataset_652_basic_12-23-2021.csv", 'r') as csvfile:
# 	csvreader = csv.reader(csvfile)
# 	for row in csvreader:
# 		if len(row)>1 and row[0] in fastqs:
# 			if len(row[54]) >= 1:
# 				filtered_list.append(row[0])

# filtered_list.remove("projid")
# folders = [item for item in filtered_list if os.path.isdir(os.path.join("/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/Preprocessed_Counts/Resilient", item))]
# folders2 = [item for item in filtered_list if os.path.isdir(os.path.join("/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/Preprocessed_Counts/ACE", item))]



# total_samples=[]

# with open("WGS_sample_QC_info.csv", 'r') as csvfile:
# 	csvreader = csv.reader(csvfile)
# 	for row in csvreader:
# 		value = row[1]
# 		checkHere = row[0]
# 		total_samples.append(checkHere)
# 		# for i in result_dict:
# 		# 	result_dict[i] = [value if x == checkHere else x for x in result_dict[i]]


# filtered_list=[]

# with open("dataset_652_basic_12-23-2021.csv", 'r') as csvfile:
# 	csvreader = csv.reader(csvfile)
# 	for row in csvreader:
# 		if len(row)>1 and row[0] in total_samples:
# 			if len(row[54]) >= 1:
# 				filtered_list.append(row[0])

# with open("ROSMAP_clinical.csv", 'r') as csvfile:
# 	csvreader = csv.reader(csvfile)
# 	for row in csvreader:
# 		value = row[0]
# 		checkHere = row[len(row)-3]
# 		print(row)
# 		filtered_list = [checkHere if x == value else x for x in filtered_list]
# 		# for i in result_dict:
# 		# 	result_dict[i] = [value if x == checkHere else x for x in result_dict[i]]

# filtered_list.pop(0)

# result_dict={}

# with open("cell-annotation.csv", 'r') as csvfile:
# 	csvreader = csv.reader(csvfile)
# 	for row in csvreader:
# 		value = row[1]
# 		key = row[0]
# 		key = key.split('_')[0]
# 		print(value)
# 		if value in filtered_list:
# 			if key in result_dict:
# 				if value not in result_dict[key]:
# 					result_dict[key].append(value)
# 			else:
# 				result_dict[key]=[]
# 				result_dict[key].append(value)

# # count=0
# # with open("dataset_652_basic_03-23-2022.csv", 'r') as csvfile:
# #	 csvreader = csv.reader(csvfile)
# #	 for row in csvreader:
# #		 if len(row) > 54 and row[54].strip() != "":
# #			 if any(row[0] in v for v in result_dict.values()):
# #				 count=count+1
# #				 print(row[54])

# # print(count)


# result_dict.keys()

# keyList = result_dict.keys()
# remove_items = [
# 	'200916-B54-B', '200313-B22-B', '200225-B10-B', '200316-B24-B', '200810-B45-B',
# 	'200305-B15-B', '190409-B5-B', '200317-B26-B', '200306-B16-B', '191219-B9-B',
# 	'200713-B33-B', '200312-B20-B', '200309-B17-A', '201007-B57-B', '200810-B47-B',
# 	'200317-B27-B', '190403-B4-B', '200707-B30-B', '200226-B11-B', '200730-B41-B',
# 	'201007-B58-B', '200303-B14-B', '200701-B28-B', '200313-B23-B', '200715-B35-B',
# 	'200804-B42-B', '200708-B31-B', '200310-B18-B', '200810-B46-B', '200930-B55-B',
# 	'201022-B61-B', '201024-B59-B', '200728-B39-B', '201002-B56-B', '200702-B29-B',
# 	'200806-B44-B', '200311-B19-B', '191122-B6-R7090624-alone'
# ]

# remove_items add six_or_more

# import os

# filtered_list = [item for item in keyList if item not in remove_items]

# filtered_list = [f for f in filtered_list if not f.startswith("MAP")]
# filtered_folders = [
# 	f for f in filtered_list
# 	if len(os.listdir(os.path.join("/om/scratch/Tue/shared_folder/WGS", f))) < 6
# ]


# def get_prefix(name):
# 	return "-".join(name.split("-")[:-1])

# remove_prefixes = {get_prefix(r) for r in remove_items}

# filtered_folders = [
# 	f for f in filtered_folders
# 	if get_prefix(f) not in remove_prefixes
# ]

# base_dir = "/om/scratch/Tue/shared_folder/WGS"

# less_than_6 = []
# six_or_more = []

# for f in filtered_list:
# 	num_files = len(os.listdir(os.path.join(base_dir, f)))
# 	if num_files < 6:
# 		less_than_6.append(f)
# 	else:
# 		six_or_more.append(f)

# merged = remove_items + six_or_more
# remove_prefixes2 = {get_prefix(r) for r in six_or_more}


# filtered_folders = [
# 	f for f in filtered_folders
# 	if get_prefix(f) not in remove_prefixes2
# ]


# for key in result_dict:
# 	print(key)
# 	with open('/om/scratch/Tue/mabdel03/WGS/individ/individPat'+key+'.txt', 'w') as f:
# 		for item in result_dict[key]:
# 			f.write(item + '\n')


#WGS - find all the samples in the VCF
#convert them to proj IDs, filter for which have soc isl data - filter the sample list as well
#list the individual IDs and the corresponding unique libraries with soc isl data



