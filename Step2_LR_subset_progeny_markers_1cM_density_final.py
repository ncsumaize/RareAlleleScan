# -*- coding: utf-8 -*-
import os as os
import pandas as pd
import numpy as np

#Get the file with positions of markers every 1 cM
os.chdir("Q:/My Drive/Teo and landrace/Landrace/Rare_Allele_Scan") # Jim's path
#os.chdir("G:/My Drive/Teo and landrace/Landrace/Rare_Allele_Scan") # LFS's path
LRmap = pd.read_csv("LR_genetic_map.csv")

#Get the file with haplotype calls from all progrenies, this is a BIG file
#hap_calls = pd.read_table("C:/Users/jholland/Documents/landrace/Imputed Tuxpeno GBS 2-19-17 from PJB/parentcalls_itr2_Combined_LR13_LR14_parents2_AGPv4_filter0219.hmp.txt.gz", compression='gzip', dtype = str, header = 0, na_values = "N")
hap_calls = pd.read_table("Q:/My Drive/Teo and landrace/Inbreeding Depression ms/Final data files/parentcalls_itr2_Combined_LR13_LR14_parents2_AGPv4_filter0219.hmp.txt.gz", compression='gzip', dtype = str, header = 0, na_values = "N") # Jim's path
#hap_calls = pd.read_table("G:/My Drive/Teo and landrace/Inbreeding Depression ms/Final data files/parentcalls_itr2_Combined_LR13_LR14_parents2_AGPv4_filter0219.hmp.txt.gz", compression='gzip', dtype = str, header = 0, na_values = "N") # LFS's path
hap_calls.iloc[:10,:20]

colnames = hap_calls.columns.tolist()
IDs = colnames[11:] #get the individual progeny IDs, first 11 columns are marker metadata
print(IDs[:10])
marker_list = hap_calls['rs#'].tolist()
marker_list[:10]

#subset the markers to include only the markers at every 0.2 cM positions
#note that the linkage map positions are Agpv4
hap_calls_sub = hap_calls.loc[hap_calls["rs#"].isin(LRmap.marker),:]
print(hap_calls_sub.shape)


#transpose the hap_calls data frame so progenies are in rows, markers in columns
hap_callsT = hap_calls_sub.transpose()
hap_callsT.columns = hap_callsT.iloc[0]
hap_callsT = hap_callsT.drop(hap_callsT.index[0:11])
hap_callsT.head()

#Get the list of IDs from hap_callsT in same format as trait data, and assign that as index
hap_calls_IDs = list()
for id in hap_callsT.index:
    hap_calls_IDs.append(id[:10])

hap_callsT.index = hap_calls_IDs
hap_callsT.index.name  = "PlantID"
hap_callsT.head()

#Now subset the individuals to include only those also in the GWAS data set
#we could read in the complete GWAS marker file, which is huge, or the covariate file, which is small
#phenos = pd.read_csv("C:/Users/jholland/Google Drive/Teo and landrace/Landrace/LR-teosinte_datasets_for_GWAS/ASReml_jobs/Landraces_FL1314_all_data_updated_4398taxa_Fest_2.csv", index_col = 0, na_values = '.')
phenos = pd.read_csv("Q:/My Drive/Teo and landrace/Inbreeding Depression ms/Final data files/Landraces_FL1314_all_data_updated_4398taxa_Fest_2_ed.csv", index_col = 0, na_values = '.') # Jim's path
#phenos = pd.read_csv("G:/My Drive/Teo and landrace/Inbreeding Depression ms/Final data files/Landraces_FL1314_all_data_updated_4398taxa_Fest_2_ed.csv", index_col = 0, na_values = '.') # LFS's path
print(phenos.head())
print(phenos.shape)

IDsubset = list(phenos.index) #will get the taxa IDs that are actually in the phenotype

#Subset the crossover data set to include only the individuals in IDsubset
hap_calls_final = hap_callsT.loc[IDsubset,:]
print(hap_calls_final.iloc[:5,:5])
print(hap_calls_final.shape)

#hap_calls_final.to_csv("C:/Users/jholland/Google Drive/Teo and landrace/Landrace/Rare_Allele_Scan/Hap_calls_LR.csv")
hap_calls_final.to_csv("final_outputs/Hap_calls_LR.csv")
##############
