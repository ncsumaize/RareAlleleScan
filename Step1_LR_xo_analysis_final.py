# -*- coding: utf-8 -*-
import os as os
import pandas as pd
import numpy as np
#Get the file with crossover positions in each progeny of Tuxpeno population
#os.chdir("C:/Users/jholland/Documents/landrace/Imputed Tuxpeno GBS 2-19-17 from PJB")
os.chdir("G:/My Drive/Teo and landrace/landrace/Rare_Allele_Scan")  # LFS's path
#os.chdir("Q:/My Drive/Teo and landrace/landrace/Rare_Allele_Scan")  # Jim's path

data_types =  {'taxon': str, 'chr': int, 'start': int, 'end': int}
#xo = pd.read_table("Q:/My Drive/Teo and landrace/Inbreeding Depression ms/Final data files/xo_Combined_LR13_LR14_parents2_AGPv4_filter0219.txt", delim_whitespace = True, dtype=data_types) # Jim's path
xo = pd.read_table("G:/My Drive/Teo and landrace/Inbreeding Depression ms/Final data files/xo_Combined_LR13_LR14_parents2_AGPv4_filter0219.txt", delim_whitespace = True, dtype=data_types) # LFS's path
print(xo.head)
print(xo.columns)

#slice out the progeny ID
xo = xo.assign(ID = xo.taxon.str.slice(start = 0, stop = 10))
print(xo.head)

#some of the progenies are repeated, but a quick check shows that Peter Bradbury already merged duplicates
taxaset = set(xo.taxon)
print(len(taxaset))

IDset = set(xo.ID)
print(len(IDset)) #same as taxaset, so we are cool, we will not double count XOs in some individuals


#Also read in the list of progenies actually used in GWAS, there are 4398 of them, these exclude individuals with high crossover rates and other problems
#we want to only keep the same set of individuals in this set for counting crossovers.

#we could read in the complete GWAS marker file, which is huge, or the covariate file, which is small
#GWAScovs = pd.read_table("C:/Users/jholland/Google Drive/Teo and landrace/Landrace/GWAS_input_files_tassel5/LR1314_4398_covs.txt", skiprows=2, header=None)
#phenos = pd.read_csv("Q:/My Drive/Teo and landrace/Inbreeding Depression ms/Final data files/Landraces_FL1314_all_data_updated_4398taxa_Fest_2_ed.csv", index_col = 0, na_values = '.') # Jim's path
phenos = pd.read_csv("G:/My Drive/Teo and landrace/Inbreeding Depression ms/Final data files/Landraces_FL1314_all_data_updated_4398taxa_Fest_2_ed.csv", index_col = 0, na_values = '.') # LFS's path


IDsubset = phenos.index #phenos.iloc#[:,0].str.slice(start = 0, stop = 10)

#Subset the crossover data set to include only the individuals in IDsubset
xo = xo.loc[xo['ID'].isin(IDsubset),:]

#Also sort the xo data frame by the crossover start position
xo = xo.sort_values(by = 'start')


#Now for each chromosome, scan across sequence, cumulating crossovers as we go
#for N progenies * 2 chromosomes each , 
#we need NXO = 0.01*2N xos to equal 1 cM 
#We will use the start positions of the xos to equal the cM position
#Since the list is ordered by xo position, this is easy, just take every NXO row!
NXO = int(round(2*0.01*len(IDsubset),0))
result = {}

for chrom in range(1,11):
    print("chrom: " + str(chrom))
    xox = xo.copy()
    xox = xox.loc[xox['chr'] == chrom,]
    xox.reset_index(inplace=True) #reset the index of the subset df
    nrow = xox.shape[0]
    print(nrow)
    
    #generate a sequence of every tenth index value 
    selex = range(0,nrow,NXO)
    #print(list(selex)[-1])
    xoxmap = pd.DataFrame(xox.loc[selex,"start"])
    #print(xoxmap.head)
    #print(xoxmap.tail)
    map_length = (nrow)#one marker per cM
    print(map_length)

    #it's possible that the same position appears twice
    xoxmap = xoxmap.assign(chrom = chrom)
    xoxmap = xoxmap.assign(cM = list(np.arange(0, np.ceil(map_length)))[:len(selex)])
    
    #put the map df into result list
    result[chrom] = xoxmap
    

#Now push contents of result into a data frame
finaldf = result[1]
for chrom in range(2,11): finaldf = finaldf.append(result[chrom])
print(finaldf.head)

#print out the dimensions of the map before filtering, this will be equal to total map length
print("Total Map Length:")
print(finaldf.shape)

#drop duplicated markers (one marker representing > 1 cM interval) but keeping proper linkage map positions
finaldf = finaldf.drop_duplicates(subset=['chrom', 'start'], keep = 'last')

#Now merge the linkage map to the marker info
#marker_info = pd.read_table("Q:/My Drive/Teo and landrace/Inbreeding Depression ms/Final data files/parentcalls_itr2_Combined_LR13_LR14_parents2_AGPv4_filter0219.hmp.txt.gz",  delim_whitespace = True, usecols = [0,2,3]) # # Jim's path
marker_info = pd.read_table("G:/My Drive/Teo and landrace/Inbreeding Depression ms/Final data files/parentcalls_itr2_Combined_LR13_LR14_parents2_AGPv4_filter0219.hmp.txt.gz",  delim_whitespace = True, usecols = [0,2,3]) # # LFS's path
marker_info = marker_info.rename(columns = {'rs#':'marker', 'pos':'AGPv4_Pos'})

finaldf_map = marker_info.merge(finaldf, how = "right", left_on = ['chrom', 'AGPv4_Pos'], right_on = ['chrom', 'start'])

#keep only the columns of interest
finaldf_map = finaldf_map.loc[:,['marker', 'chrom', 'AGPv4_Pos', 'cM']] 
 
#finaldf_map.to_csv("Q:/My Drive/Teo and landrace/Inbreeding Depression ms/Final data files/LR_genetic_map.csv", index = False)
finaldf_map.to_csv("final_outputs/LR_genetic_map.csv", index = False)
