# -*- coding: utf-8 -*-
'''
A script to take in parentage information and haplotype calls from Peter Bradbury's imputed progeny haplotypes and:

(I already reduced the input marker file to 1 marker per ~ 0.2 cM)
 
1. Recode haplotype calls for each parent in turn to indicate number of copies of a particular haplotype from that parent
2. Output the recoded parent haplotype calls to a file

A subsequent script will be used to implement the scan on marker alleles 
'''
import os as os
import pandas as pd
import numpy as np
import time as time

#Get the reduced haplotype call file
#os.chdir("G:/My Drive/Teo and landrace/landrace/Rare_Allele_Scan")  # LFS's path
os.chdir("Q:/My Drive/Teo and landrace/landrace/Rare_Allele_Scan")  # Jim's path
hap_callsT = pd.read_csv("final_outputs/Hap_calls_LR.csv", index_col = 0)
print(hap_callsT.shape)

#Are any progeny IDs repeated?
IDs = list(hap_callsT.index)
print(IDs[:5])
print(len(IDs))

#make a set of progeny IDs and compare its length
IDset = set(IDs)
print(len(IDset))

#OK, all is good, none of the progenies are duplicated,  data are ready to go!

#Get trait data
#traitdf = pd.read_csv("C:/Users/jholland/Google Drive/Teo and landrace/Landrace/LR-teosinte_datasets_for_GWAS/ASReml_jobs/Landraces_FL1314_all_data_updated_4398taxa_Fest_2.csv", index_col = 0, na_values = '.')
#traitdf = pd.read_csv("G:/My Drive/Teo and landrace/Inbreeding Depression ms/Final data files/Landraces_FL1314_all_data_updated_4398taxa_Fest_2_ed.csv", index_col = 0, na_values = '.') # LFS's path
traitdf = pd.read_csv("Q:/My Drive/Teo and landrace/Inbreeding Depression ms/Final data files/Landraces_FL1314_all_data_updated_4398taxa_Fest_2_ed.csv", index_col = 0, na_values = '.') # Jim's path
hap_trait = pd.merge(traitdf, hap_callsT, how='left', left_index=True, right_index=True)

print(traitdf.shape)
print(hap_callsT.shape)
print(hap_trait.shape)

#4398 individuals were used in final data set for GWAS, and included in the trait data here
#After merging, should be down to 4398 individuals

parent_list = list(set(hap_trait["Mother"].append(hap_trait["Father"])))
all_markers = list(hap_callsT.columns)

#Make dummy variables for parents, this will be used later to compute haplotype calls specific to each parent
def parent_index(x, par): #make a function to return 1 when parent info matches current parent
    if x == par: return 1
    else: return 0

for parent in parent_list:
    mom_index = "mom_" + parent #create the new column name
    dad_index = "dad_" + parent
    hap_trait[mom_index] = hap_trait["Mother"].apply(parent_index, args=(parent,))
    hap_trait[dad_index] = hap_trait["Father"].apply(parent_index, args=(parent,))

#first check if all data are missing at any markers, remove them from data set if so
for marker in all_markers:
    if hap_trait[marker].isnull().all():
        hap_trait = hap_trait.drop(marker, index = 1)
        all_markers.remove(marker)

#This may not be necessary, but we will break the marker calls into separate files for each chromosome
#This will reduce amount of RAM allocated to holding data while processing, but at a cost of more file I/O

markers_by_chrom = dict()
for marker in all_markers:
    chrom = marker[1:3]
    if chrom[1] == '_': chrom = chrom[0]
    if chrom not in markers_by_chrom:
        markers_by_chrom[chrom] = [marker]
    else: markers_by_chrom[chrom].append(marker)

# Now recode each marker for individual parent allele effects.
# Here is how the data are coded in the hap_calls data frame to indicate parental haplotype inheritance, this is from an email from Peter Bradbury 1/5/2017:
# "A = one allele from parent 1, chromosome 1 and one from parent 2 chromosome 1 C = parent 1 chromosome 1, parent 2 chromosome 2 C = parent 1 chromosome 2, parent 2 chromosome 1 C = parent 1 chromosome 2, parent 2 chromosome 2
# Where flanking sites were identical the missing data between them was filled in with the flanking values, which is why there is not much missing. This was an intermediate data set that was generated. The imputation was done by replacing those with parent haplotypes. The non-missing data in the imputed set are the sites that were phased in both parents. The main value of the parent calls after the imputation is that a change in state, A to C, etc, indicates a crossover."
# "To be sure, the last two cases are "G" and "T", right (you have "C" and "C" above) - jim"
# "Yes, A,C,G,T not ACCC. The only thing missing for linkage is the identity of parent1 and parent2. That would come from the parentage file used for the analysis.
# Peter"
# "So, to be sure I understand:
# parent calls 'A' and 'C' when parent 1 is used as female involves the same chromosome segment (haplotype) as parent calls 'A' and 'G' when parent 1 is used as male in another cross?
# or is the coding arbitrary within families, such that the chromosome segments from parent 1 involved in 'A' and 'C' calls in one family may involve the same chromosome as 'C' and 'T' in another family from the same parent (in which case we would need to do some work to sort out identity of haplotypes across families involving the same parent)?
# thanks - jim"
# "Jim,
# The first case. The chromosome order within parent is determined by the parent phasing data, which is extracted in a previous step. The assignment of chromosomes is arbitrary at that point but that assignment gets applied to all the parent calls in the imputation step.
# Peter"
# 

MinCount = 20 #Set the minimum count for individuals with at least one copy of allele being tested
            #I am not sure what is most appropriate, but I think we will have trouble for cases with 
            #very low sample sizes in the allele class

chrom_list = [str(i) for i in range(1,11)]
#chrom_list = ['1'] #for testing

for chrom in chrom_list:
    print('Chromosome ' + chrom)
    start_time = time.perf_counter()
    print('Starting time: ' + str(start_time))
    
    res = [] #list to hold results for current chrom

    marker_list = markers_by_chrom[chrom]

    #testing code   
    #parent_list = ["164_2"]
    #marker_list = ["S1_222588"]
    #parent_list = parent_list[:10]
    #marker_list = marker_list[:20]

    for marker in marker_list:
    
        #print(marker)
    
        for parent in parent_list:
        
            #recode the marker calls
        
            #for this parent and marker combination, we have to test two haplotypes
            #geno codes are A = mother allele 1/dad allele 1, C = mother allele 1/dad allele 2, 
            #G = mother allele 2/dad allele 2, T = mother allele 2/dad allele 2
            #for each individual, we make two new columns
            #hap1 = 0, 1, 2 counts of haplotypes within each individual for haplotype 1 from parent
            #hap2 = 0, 1, 2 counts of haplotypes within each individual for haplotype 2 from parent
            
            mom_index = "mom_" + parent            
            dad_index = "dad_" + parent
            
            get_columns = ["Mother", "Father", marker, mom_index, dad_index]
            tempdf = hap_trait[get_columns]
            
            #Map the letter calls to four numeric columns representing the two haplotypes of each parent
            #create a dictionary that maps alphabetic calls to the maternal haplotype 1 (either 0 or 1)
            mom_dict = {"A":1, "C": 1, "G":0, "T":0}
            #do the same for paternal haplotype 1
            dad_dict = {"A":1, "C": 0, "G":1, "T":0}            
            #note that the numeric indicator of the 2nd haplotype is just 1 - dad1 and 1 - dad1 (so we don't need to explicitly create them)
            tempdf["mom1"] = tempdf[marker].map(mom_dict)
            tempdf["dad1"] = tempdf[marker].map(dad_dict)
            
            #Haplotype calls are the sum of parental indices times numeric haplotype calls
            hap1 = marker + "_" + parent + "hap1" #create the column name for haplotype calls in the final data set
            hap2 = marker + "_" + parent + "hap2"
            tempdf[hap1] = (tempdf.mom1*tempdf[mom_index]) + (tempdf.dad1*tempdf[dad_index])
            tempdf[hap2] = ((1 - tempdf.mom1)*tempdf[mom_index]) + ((1 - tempdf.dad1)*tempdf[dad_index])          
    
            if tempdf[hap1].dropna().astype(bool).sum() < MinCount: 
                tempdf = tempdf.drop(hap1, axis = 1)          
            if tempdf[hap2].dropna().astype(bool).sum() < MinCount: 
                tempdf = tempdf.drop(hap2, axis = 1)          
                              
            #drop unnecessary columns
            tempdf = tempdf.drop(get_columns + ["mom1", "dad1"], axis = 1)  
            #append the new columns as dataframe to res list
            if tempdf.shape[1] >= 1:
                res.append(tempdf)
            #tempdf.to_csv("tempdf.csv")  #for checking results 
                
    resdf = pd.concat(res, axis = 1) #, columns = ["Mother", "Father", "marker", "hap1_A","hap1_D", "hap2_A", "hap2_D"])
    
    filename = "final_outputs/LR_rare_allele_hap_calls_chr" + chrom + ".csv"
    resdf.to_csv(filename, index = True)

    end_time = time.perf_counter()
    elapsed = round(end_time - start_time)
    print('Elapsed Time for Chromosome ' + chrom + ": "  + str(elapsed))
    



