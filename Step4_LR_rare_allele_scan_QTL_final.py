# -*- coding: utf-8 -*-
'''
A script to take in recoded haplotype calls for each parent (indicating number of copies of a particular haplotype from that parent)

Generate dummmy variables to serve as covariates for parentage effects

Check for each marker if it meets minimum count of rare allele: if no, return nothing
Check for each marker if there are homozygous rare types
    if so test a + d effects together
    if not, test a effect only
    
Return additive and dominance effect estimates and p-values into a data frame
'''
import os as os
import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
import time as time

#set working directory
#os.chdir("G:/My Drive/Teo and landrace/Landrace/Rare_Allele_Scan") # LFS's path
os.chdir("Q:/My Drive/Teo and landrace/Landrace/Rare_Allele_Scan") # Jim's path

#set global values
trait_list= ['DTA','DTS','PltHt', 'LfLen', 'LfWth', 'TillN', 'HkLfN', 'pEarN', 'ShkLen', 'ShkIntLen', 'CobLen', 'K1rowN', 'CobDia', 'pKN', 'EarInternodeLen', 'tKN', 'tKW', 'GW']
chrom_list = list(range(1,11))

#make a dictionary that will translate trait names from the data file to the acronyms used in the paper
trait_dict = {'DTA':'DTA', 'DTS':'DTS','PltHt':'PLHT', 'LfLen':'LFLN', 'LfWth':'LFWD', 
              'TillN':'TILN', 'pEarN':'PROL', 'HkLfN':'LBNN', 'ShkLen':'LBLN', 
              'ShkIntLen':'LBIL', 'CobLen':'EL', 'K1rowN':'CUPR', 'CobDia':'ED', 'pKN':'GE',
              'EarInternodeLen':'EILN', 'tKN':'TGPP', 'tKW':'TGWP', 'GW':'GW'}

#get trait data and covariates
traitdf = pd.read_csv("Q:/My Drive/Teo and landrace/Inbreeding Depression ms/Final data files/Landraces_FL1314_all_data_updated_4398taxa_Fest_2_ed.csv", index_col = 0, na_values = '.') # Jim's path
#traitdf = pd.read_csv("G:/My Drive/Teo and landrace/Inbreeding Depression ms/Final data files/Landraces_FL1314_all_data_updated_4398taxa_Fest_2_ed.csv", index_col = 0, na_values = '.') # Jim's path
traitdf = traitdf.rename(columns = {'F.est':'F_est'})
print(traitdf.shape)

parent_list = list(set(traitdf["Mother"].append(traitdf["Father"])))

#==================================================================================================================>>>>>>
#Create dummy variables for parent effects on each progeny
#GCA_columns = []

#Create dummy variable columns for parent effects
#for parent in parent_list:
#    col_name = "GCA_" + parent 
#    GCA_columns.append(col_name)        
#    traitdf[col_name] = 0  
#    traitdf.loc[(traitdf.Mother == parent), col_name] = traitdf.loc[(traitdf.Mother == parent), col_name] + 1
#    traitdf.loc[(traitdf.Father == parent), col_name] = traitdf.loc[(traitdf.Father == parent), col_name] + 1

#make a string of all the elements of GCA_columns
#GCA_columns_str = " + ".join(GCA_columns)
#print(GCA_columns_str)
#===================================================================================================================>>>>>

#here's an example of how to make a string that represents 50 PC variables that are named in the data frame as PC1, PC2, ... PC50
PCnums = list(range(1,51))
PCstring = " + PC".join(map(str, PCnums))
PCstring = "PC" + PCstring
#print(PCstring)

#Define a function to get dominance coefficients for a marker
def dominance(marker):
    return marker%2

#Define function to test one parental haplotype against all others, adjusting for parental effects and non-genetic effects
def rareQTL(col):
 
    '''
    This is the function that will be applied to each column (col) of the haplotype calls data frame.
    It takes one column of haplotype additive copy calls for one parent and one haplotype one marker
    Creates a second column variable for the dominance effects
    Then fits the model with all the covariates plus the A and D effects for the marker
    NOTE THAT trait IS SET AS A GLOBAL VARIABLE OUTSIDE OF THIS FUNCTION
    '''   
         
    #check if there are homozygotes for rare allele, if not, we cannot estimate dominance
    if col[col == 2].count() == 0:
        #print('no rare homozygotes')
        form = trait + ' ~ C(Year) + Tire_track:C(Year) + Shading:C(Year) + F_est:C(Year) + x1_1 + x1_2 + x1_3 + x1_4 + x2_1 + x2_2 + x2_3 + x2_4 +	y1_1 +	y1_2 +	y1_3 +	y1_4 +	y2_1 +	y2_2 +	y2_3 +	y2_4 + ' + PCstring + ' + ' + col.name 
        #print(form)
        results = smf.ols(formula = form, data=hap_trait).fit()
        RSSmark = results.ssr
        #LOD = (n/2)log10(RSS0/RSSmark)
        LOD = (n/2)*np.log10(RSS0/RSSmark)
        #print(results.summary())
        return pd.Series({'marker' : col.name, 'a': results.params[col.name], 'p_value_add': results.pvalues[col.name], 'd':np.nan, 'p_value_d': np.nan, 'LOD': LOD})
      
       
    #otherwise, include the dominance effect
    form = trait + ' ~ C(Year) + Tire_track:C(Year) + Shading:C(Year) + F_est:C(Year) + x1_1 + x1_2 + x1_3 + x1_4 + x2_1 + x2_2 + x2_3 + x2_4 +	y1_1 +	y1_2 +	y1_3 +	y1_4 +	y2_1 +	y2_2 +	y2_3 +	y2_4 +' + PCstring + ' + ' + col.name + ' + ' + 'dominance(' + col.name + ')'
    #print(form)
    #print(form)
    results = smf.ols(formula = form, data=hap_trait).fit()
    RSSmark = results.ssr
    LOD = (n/2)*np.log10(RSS0/RSSmark)
    
    #print('RSSmark: ' + str(RSSmark))
    #print('No. obs: ' + str(results.nobs))
    #print('LOD score: ' + str(LOD))

    #print(results.summary())
    #print(results.params)
    return pd.Series({'marker' : col.name, 'a': results.params[col.name], 'p_value_add': results.pvalues[col.name], 'd':results.params['dominance(' + col.name + ')'], 'p_value_d': results.pvalues['dominance(' + col.name + ')'], 'LOD': LOD})

#############################################################
#OUTERMOST LOOP STARTS HERE
#EXECUTE A COMPLETE ANALYSIS FOR ONE TRAIT AT A TIME
############################################################

for trait in trait_list:
    print("Starting analysis of trait " + trait)
    
#############################################################
#FIT BASE MODEL WITH COVARIATES, NO MARKERS
#Get number of obs and SSresidual to use in LOD computations
############################################################
    base_form = trait + ' ~ C(Year) + Tire_track:C(Year) + Shading:C(Year) + F_est:C(Year) + x1_1 + x1_2 + x1_3 + x1_4 + x2_1 + x2_2 + x2_3 + x2_4 +	y1_1 +	y1_2 +	y1_3 +	y1_4 +	y2_1 +	y2_2 +	y2_3 +	y2_4 +' + PCstring 
    base_results = smf.ols(formula = base_form, data=traitdf).fit()
    n = base_results.nobs
    RSS0 = base_results.ssr

    print('BASE MODEL no. obs: ' + str(n))
    print('BASE MODEL RSS0:' + str(RSS0))

#############################################################
          
            
 #############################################################
 #SECOND LOOP STARTS HERE
 #EXECUTE A COMPLETE ANALYSIS FOR ONE CHROMOSOME AT A TIME
 ############################################################
                                                  
    #Execute the functions and save results for each chromosome separately
    for chrom in chrom_list:

    #Get the haplotype calls file for current chromosome
    #I read in the row names in a separate file, the reason for this is that rownames are str, whereas all marker data are float
    #(another issue, can't use int because there are NaN values)
    #setting dtype in read_csv should make it faster
        print('Chromosome ' + str(chrom))
        start_time = time.perf_counter()

    #to get the whole file, use this:
    #for pandas version 20:
    #hap_calls = pd.read_csv("Teo_rare_allele_hap_calls_chr" + str(chrom) + ".csv", index_col = None, usecols=lambda x: x not in ['PlantID'] , dtype = float)

    #for older versions:
    #first read the header line to count columns
                                        
#        with open("J:/LFS_2/LANDRACES_2/Rare_alleles_scan/LR/Rare_Allele_Scan/LR_rare_allele_hap_calls_chr" + str(chrom) + ".csv", 'r') as f:
        with open("final_outputs/LR_rare_allele_hap_calls_chr" + str(chrom) + ".csv", 'r') as f:
            first_line = f.readline()

        ncols = first_line.count(',') + 1  
    
        #hap_calls = pd.read_csv("J:/LFS_2/LANDRACES_2/Rare_alleles_scan/LR/Rare_Allele_Scan/LR_rare_allele_hap_calls_chr" + str(chrom) + ".csv", index_col = None, usecols = list(range(1,ncols)), dtype = float)
        hap_calls = pd.read_csv("final_outputs/LR_rare_allele_hap_calls_chr" + str(chrom) + ".csv", index_col = None, usecols = list(range(1,ncols)), dtype = float)


        #fill in missing calls with zeroes!
        hap_calls.fillna(0, inplace = True)

        after_hap = time.perf_counter()
        elapsed = round(after_hap - start_time)
        print('Time to read hap calls for Chromosome ' + str(chrom) + ": "  + str(elapsed))

        #now separately pull in the index
        #hap_calls.index = pd.read_csv("J:/LFS_2/LANDRACES_2/Rare_alleles_scan/LR/Rare_Allele_Scan/LR_rare_allele_hap_calls_chr" + str(chrom) + ".csv", index_col = 0, usecols = [0,1]).index #read in one extra column then slice it out!
        hap_calls.index = pd.read_csv("final_outputs/LR_rare_allele_hap_calls_chr" + str(chrom) + ".csv", index_col = 0, usecols = [0,1]).index #read in one extra column then slice it out!

        marker_list = list(hap_calls.columns)
  
        #would be smart to time this vs. single read in approach on a large file
        after_index = time.perf_counter()
        elapsed = round(after_index - after_hap)
        print('Time to read index for Chromosome ' + str(chrom) + ": "  + str(elapsed))
 
 
        hap_trait = pd.merge(traitdf, hap_calls, how='left', left_index=True, right_index=True)

        before_apply = time.perf_counter()

        QTL_tests = hap_trait.loc[:,marker_list].apply(rareQTL) #here is where we execute the analysis function

        after_apply = time.perf_counter()
        elapsed = round(after_apply - before_apply)
        print('Elapsed Time for Scan on Chromosome ' + str(chrom) + ": "  + str(elapsed))

        resultsT = QTL_tests.transpose() #transpose results to make it long rather than wide
        trait_code = trait_dict[trait]
        resultsT['Trait'] = trait_code

        resultsT.to_csv("final_outputs/LR_Rare_allele_scan_output_" + trait_code + "_chr_" + str(chrom) + ".csv")
        print('time to write results for chrom ' + str(chrom) + ':' + str(round(time.perf_counter() - after_apply)))

 #############################################################
 #END SECOND LOOP FOR CHROMOSOME 
 ############################################################
  