# -*- coding: utf-8 -*-
'''
A script to read in QTL peak positions for one trait

Fit all QTL simultaneously in a linear model with all covariates to a subset of the data, holding out a validation set of individuals' phenotypes

Rank markers based on combined additive + dominance SS

If lowest ranking marker has combined F test p-value > 0.001, drop that marker and re-run model

Iterate until lowest ranking marker has combined A + D p-value < 0.001

Predict phenotypic values from final selected QTL model

Also predict phenotypic values from model with covariates but no QTL

Compare the predicted values to observed values in validation set to see how much the QTL improve prediction ability

'''
import os as os
import pandas as pd
import numpy as np
import re as re
import statsmodels.formula.api as smf
from sklearn.model_selection import KFold
from patsy import dmatrices

#set working directory
#os.chdir("G:/My Drive/Teo and landrace/Landrace/Rare Allele Scan/final_outputs") # LFS's path
os.chdir("Q:/My Drive/Teo and landrace/Landrace/Rare_Allele_Scan/final_outputs") # Jim's path

#set global values
trait_list= ['DTA', 'DTS','PLHT', 'LFLN', 'LFWD', 'TILN', 'PROL', 'LBNN', 'LBLN', 
              'LBIL', 'EL', 'CUPR', 'ED', 'GE', 'EILN', 'TGPP', 'TGWP', 'GW']

#make a dictionary that will translate trait names from the data file to the acronyms used in the paper
trait_dict = {'DTA':'DTA', 'DTS':'DTS','PltHt':'PLHT', 'LfLen':'LFLN', 'LfWth':'LFWD', 
              'TillN':'TILN', 'pEarN':'PROL', 'HkLfN':'LBNN', 'ShkLen':'LBLN', 
              'ShkIntLen':'LBIL', 'CobLen':'EL', 'K1rowN':'CUPR', 'CobDia':'ED', 'pKN':'GE',
              'EarInternodeLen':'EILN', 'tKN':'TGPP', 'tKW':'TGWP', 'GW':'GW'}

chrom_list = list(range(1,11))

#get trait data and covariates
#traitdf = pd.read_csv("G:/My Drive/Teo and landrace/Inbreeding Depression ms/Final data files/Landraces_FL1314_all_data_updated_4398taxa_Fest_2_ed.csv", index_col = 0, na_values = '.') # LFS's path
traitdf = pd.read_csv("Q:/My Drive/Teo and landrace/Inbreeding Depression ms/Final data files/Landraces_FL1314_all_data_updated_4398taxa_Fest_2_ed.csv", index_col = 0, na_values = '.') # Jim's path
traitdf = traitdf.rename(columns = {'F.est':'F_est', 'GBSMom':'Mother', 'GBSDad': 'Father'})
traitdf = traitdf.rename(columns = trait_dict)
#print(traitdf.shape)

parent_list = list(set(traitdf["Mother"].append(traitdf["Father"])))

#==================================================================================================================>>>>>>
#Create dummy variables for parent effects on each progeny
#GCA_columns = []
#for parent in parent_list:
#    col_name = "GCA_" + parent 
#    GCA_columns.append(col_name)        
#    traitdf[col_name] = 0  
#    traitdf.loc[(traitdf.Mother == parent), col_name] = traitdf.loc[(traitdf.Mother == parent), col_name] + 1
#    traitdf.loc[(traitdf.Father == parent), col_name] = traitdf.loc[(traitdf.Father == parent), col_name] + 1

#make a string of all the elements of GCA_columns to use in model formula
#GCA_columns_str = " + ".join(GCA_columns)
#==================================================================================================================>>>>>>

#here's an example of how to make a string that represents 50 PC variables that are named in the data frame as PC1, PC2, ... PC50
PCnums = list(range(1,51))
PCstring = " + PC".join(map(str, PCnums))
PCstring = "PC" + PCstring
#print(PCstring)

#Get the haplotype calls for all marker-parent-haplotype combinations
haps = [] #empty list to store the data frames of hap calls from each chromosome
for chrom in chrom_list:

    #with open("C:/Users/lfsamayo/Documents/Downstairs/LR/Rare_Allele_Scan/LR_rare_allele_hap_calls_chr" + str(chrom) + ".csv", 'r') as f:
    with open("LR_rare_allele_hap_calls_chr" + str(chrom) + ".csv", 'r') as f:
        first_line = f.readline()

    ncols = first_line.count(',') + 1  
    #hapsi = pd.read_csv("C:/Users/lfsamayo/Documents/Downstairs/LR/Rare_Allele_Scan/LR_rare_allele_hap_calls_chr" + str(chrom) + ".csv", index_col = None, usecols = list(range(1,ncols)), dtype = float)
    hapsi = pd.read_csv("LR_rare_allele_hap_calls_chr" + str(chrom) + ".csv", index_col = None, usecols = list(range(1,ncols)), dtype = float)

    #now separately pull in the index
    #hapsi.index = pd.read_csv("C:/Users/lfsamayo/Documents/Downstairs/LR/Rare_Allele_Scan/LR_rare_allele_hap_calls_chr" + str(chrom) + ".csv", index_col = 0, usecols = [0,1]).index #read in one extra column then slice it out!
    hapsi.index = pd.read_csv("LR_rare_allele_hap_calls_chr" + str(chrom) + ".csv", index_col = 0, usecols = [0,1]).index #read in one extra column then slice it out!
    haps.append(hapsi)
    
hap_calls = pd.concat(haps, axis = 1)

##########################################################################################
#DEFINE FUNCTIONS
##########################################################################################

#Define a function to get dominance coefficients for a marker
def dominance(marker):
    return marker%2
    
def check_geno_class(df, marker):
    '''
    A helper function to take a marker and check if it has all three geno classes
    Returns True or False
    '''
    return(len(df[marker].value_counts()) == 3) 

#do forward selection to fit A and D effects for all QTL detected
#starting from base model (covariates plus PCs), try each marker one at a time, adding both A + D effects for the marker
#select marker that reduces residual SS most when added
#check BIC of model, if it improves, rescan the remaining markers and continue iterating
def forward_select(indata, yVar, covars, markers, markers2class):
    '''
    Take in a list of markers at QTL peaks, fit additive + dominance effects simultaneously for each marker one at a time
    Find marker that improves BIC most
    Check BIC of model, if it's worse than previous step, return list of markers in list and final model fit object
    If it's better than previous step, update the model with marker, and check remaining markers one a time
    markers2class is a list of markers within markers list that have only two geno classes. it may be an empty list.
    '''
    markers_to_test = markers[:]
    markers_in_model = []
    
    previous_form = yVar + ' ~ ' + covars
    previous_results = smf.ols(formula = previous_form, data=indata).fit()
    previous_BIC = previous_results.bic
    print("Starting forward selection, base model BIC is:" + str(previous_BIC))
    
    while len(markers_to_test) > 0:
        print("Number of markers in model: " + str(len(markers_in_model)))
        #print("Base formula is updated to:" + previous_form)
        for marker in markers_to_test:
            if marker in markers2class:
                new_form = previous_form + ' + ' + marker #fit only additive effect if we don't have all three geno classes
            else:
                new_form = previous_form + ' + ' + marker + ' + ' + marker + 'D'
            new_results = smf.ols(formula = new_form, data=indata).fit()
    
            #get BIC of model
            new_BIC = new_results.bic
            
            #compare BIC to previous BIC, if better, then current marker is selected as current best and previous versions are all updated
            if new_BIC < previous_BIC:
                print("Found a better model with marker " + marker)
                best_marker = marker
                best_results = new_results
                previous_BIC = new_BIC
        
        #after checking all markers in marker_to_test, check if marker is significant according to threshold,select best marker, add to markers in model and drop from markers_to_test
        if not 'best_marker' in locals():
            print('No QTL could be fit for this trait')
            return([], previous_results, previous_form) #return an empty list and base model fit object
        elif best_marker in markers_in_model:
            #this means best marker in this iteration is same as best marker previously, so no new marker was selected, and we are done
            print("No new markers to add, done with forward selection")
            return(markers_in_model, previous_results, previous_form)
        else:
            print("Best marker in this round is being added to model: " + best_marker)

            markers_in_model.append(best_marker)
            markers_to_test.remove(best_marker)
            if best_marker in markers2class:
                previous_form = previous_form + ' + ' + best_marker  #fit only additive effect if we don't have all three geno classes
            else:
                previous_form = previous_form + ' + ' + best_marker + ' + ' + best_marker + 'D'
            previous_results = best_results
            #print('len(markers_to_test) = ' + str(len(markers_to_test)))
            #print('new_BIC = ' + str(new_BIC))
            #print('previous_BIC = ' + str(previous_BIC))
            #return back to top of while loop

    #if we get through this whole process without returning, and there are no markers left to test, it means we fit all markers in model
    return(markers_in_model, previous_results, previous_form)


#define a function to take the final set of markers, drop one at a time to measure the change in model SS compared to full model
#as a way to compute the SS accounted for by that marker
#have to do this the hard way because statsmodels OLS method fails with high colinearity (type III SS are not correct!)
def markerSS(indata, yVar, covars, markers, markers2class):
    '''
    Take in a list of markers from final model selected using forward_selection()
    Fit additive + dominance effects simultaneously for all markers together as full model
    Next fit only additive effects for all markers together, compute the change in resid SS as the SS due to dominance effects combined
    Next fit only dominance effects for all markers together, compute the change in resid SS as the SS due to additive effects combined
    Then, for each marker one at a time, drop the marker's A effect, then it's D effect, then A + D together,
    each time record the change in residual SS compared to the full model, and report this as the SSA, SSD, and SSA+D for each marker
    Return a data frame with the combined and individual SS for A, D, and A+ D effects
    markers2class is a list of markers within markers list that have only two geno classes
    Do not include dominance tests for these markers
    '''
    markers2classUpdated = markers2class[:] #make a shallow copy that we can update with any additional markers that 'effectively' have only one estimable effect
    marker_list_A = " + ".join(markers) 
    marker_list_D = "D + ".join([m for m in markers if m not in markers2class]) + "D" #exclude markers in markers2class
    marker_list_AD = marker_list_A + " + " + marker_list_D
    full_form = yVar + ' ~ ' + covars + ' + ' + marker_list_AD
    full_results = smf.ols(formula = full_form, data=indata).fit()
    full_ssr = full_results.ssr #residual SS from full model
    full_sstot = full_results.ess + full_ssr
    full_r2 = full_results.rsquared
    print("Full model total SS: " + str(full_sstot))
    print("Full model model SS: " + str(full_results.ess))
    print("Full model residual SS: " + str(full_results.ssr))
    #print("summary of full model:")
    #print(full_results.summary())
    base_results = smf.ols(formula = yVar + ' ~ ' + covars, data=indata).fit()
    #print("summary of base model, no markers:")
    #print(base_results.summary())

    r2_D = full_r2 - smf.ols(formula = yVar + ' ~ ' + covars + ' + ' + marker_list_A, data=indata).fit().rsquared #SSD is diff in ssr when D effects are dropped from model
    r2_A = full_r2 - smf.ols(formula = yVar + ' ~ ' + covars + ' + ' + marker_list_D, data=indata).fit().rsquared
    r2_AD = full_r2 - base_results.rsquared

    dfAD = pd.DataFrame({'Trait': yVar, 'Marker': 'All markers', 'r2_A':r2_A, 'r2_D':r2_D, 'r2_AD': r2_AD}, index = ['ALL_'+yVar])

    results_list = []    
    for marker in markers:
        print("Testing marker " + marker)
        others = [x for x in markers if x != marker] #gets a list with all the other markers except current one
        new_list_A = " + ".join(others) 
        #print("new_list_A: " + new_list_A) 
        new_list_D = "D + ".join([m for m in others if m not in markers2class]) + "D"
        #print("new_list_D: " + new_list_D) 
        new_list_AD = new_list_A + " + " + new_list_D
        #print("new_list_AD: " + new_list_AD) 

        #print("Drop additive effect of marker, new residual SS = " + str(smf.ols(formula = trait + ' ~ ' + covars + ' + ' + new_list_D, data=indata).fit().ssr)) 
        #fit all markers EXCEPT test marker additive effects (in new_list_A) + all marker dominance effects (in marker_list_D)
        #the drop in R2 compared to full final model measures the additive effect of test marker alone
        if marker not in markers2class: #only compute A and D variances separately for markers with all 3 classes
            dropA_fit = smf.ols(formula = yVar + ' ~ ' + covars + ' + ' + marker_list_D + ' + ' + new_list_A, data=indata).fit()
            
            assert full_results.rsquared - dropA_fit.rsquared >= -0.0001 #give a little tolerance to rounding error when difference is actually zero
            assert dropA_fit.df_resid >= full_results.df_resid
            markerA_r2 =  full_results.rsquared - dropA_fit.rsquared 
            if dropA_fit.df_resid == full_results.df_resid: #this condition occurs when A and D effects at a marker are colinear in context of full model, we can only estimate A effect in this case
                markers2classUpdated.append(marker) #UPDATE THE markers2class list!!
                markerA_r2 =  np.nan
                markerD_r2 =  np.nan
            else:
                dropD_fit = smf.ols(formula = yVar + ' ~ ' + covars + ' + ' + marker_list_A + ' + ' + new_list_D, data=indata).fit()
                assert dropD_fit.rsquared < full_results.rsquared
                assert dropD_fit.df_resid == full_results.df_resid + 1
                markerD_r2 =  full_results.rsquared - dropD_fit.rsquared 
                
        else: #for markers with only two classes, we do not report A and D r2 separately
            markerA_r2 = np.nan
            markerD_r2 = np.nan
            
        #fit all marker effects, minus additive and dominance effects of test marker (new_list_AD)
        #the drop in R2 compared to full final model measures the combined additive + dominance effects of test marker
        dropAD_fit = smf.ols(formula = yVar + ' ~ ' + covars + ' + ' + new_list_AD, data=indata).fit()
        assert dropAD_fit.df_resid >= full_results.df_resid + 1 #markers with 1 or more df will pass this check
        assert dropAD_fit.rsquared < full_results.rsquared
        markerAD_r2 =  full_results.rsquared - dropAD_fit.rsquared 
        df = pd.DataFrame({'Trait':yVar,'Marker': marker, 'r2_A':markerA_r2, 'r2_D': markerD_r2, 'r2_AD': markerAD_r2}, index = [marker + "_" + trait])
        #print(df)
        results_list.append(df)

    dfmarks = pd.concat(results_list).append(dfAD)
    return(markers2classUpdated, dfmarks)
        
##########################################################################################
# END DEFINE FUNCTIONS
##########################################################################################

#create the generator for training and test splits.
#this only needs to be done once and can be re-used for different traits
kf = KFold(n_splits=10, shuffle=True, random_state=1)

#make a string of all the covariates to include in all models
#covar_str = 'C(Year) + Tire_track:C(Year) + Shading:C(Year) + F_est:C(Year) + ' + GCA_columns_str
covar_str_noPCs = 'C(Year) + Tire_track:C(Year) + Shading:C(Year) + F_est:C(Year) + x1_1 + x1_2 + x1_3 + x1_4 + x2_1 + x2_2 + x2_3 + x2_4 +    y1_1 +    y1_2 +    y1_3 +    y1_4 +    y2_1 +    y2_2 +    y2_3 +    y2_4 '  
covar_str = covar_str_noPCs + ' + ' + PCstring 

##########################################################################################
#loop through the traits
##########################################################################################
cv_results_list = [] #empty list to collect results from each CV fold and trait

for trait in trait_list:
    print("Current trait: " + trait)

    #Get the QTL peaks from the one-at-a-time marker scan
    peaks = pd.read_csv("QTL peaks LR " + trait + ".csv")

    #Make a list of the marker_parent_hap combinations that were declared as peaks
    peak_list = list(set(peaks.marker))
    
    #filter the hap calls data frame to only include markers that are QTL peaks
    hap_calls_filt = hap_calls.iloc[:,hap_calls.columns.isin(peak_list)]    
    
    #replace missing marker calls with zeroes
    hap_calls_filt.fillna(0, inplace = True)

    #merge hap_calls with trait data frame    
    hap_trait =  traitdf.merge(hap_calls_filt, how='left', left_index=True, right_index=True)

    #make a list of any markers that have only two genotypic classes
    markers_2class = [m for m in peak_list if not check_geno_class(hap_trait, m) ]
  
    #create the dominance coefficients for each marker with 3 geno classes and name them as markerD
    for marker in peak_list:
        if marker not in markers_2class:
            hap_trait[marker + "D"] = dominance(hap_trait[marker])


    #TESTING
    ####################################################################################
    for train_index, test_index in kf.split(hap_trait):
    
        train_set = hap_trait.iloc[train_index,]
        test_set = hap_trait.iloc[test_index,]
        null_form = trait + ' ~ ' + covar_str_noPCs
        null_results = smf.ols(formula = null_form, data=train_set).fit()
        preds_null = null_results.predict(test_set)
        pc_form = trait + ' ~ ' + covar_str
        pc_results = smf.ols(formula = pc_form,  data=train_set).fit()
        preds_pc = pc_results.predict(test_set)
        select_markers, final_fit, final_form = forward_select(train_set, trait, covar_str, peak_list, markers_2class)        
        
        #predict values in the validation set
        preds_full = final_fit.predict(test_set)
        
        #predict values in valiation set with QTLs but without the PCs
        qtl_form = re.sub("\+ PC[0-9]+", "", final_form)
        qtl_results = smf.ols(formula = qtl_form,  data=train_set).fit()
        preds_qtl = qtl_results.predict(test_set)
    
        #check correlations between observed and predicted values
        r_null = test_set[trait].corr(preds_null) #environment covariates only
        r_pc = test_set[trait].corr(preds_pc) #environment covariates + PCs
        r_qtl = test_set[trait].corr(preds_qtl) #environment covariates + qtl
        r_full = test_set[trait].corr(preds_full)  #environment covariates + PCs + qtl
        
        #compute change in correlation when adding qtl to null or pc model
        r_qtl_vs_null = r_qtl - r_null
        r_full_vs_pc = r_full - r_pc
        
        #regress observed values on predicted values
        testdf = pd.DataFrame({'obs':test_set[trait], 'preds_null':preds_null, 'preds_pc':preds_pc, 'preds_qtl':preds_qtl, 'preds_full': preds_full}).dropna()
        
        r2_null = smf.ols(formula = 'obs ~ preds_null', data = testdf).fit().rsquared
        r2_pc = smf.ols(formula = 'obs ~ preds_pc', data = testdf).fit().rsquared
        r2_qtl = smf.ols(formula = 'obs ~ preds_qtl', data = testdf).fit().rsquared       
        r2_full = smf.ols(formula = 'obs ~ preds_full', data = testdf).fit().rsquared

        #compute change in R2 when adding qtl to null or pc model
        r2_qtl_vs_null = r2_qtl - r2_null
        r2_full_vs_pc = r2_full - r2_pc
    
        #combine CV results into a data frame and save to results list
        cv_results = pd.DataFrame({'Trait': trait, 'r_null':r_null, 'r_pc':r_pc, 'r_qtl':r_qtl, 'r_full':r_full, 'r_qtl_vs_null':r_qtl_vs_null, 'r_full_vs_pc':r_full_vs_pc,'r2_null':r2_null, 'r2_pc':r2_pc, 'r2_qtl':r2_qtl, 'r2_full':r2_full, 'r2_qtl_vs_null':r2_qtl_vs_null, 'r2_full_vs_pc':r2_full_vs_pc,}, index = [0])
        
        cv_results_list.append(cv_results)
    

cv_r_and_r2 = pd.concat(cv_results_list)

#compute the means over folds
cv_r_r2_means = cv_r_and_r2.groupby('Trait').mean()

cv_r_r2_means.to_csv("Q:/My Drive/Teo and landrace/Landrace/Rare_Allele_Scan/Cross-validation/LR_RAS_CV_results.csv")
