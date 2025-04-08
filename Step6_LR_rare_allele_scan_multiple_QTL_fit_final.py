# -*- coding: utf-8 -*-
'''
A script to read in QTL peak positions for one trait

Fit all QTL simultaneously in a linear model with all covariates

Rank markers based on combined additive + dominance SS

If lowest ranking marker has combined F test p-value > 0.001, drop that marker and re-run model

Iterate until lowest ranking marker has combined A + D p-value < 0.001

Make a data frame with marker effects and type 1 and 3 SS estimated from final model and write to a csv file
'''
import os as os
import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
#import time as time
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib

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
            return([], previous_results.fit()) #return an empty list and base model fit object
        elif best_marker in markers_in_model:
            #this means best marker in this iteration is same as best marker previously, so no new marker was selected, and we are done
            print("No new markers to add, done with forward selection")
            return(markers_in_model, previous_results)
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
    return(markers_in_model, previous_results)


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

##########################################################################################
#loop through the traits
##########################################################################################
marker_summary_all_traits_list = [] #empty list to collect marker effects from each trait
peak_info_all_traits_list = [] #empty list to retain info on QTL intervals for each trait
summary_all_traits_list = [] #empty list to collect summary results from each trait~\Google Drive\Teo and landrace\Teo\Rare Allele Scan

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

    #make a string of all the covariates to include in all models
    #covar_str = 'C(Year) + Tire_track:C(Year) + Shading:C(Year) + F_est:C(Year) + ' + GCA_columns_str
    covar_str = 'C(Year) + Tire_track:C(Year) + Shading:C(Year) + F_est:C(Year) + x1_1 + x1_2 + x1_3 + x1_4 + x2_1 + x2_2 + x2_3 + x2_4 +    y1_1 +    y1_2 +    y1_3 +    y1_4 +    y2_1 +    y2_2 +    y2_3 +    y2_4 + ' + PCstring 
    covar_str_noPCs = 'C(Year) + Tire_track:C(Year) + Shading:C(Year) + F_est:C(Year) + x1_1 + x1_2 + x1_3 + x1_4 + x2_1 + x2_2 + x2_3 + x2_4 +    y1_1 +    y1_2 +    y1_3 +    y1_4 +    y2_1 +    y2_2 +    y2_3 +    y2_4 '  

    #test the null model with no QTL and no PCs but all other non-genetic ovariates
    null_form = trait + ' ~ ' + covar_str_noPCs
    null_results = smf.ols(formula = null_form, data=hap_trait).fit()
    #print(null_results.summary())

    #test the model with covariates plus PCs but no QTL
    pc_form = trait + ' ~ ' + covar_str
    pc_results = smf.ols(formula = pc_form, data=hap_trait).fit()
    
    #run the model selection function    
    select_markers, final_fit = forward_select(hap_trait, trait, covar_str, peak_list, markers_2class)        

    print("Number of markers tested: " + str(len(peak_list)))
    print("Number of markers in final model: " + str(len(select_markers)))
    print("R2 of null model, no PCs, no QTL: " + str(null_results.rsquared))
    print("R2 of PC model, no QTL: " + str(pc_results.rsquared))
    print("R2 of PC model plus markers: " + str(final_fit.rsquared))
    print("Variance accounted for by PCs: " + str(pc_results.rsquared - null_results.rsquared))
    print("Variance accounted for by markers: " + str(final_fit.rsquared - pc_results.rsquared))
    
    if len(select_markers) > 0:
        #get the combined and individual marker Sums of Squares
        markers_2class_updated, markerSSdf = markerSS(hap_trait, trait, covar_str, select_markers, markers_2class)
        #print(markerSSdf)    
        R2_A = round(markerSSdf.loc['ALL_'+trait, 'r2_A'],3)
        R2_D = round(markerSSdf.loc['ALL_'+trait, 'r2_D'],3)
        R2_AD = round(markerSSdf.loc['ALL_'+trait, 'r2_AD'],3)
        
        print('R^2 for all selected additive + dominance effects combined = ' + str(R2_AD))
        print('R^2 for all selected additive effects combined = ' + str(R2_A))
        print('R^2 for all selected dominance effects combined = ' + str(R2_D))

        if len(markers_2class_updated) > len(markers_2class): #if we caught some additional cases of a/d colinearity, we need to refit final model, dropping the extra D effects to avoid reporting the colinear effect estimates
            marker_effects_A = " + ".join(select_markers) 
            marker_effects_D =  "D + ".join([m for m in select_markers if m not in markers_2class_updated]) + "D"
            #final_form =  pc_form + '+' + marker_effects_A + '+' + marker_effects_D
            
            # I added this if/then to fix the problem of the empty 'D' list, based on the Jim indication
            if marker_effects_D == 'D': final_form =  pc_form + '+' + marker_effects_A
            else: final_form =  pc_form + '+' + marker_effects_A + '+' + marker_effects_D 
            
            final_fit = smf.ols(formula = final_form, data=hap_trait).fit()
            
        #put together the individual additive and dominance effects for each marker in final fit
        final_ests = pd.DataFrame({'effect':final_fit.params, 'pval':final_fit.pvalues, 'tval':final_fit.tvalues})
    
        #extract the additive and dominance effects separately, and rename them
        add_fx = final_ests.loc[select_markers,]
        add_fx.rename(columns = {'effect':'a_effect', 'pval':'a_pval', 'tval':'a_tval'}, inplace = True)
        #some markers do not have dominance effects, so we need to remove any inded labels for such markers
        dom_indices = [m + "D" for m in select_markers]
        dom_indices = [x for x in dom_indices if x in final_ests.index]
        dom_fx = final_ests.loc[dom_indices,]
        dom_fx.rename(columns = {'effect':'d_effect', 'pval':'d_pval',  'tval':'d_tval'}, inplace = True)
        dom_fx.index = [i[:-1] for i in dom_fx.index] #strip the "D" from index values, so we can merge a and d values on common index
        add_dom_fx = add_fx.merge(dom_fx, left_index = True, right_index = True, how = 'left')
        add_dom_fx.index = [m + "_" + trait for m in add_dom_fx.index] #make the index same as in the markerSSdf so we can merge below 
    
        #count the number of het and homozygous rare allele individuals in data set
        counts = hap_trait[select_markers].apply(pd.value_counts, axis = 0).transpose()
        counts.drop(labels = [0.0], axis = 1, inplace = True)
        counts.rename(columns = {1.0: 'N_het', 2.0: 'N_hom'}, inplace = True)
    
        #combine marker effect estimates from final model with R2 values and het/homoz calls
        markerSSfx = markerSSdf.merge(add_dom_fx, how = "left", left_index = True, right_index = True)
        markerSSfx = markerSSfx.merge(counts, how = "left", left_on = "Marker", right_index = True)
        
        #get the support interval information from the original scans
        peak_info = peaks[['marker', 'SNP', 'Parent', 'Hap', 'chrom', 'pos_Agpv4', 'cM', 'SNP_start', 'cM_start', 'Pos_start', 'SNP_end', 'cM_end', 'Pos_end']]
        peak_info = peak_info.rename(columns = {'marker':'Marker'})
        markerSSfx = markerSSfx.merge(peak_info, how = "left", on = "Marker")
        
        marker_summary_all_traits_list.append(markerSSfx)
    
        #make an overall summary for each trait and combine into a common dataframe
        summdf = pd.DataFrame({'Trait': trait, 'QTL_in_final_mod':len(select_markers), 'R2_null_model': null_results.rsquared, 'R2_PC_model': pc_results.rsquared,'R2_w_QTL': final_fit.rsquared, 'R2_due_to_QTL': R2_AD, 'R2_QTL_Add': R2_A, 'R2_QTL_Dom': R2_D}, index = [trait])

        #plot additive versus dominance effects
        plt.clf() #clear any previous plot
        plt.scatter(markerSSfx.a_effect, markerSSfx.d_effect)
        plt.title("Tuxpeno Additive vs Dominance Effects at Rare QTL for " + trait)
        plt.xlabel("Additive Effect")
        plt.ylabel("Dominance Effect")
        plt.text(0, markerSSfx.d_effect.max(), 'Variance due to rare QTL additive = ' + str(round(R2_A,2)))
        plt.text(0, markerSSfx.d_effect.max()*0.85, 'Variance due to rare QTL dominance = ' + str(round(R2_D,2)))
        plt.savefig("fig_LR_" + trait + " A vs D effects.png")
        plt.show()

        #plot additive versus dominance effects scaled to their standard errors (t-values)
        plt.clf() #clear any previous plot
        plt.scatter(markerSSfx.a_tval, markerSSfx.d_tval)
        plt.title("Tuxpeno Additive vs Dominance Scaled Effects at Rare QTL for " + trait)
        plt.xlabel("Scaled Additive Effect")
        plt.ylabel("Scaled Dominance Effect")
        plt.text(0, markerSSfx.d_effect.max(), 'Variance due to rare QTL additive = ' + str(round(R2_A,2)))
        plt.text(0, markerSSfx.d_effect.max()*0.85, 'Variance due to rare QTL dominance = ' + str(round(R2_D,2)))
        plt.savefig("fig_LR_" + trait + " A vs D effects scaled.png")
        plt.show()
    else: #case where no markers were selected in final model
        summdf = pd.DataFrame({'Trait': trait, 'QTL_in_final_mod':0, 'R2_null_model': null_results.rsquared, 'R2_PC_model': pc_results.rsquared, 'R2_w_QTL':pc_results.rsquared, 'R2_due_to_QTL': 0, 'R2_QTL_Add': 0, 'R2_QTL_Dom': 0}, index = [trait])

    summary_all_traits_list.append(summdf)

##########################################################################################
#END loop through traits
##########################################################################################

marker_summary_all_traits_df = pd.concat(marker_summary_all_traits_list)
#marker_summary_all_traits_df = marker_summary_all_traits_df[['Trait', 'Marker', 'r2_A', 'r2_D', 'r2_AD', 'a_effect',  'a_tval','a_pval', 'd_effect',  'd_tval', 'd_pval']]

#write out the two summary data frames to files
#write to a file, #append without header if it already exists with open(fname, 'a') if needed
fname = "LR_rare_QTL_final_marker_effects.csv"
with open(fname, 'w') as f:
    marker_summary_all_traits_df.to_csv(fname, mode='a', index = False, header=f.tell()==0)
f.close()

summary_all_traits_df = pd.concat(summary_all_traits_list)
#summary_all_traits_df = summary_all_traits_df[['Trait', 'R2_null_model', 'R2_PC_model', 'QTL_in_final_mod', 'R2_w_QTL', 'R2_QTL_Add', 'R2_QTL_Dom']]

#write to a file, append without header if it already exists with open(fname, 'a') if needed
fname2 = "LR_rare_QTL_summaries.csv"
with open(fname2, 'w') as f:
    summary_all_traits_df.to_csv(fname2, mode='a', index = False, header=f.tell()==0)
f.close()
    