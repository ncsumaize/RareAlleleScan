# -*- coding: utf-8 -*-
'''
A script to take in output of rare allele QTL scan and:
    
    1. group continuous markers above some threshold as potentially within common QTL support intervals
    2. identify the peak marker for each QTL
    3. identify the markers that define the ends of the QTL support interval
    4. output a file that identifies the QTL peaks and support intervals for each combination of trait, parent, and haplotype
'''
import os as os
import pandas as pd

#set working directory
#os.chdir("G:/My Drive/Teo and landrace/Landrace/Rare_Allele_Scan/final_outputs") # LFS's path
os.chdir("Q:/My Drive/Teo and landrace/Landrace/Rare_Allele_Scan/final_outputs") # Jim's path

#define trait names
#old names
#trait_list= ['DTA', 'DTS', 'PltHt','LfLen', 'LfWth', 'TillN', 'pEarN', 'HkLfN', 'ShkLen', 'ShkIntLen', 'EarInternodeLen', 'CobLen', 'K1rowN', 'CobDia', 'pKN', 'tKN', 'tKW', 'pKW50']
#new names
trait_list= ['DTA', 'DTS','PLHT', 'LFLN', 'LFWD', 'TILN', 'PROL', 'LBNN', 'LBLN', 
              'LBIL', 'EL', 'CUPR', 'ED', 'GE', 'EILN', 'TGPP', 'TGWP', 'GW']

#DEFINE threshold values for declaring QTL as global values here
LOD_thresh1 = 4 #A QTL must have LOD threshold greater than this to be a QTL peak
LOD_drop_thresh = 2 #left and right support intervals of QTL are within this LOD threshold of the peak QTL

#get marker position information to get them ordered properly
#marker_info = pd.read_csv("G:/My Drive/Teo and landrace/Inbreeding Depression ms/Final data files/LR_genetic_map.csv") # LFS's path
marker_info = pd.read_csv("Q:/My Drive/Teo and landrace/Inbreeding Depression ms/Final data files/LR_genetic_map.csv") # Jim's path
marker_info = marker_info.rename(columns = {'marker': 'SNP', 'AGPv4_Pos': 'pos_Agpv4'})

#Make a function to split the 'marker' information into marker/parent/haploty
#for LR, pattern is 'S#_####_###_#hap#'
def marker_split(x):
    split1 = x.split("hap")
    hap = split1[-1]
    rest = split1[0]
    split2 = rest.split("_")
    marker = split2[0] + "_" + split2[1]
    parent = split2[2] + "_" + split2[3]
    return(pd.Series({'SNP': marker, 'Parent' : parent, 'Hap' : hap}))



#define a function to identify QTL peaks and left and right support interval boundaries
#this will operate on a data frame containing QTL scan results from one parent and haplotype and chromosome
#it will find the peak (if any) QTL on the chromosome, identify the left and right support boundaries,
#then recursively process sets of markers to left and right of the QTL boundaries on the same chromosome
def recursive_qtl(df, Trait, results_list):
    #dfQ is a subset of the scan_sub2 data frame containing only markers in current QTL group: could be a whole chromosome, or a contiguous set of markers that are part of a chromosome
    #get the row with maximum LOD score in this set
    #make a deep copy of the data frame originally passed in, so we don't mutate it
    dfQ = df.copy()
    peakIndex = dfQ.LOD.idxmax()
    scan_qt_max = dfQ.loc[peakIndex,:]
    #get the integer row number values of the adjacent left and right markers
    #I name these as Index, but they are in fact integers to use with .iloc[]
    peakInteger = dfQ.index.get_loc(peakIndex)
    leftIndex = peakInteger - 1
    rightIndex = leftIndex + 2
    
    #BASE CASE is that there is no marker passing qtl peak, end recursion and return None
    if scan_qt_max.LOD < LOD_thresh1: return None
    
    #compute LOD drop for every marker within the QTL region
    dfQ["Lod_drop"] = scan_qt_max.LOD - dfQ.LOD
    
    #get the left-hand marker cM position for the peak QTL
    #edge case when peak is at left end of marker group, peak is also left marker
    if peakInteger == 0: 
        leftSide = None
        leftMarker = scan_qt_max 
    #edge case when peak is one marker from left end of chromosome
    elif peakInteger == 1: 
        leftSide = None
        leftMarker = dfQ.iloc[leftIndex,:]
            
    else:
        leftSide = dfQ.iloc[0:peakInteger,:] 
        leftSide = leftSide.loc[leftSide.Lod_drop < LOD_drop_thresh,:] #keep only markers that exceed the LOD drop threshold

        #case of no markers within threshold on left side, extend to one marker left of peak
        if leftSide.empty:
            leftMarker = dfQ.iloc[leftIndex,:]
            leftSide = dfQ.iloc[0:leftIndex,:] #repopulate leftSide with all left markers even if below LOD drop, so they can be recursively processed
            if leftSide.empty: leftSide = None
            
        else:        
            leftSide["gap"] = leftSide.cM.diff(periods = 1)
            leftSide.gap.fillna(value = 1, inplace = True) #replace na at start or end with 1
           #identify the left-hand marker with gap > 1 and highest cM value
            if sum(leftSide.gap > 1) == 0: #case of no gaps to left-most marker within LOD drop threshold, QTL interval ends one marker to left of left-most marker within this group
                if dfQ.index[0] == leftSide.index[0]: #case where left marker inside current set is also first marker of initial QTL group, we cannot extend any farther left
                    leftMarker = leftSide.iloc[0,:]
                    leftSide = None
                else: #case where left side support interval does not include all the potential left side markers we started with, they are separate by a gap but need to be recursively processed
                    leftSide = dfQ.loc[dfQ.index[0]:leftSide.index[0],:]
                    leftMarker = leftSide.iloc[-2,:] #leftMarker is one marker to LEFT of last marker inside group
                    leftSide.drop(leftSide.tail(2).index,inplace=True) #drop the rightmost 2 markers , are inside the supporrt interval of current QTL
                    if leftSide.empty: leftSide = None
            else: #if there are gaps, the right-most marker before last gap is adjacent to left-hand marker, other markers get recursed as possible second QTL
                leftSideGaps = leftSide.loc[leftSide.gap > 1, :]
                leftSide = dfQ.loc[dfQ.index[0]:leftSideGaps.cM.idxmax(),:]
                leftMarker = leftSide.iloc[-2,:]
                leftSide.drop(leftSide.tail(2).index,inplace=True) #drop the rightmost2 markers, which inside support interval of current QTL
                if leftSide.empty: leftSide = None  
                
    #do the same for markers to right of peak
    #edge case when peak is at right end of marker group, peak is also right marker
    if peakIndex == dfQ.index[-1]: 
        rightSide = None
        rightMarker = scan_qt_max 
    #edge case when peak is one marker from right end of chromosome
    elif peakIndex == dfQ.index[-2]: 
        rightSide = None
        rightMarker = dfQ.iloc[rightIndex,:]
        
    else:
        rightSide = dfQ.iloc[rightIndex:,:] 
        rightSide = rightSide.loc[rightSide.Lod_drop < LOD_drop_thresh,:] #keep only markers within range of the QTL peak by LOD drop

        #case of no markers within threshold on right side, extend to one marker right of peak
        if rightSide.empty:
            rightMarker = dfQ.iloc[rightIndex,:]
            rightSide = dfQ.iloc[(rightIndex+1):,:] #repopulate rightSide with all right markers even if below LOD drop, so they can be recursively processed
            if rightSide.empty: rightSide = None
        else:        
            rightSide["gap"] = abs(rightSide.cM.diff(periods = -1))
            rightSide.gap.fillna(value = 1, inplace = True) #replace na at start or end with 1
            #identify the right-hand marker with gap > 1 and lowest cM value
            if sum(rightSide.gap > 1) == 0: #case of no gaps to right-most marker within LOD threshold of peak, QTL interval ends at right-most marker in this remaining group
                if dfQ.index[-1] == rightSide.index[-1]: #case where right marker is also right marker of initial QTL group
                   rightMarker = rightSide.iloc[-1,:]
                   rightSide = None
                else: #case where right side support interval does not include all the potential right side markers we started with, they are separated by a gap but need to be recursively processed
                    rightSide = dfQ.loc[rightSide.index[-1]:,:]
                    rightMarker = rightSide.iloc[1,:] #rightMarker is one marker to RIGHT of last marker inside group
                    rightSide.drop(rightSide.head(2).index,inplace=True) #drop the two leftmost markers, which inside support interval of current QTL
                    if rightSide.empty: rightSide = None
            else: #if there are gaps, the right-most marker before first gap is adjacent to right-hand marker, other markers get recursed as possible second QTL
                rightSideGaps = rightSide.loc[rightSide.gap > 1, :]
                rightSide = dfQ.loc[rightSideGaps.cM.idxmin():,:]
                rightMarker = rightSide.iloc[1,:]
                rightSide.drop(rightSide.head(2).index,inplace=True) #drop the leftmost marker, which is actually the right-flanking marker of current QTL
                if rightSide.empty: rightSide = None 
    #merge left and right marker info with the peak, then append result to result2 list
    QTL_starts = leftMarker[['Parent', 'Hap', 'chrom', 'SNP', 'cM', 'pos_Agpv4']].to_frame().T #convert to a data frame and transpose
    QTL_starts.rename(columns = {'SNP': 'SNP_start', 'cM': 'cM_start', 'pos_Agpv4': 'Pos_start'}, inplace = True)
    
    QTL_ends = rightMarker[['Parent', 'Hap', 'chrom', 'SNP', 'cM', 'pos_Agpv4']].to_frame().T #convert to a data frame and transpose
    QTL_ends.rename(columns = {'SNP': 'SNP_end', 'cM': 'cM_end', 'pos_Agpv4': 'Pos_end'}, inplace = True)
    
    QTL_final = scan_qt_max.to_frame().T #get the peak QTL as a one-row data frame
    QTL_final = QTL_final.merge(QTL_starts, how = 'inner', on = ['Parent', 'Hap', 'chrom'])
    QTL_final = QTL_final.merge(QTL_ends, how = 'inner', on = ['Parent', 'Hap', 'chrom'])

    QTL_final["Trait"] = Trait
    #print(QTL_final)
    
    results_list.append(QTL_final)
        
    #if any left side markers, recurse on them
    if leftSide is not None:
        recursive_qtl(leftSide, Trait, results_list)        
    #if any right side markers, recurse on them
    if rightSide is not None:
        recursive_qtl(rightSide, Trait, results_list)        
    
    #end recursive_qtl()

def pickPeaks(trait):
    
    scan = pd.read_csv("LR_Rare_allele_scan_output_" + trait + "_chr_1.csv" , index_col=0)
    #print(scan.shape)
    #print(scan.head())
    for chr in list(range(2,11)):
        scani = pd.read_csv("LR_Rare_allele_scan_output_" + trait + "_chr_" + str(chr) + ".csv" , index_col=0)
        scan = pd.concat([scan, scani], sort = True)
    
    #Apply the splitting function to the 'marker' column of scan output, returning three new columns
    hap_info = scan.marker.apply(marker_split)
    
    scan2 = scan.merge(hap_info, left_index = True, right_index = True)
    
    #merge the QTL results with AGPv4 position info and order the markers by genetic map position
    scan3 = marker_info.merge(scan2, how = "right") 
    scan3 = scan3.sort_values(by = ["chrom", "cM"]) 
    
    #we seem to have some markers with no map positions in this version, will fix in future versions
    #for now, just drop them
    scan3 = scan3.loc[pd.notnull(scan3['cM'])]
    scan4 = scan3
    
    #replace any NaN LOD values with zeroes
    scan4.fillna({'LOD':0}, inplace = True)
    
                                        
    parent_list = list(set(scan4.Parent))
    
    results = [] #for each trait make a new results list
    tested = []
    for parent in parent_list:
        for hap in ['1','2']:
            scan_sub = scan4.loc[(scan4.Parent == parent) & (scan4.Hap == hap),:]
            chroms = scan_sub["chrom"].unique() #get a vector with unique chrom numbers for this parent/haplotype
            for chromo in chroms:
                scan_sub2 = scan_sub.loc[scan_sub.chrom == chromo,:]
                print("Parent, Hap, Chrom:")
                print(parent)
                print(hap)
                print(chromo)
                tested.append(parent + "_" + hap + "_" + str(chromo))
                recursive_qtl(df = scan_sub2, Trait = trait, results_list = results)
    
    currentTraitQTL = pd.concat(results, sort = True)
    if 'Lod_drop' in currentTraitQTL.columns:
        currentTraitQTL.drop(columns = ['Lod_drop'], inplace = True)
    currentTraitQTL.to_csv("QTL peaks LR " + trait + ".csv", index = False)

for tr in trait_list: pickPeaks(tr)
#pickPeaks("PLHT")

# This script produces one new data files as output per trait:
# QTL peaks LR TRAIT.csv