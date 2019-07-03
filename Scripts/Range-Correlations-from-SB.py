# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 14:22:00 2019

@author: mjrubino
"""







import sys, sciencebasepy
import pandas as pd
import numpy as np
from datetime import datetime
from io import StringIO

starttime = datetime.now()
timestamp = starttime.strftime('%Y-%m-%d')

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#            ++++ Directory Locations ++++
analysisDir = 'C:/Data/USGS Analyses/'
workDir = analysisDir + 'Range-vs-Habitat/'
tempDir = workDir + 'downloadtemp/'

pd.set_option('display.max_columns', 10)

# Total numbers of species in taxa  # including subspecies
totAmph = 282 #284
totBird = 621 #649
totMamm = 365 #459
totRept = 322 #327


# Make an empty master dataframe
dfMaster = pd.DataFrame()

'''
    Connect to ScienceBase to pull down a HUC12 range tables.
    This uses the ScienceBase item for species habitat maps
    and searches for range table files with species range data in it.
    The range maps item has a unique id (5951527de4b062508e3b1e79).
    If this changes, the code will need to be re-written.

'''
sb = sciencebasepy.SbSession()
rangeItem = sb.get_item("5951527de4b062508e3b1e79")

# Set column types for the range table dataframe
coltypes = {'strUC':str,
            'strHUC12RNG':str,
            'intGapOrigin':int,
            'intGapPres':int,
            'intGapRepro':int,
            'intGapSeas':int}

#taxa = ['Amphibians','Reptiles'] 
taxa = ['Amphibians','Birds','Mammals','Reptiles']
for t in taxa:
    print("\nWorking on Range Table for " + t + " ....")
    for file in rangeItem["files"]:
        # Search for the file name pattern in the range item files dictionary
        if file['name'].startswith('National_GAP_{0}'.format(t)):
            try:
                
                dfRangeTable = pd.read_csv(StringIO(sb.get(file['url'])), dtype=coltypes)
                # Select only those records that are for full species
                dfRT2 = dfRangeTable[dfRangeTable['strUC'].str.endswith('x')]
                # Select only known, possibly, or potentially present;
                #             year-round, winter, or summer seasons
                select={'intGapPres':[1,2,3], 'intGapSeas':[1,3,4]}
                dfS1 = dfRT2[dfRT2[list(select)].isin(select).all(axis=1)]
                
                # Group by HUC12 and get the species count for each
                dfS2 = dfS1.groupby(['strHUC12RNG'])['strUC'].count()
                dfCount = pd.DataFrame(data=dfS2)
                dfCount = dfCount.rename(columns={'strUC':'n' + t})
                
                # Concatenate to the master dataframe
                print('   Concatenating with master dataframe ....')
                dfMaster = pd.concat([dfMaster,dfCount], axis=1, sort=False)
                # Make sure any missing values are 0 and the column type is integer
                dfMaster = dfMaster.replace(np.nan, 0, regex=True)
                dfMaster = dfMaster.astype(int)
                dfMaster.index.name = 'HUC12'
            
                print("="*65)
            except:
                print('!!! Could not Find Range Table Item File. Exiting !!!')
                sys.exit()



'''

    Start manipulating the master dataframe to conduct correlation analyses

'''

# Calculate proportional richness and overall richness index
dfMaster['ori-Amph'] = ((dfMaster['nBirds']/totBird)
                        +(dfMaster['nMammals']/totMamm)
                        +(dfMaster['nReptiles']/totRept))/3
dfMaster['ori-Bird'] = ((dfMaster['nAmphibians']/totAmph)
                        +(dfMaster['nMammals']/totMamm)
                        +(dfMaster['nReptiles']/totRept))/3
dfMaster['ori-Mamm'] = ((dfMaster['nAmphibians']/totAmph)
                        +(dfMaster['nBirds']/totBird)
                        +(dfMaster['nReptiles']/totRept))/3
dfMaster['ori-Rept'] = ((dfMaster['nAmphibians']/totAmph)
                        +(dfMaster['nBirds']/totBird)
                        +(dfMaster['nMammals']/totMamm))/3

# run correlations
rcols = ['nAmphibians','nBirds','nMammals','nReptiles']
dfRichCorr = dfMaster[rcols]
# Raw pairwise correlations using Kendall's tau
pw = dfRichCorr.corr(method='kendall')

# Spearman's rho correlations of raw richness v. overall richness index
print("  Calculating Spearman's rho with overall richness index ....")
spcorrAmph = dfMaster['nAmphibians'].corr(dfMaster['ori-Amph'], method='spearman')
spcorrBird = dfMaster['nBirds'].corr(dfMaster['ori-Bird'], method='spearman')
spcorrMamm = dfMaster['nMammals'].corr(dfMaster['ori-Mamm'], method='spearman')
spcorrRept = dfMaster['nReptiles'].corr(dfMaster['ori-Rept'], method='spearman')

'''
# Pull the data of raw correlation coefficients from
# the correlation dataframe element by element
data = [pw.at['nBirds','nAmphibians'],
       pw.at['nMammals','nAmphibians'],
       pw.at['nMammals','nBirds'],
       pw.at['nReptiles','nAmphibians'],
       pw.at['nReptiles','nBirds'],
       pw.at['nReptiles','nMammals']]

# Make a static column of taxa pairs and insert
# correlation data and HUC number
pwdict = {'TaxaPair':['B-A','M-A','M-B','R-A','R-B','R-M'],
          'KendallsTau':data}
# Make the dataframe
dfCC = pd.DataFrame(pwdict)
# Reorder columns
dfCC = dfCC[['TaxaPair', 'KendallsTau']]

# Export to a table
print("  Exporting correlations to tables ....")
dfCC.to_csv(workDir + huc + "_PairwiseCorrelations.txt")
'''

# Pull Spearman's correlation data
spdata = [spcorrAmph, spcorrBird, spcorrMamm, spcorrRept]
# Make a static column of taxa and insert
# Spearman's correlation data and HUC number
spdict = {'Taxon':['Amphibians','Birds','Mammals','Reptiles'],
          'SpearmansRho':spdata}
# Make the dataframe
dfSR = pd.DataFrame(spdict)
# Reorder columns
dfSR = dfSR[['Taxon','SpearmansRho']]
#dfSR = dfSR.set_index(['HUC'])
# Export to a table
#dfSR.to_csv(corrDir + huc + "_SpearmanCorrelations.txt")




endtime = datetime.now()
delta = endtime - starttime
print("+"*35)
print("Processing time: " + str(delta))
print("+"*35)