# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 14:22:00 2019

@author: mjrubino
"""

"""

                        Range-Correlations-from-SB.py
        
    Use this to compile range richness data using range tables in ScienceBase
    and conduct correlation analyses of pairwise taxa richness comparisons and
    an overall richness index.
    
    This uses the sciencebasepy package to access taxa range files:
        National_GAP_Amphibians_Range_Table.txt
        National_GAP_Birds_Range_Table.txt
        National_GAP_Mammals_Range_Table.txt
        National_GAP_Reptiles_Range_Table.txt
    
    Each of these files contain columns with 6 letter species code, the 12-digit
    HUC code, and values for species origin, presence, reproductive use, and
    season in each HUC.
    
    Pairwise comparisons between taxa groups are made and plotted and correlation
    coefficients are calculated using Pearson and Kendall methods. The overall
    richness index is calculated and compared to each taxa group to test that
    group's richness with all other groups' richness.

    Package dependancies:
        sciencebasepy
        pyodbc
        pandas
        numpy
        SciPy
        datetime (for calculating processing time)
        StringIO
        Seaborn
        datetime
        seaborn
        matplotlib


@author: mjrubino
11 July 2019

"""

##############################################################################

import sys, sciencebasepy
import pandas as pd
import numpy as np
from datetime import datetime
from io import StringIO
import matplotlib.pyplot as plt
import matplotlib
# Seaborn for pairplots
import seaborn as sns


starttime = datetime.now()
timestamp = starttime.strftime('%Y-%m-%d')

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#            ++++ Directory Locations ++++
analysisDir = 'C:/Data/USGS Analyses/'
workDir = analysisDir + 'Taxa-Richness-Correlations/'

pd.set_option('display.max_columns', 10)

# Total numbers of species in taxa  # including subspecies
totAmph = 282 #284
totBird = 621 #649
totMamm = 365 #459
totRept = 322 #327


# Make an empty master dataframe
dfMaster = pd.DataFrame()

'''
    Connect to ScienceBase to pull down HUC-12 range tables.
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
    print("="*65)
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
print('\n\n---Calculating Taxa Richness Correlations---\n')
# Run pairwise correlations using Kendall's tau
rcols = ['nAmphibians','nBirds','nMammals','nReptiles']
dfRichCorr = dfMaster[rcols]
# Raw pairwise correlations using Kendall's tau
print("  Calculating Kendalls's tau for pairwise taxa richness comparisons ....")
pw = dfRichCorr.corr(method='kendall')

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

# Spearman's rho correlations of raw richness v. overall richness index
print("  Calculating Spearman's rho with overall richness index ....\n")
spcorrAmph = dfMaster['nAmphibians'].corr(dfMaster['ori-Amph'], method='spearman')
spcorrBird = dfMaster['nBirds'].corr(dfMaster['ori-Bird'], method='spearman')
spcorrMamm = dfMaster['nMammals'].corr(dfMaster['ori-Mamm'], method='spearman')
spcorrRept = dfMaster['nReptiles'].corr(dfMaster['ori-Rept'], method='spearman')

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

'''
    Plot pairwise data using the PairGrid class
    in the Seaborn visualization library 

'''

## First pull out only the pairwise columns from the master dataframe
dfPW = dfMaster[['nAmphibians','nBirds','nMammals','nReptiles']]


# Make a function to calculate Pearson's and Kendall's correlation coefficients
# for pairwise comparisons and make them labels for scatter plots
def corrs(x, y, **kwargs):
    from scipy.stats import kendalltau
    # Calculate Pearson coefficient value
    p = np.corrcoef(x, y)[0][1]
    # Calculate Kendall coefficient values
    t = kendalltau(x, y)[0]
    
    # Make the rho label
    label1 = r'$\rho$ = ' + str(round(p, 2))
    # Make the tau label
    label2 = r'$\tau$ = ' + str(round(t, 2))
    
    # Add each label to the plot
    ax = plt.gca()
    ax.annotate(label1, xy = (0.1, 0.9), xycoords = ax.transAxes)
    ax.annotate(label2, xy = (0.1, 0.8), xycoords = ax.transAxes)

# Plot the pairwise comparisons with PairGrid where:
#  the diagonal of the grid plot are histograms for that taxa's richness within range HUCs
#  the upper portion of the grid plot are kernel density contour plots
#  the lower of the grid plot are scatter plots with a trend line 
#  and Pearson and Kendall correlation coefficents for each comparison

grid = sns.PairGrid(data= dfPW, size = 3)
grid = grid.map_upper(sns.kdeplot,n_levels=30, cmap='jet')
grid = grid.map_diag(plt.hist, bins = 20, edgecolor='k')
grid = grid.map_lower(sns.regplot, marker='.',
                      scatter_kws={"alpha":0.3, "edgecolor":'none', "color":'olivedrab'},
                      line_kws={"color":'black',"alpha":0.5,"lw":1})
grid = grid.map_lower(corrs)

'''
    Now make scatterplots for each taxa group comparison with the
    Overall Richness Index with that group removed (I).
    The non-parametric Spearman's correlation coefficient (rho) is
    calculated for each taxa comparison because the data are not
    normally distributed.

'''
### Make individual taxa subplots on a 2x2 figure ###

fig, axs = plt.subplots(2, 2, figsize=(11, 11))
axs = axs.flatten()


# Add an iterator variable
i=0
for ax, t in zip(axs, taxa):

    # Pull out the taxa specific data from dfMaster
    ncol = 'n{0}'.format(t)
    icol = 'ori-{0}'.format(t[:4])
    dfT = dfMaster[[ncol, icol]]
    
    # Set an axes title
    ax.set_title(t)
    
    # Use the default color palette from Seaborn as variables
    c = sns.color_palette()[i]
    
    # Plot a scatter plot and add a trend line
    # Use the Seaborn colors in the plot with all taxa for individual taxa
    sns.regplot(x=dfT[ncol],y=dfT[icol],marker='.', ax=ax,
               scatter_kws={"alpha":0.3, "edgecolor":'none', "color":c},
               line_kws={"color":'black',"alpha":0.5,"lw":1})

    sr = dfSR.iloc[i,1]
    # Make the tau label
    label = r'Spearman $\rho$ = ' + str(round(sr, 4))
    ax.annotate(label, xy = (0.45, 0.9), xycoords = ax.transAxes)
    
    i+=1


endtime = datetime.now()
delta = endtime - starttime
print("+"*35)
print("Processing time: " + str(delta))
print("+"*35)