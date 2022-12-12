from os import listdir
from os.path import join, isfile, isdir, split
import os
import pandas as pd
import configparser

# ------------------------------------------------------------------------------
# ORGANIZATIONAL FUNCTIONS
# ------------------------------------------------------------------------------

class pathParser():
    def __init__(self, filename = 'paths.ini'):

        folder = os.path.abspath(os.pardir)
        inifile = os.path.join(folder, filename)
        config = configparser.ConfigParser()
        config.read(inifile)

        self.alldata= config['FOLDERS']['resultsFolder']
        self.genes = config['FILES']['genes']
        self.manifest = config['FILES']['manifest']
        self.structures = config['FILES']['structures']

        return 

def customCortexOrder():
    """
    Returns a list of cortical areas IDs in a custom order

        Parameters
        ----------
        None

        Returns
        -------
        customOrder:list
            List of cortical area IDs in a curated order so that areas with similar
            neurobiological function are grouped together.
    """
    customOrder = [
        184,985,993,353,329,337,345,369,361,182305689,378,              # SS and motor areas
        402,394,409,385,425,533,312782574,312782628,312782546,417,      # Visual Areas
        1011,1002,1027,1018,                                            # Auditory
        44,972,731,                                                     # Media Prefrontal
        723,746,39,48,894,879,886,                                      # Medial Association
        1057,677,104,111,119,541,895,922,918,926,                       # Lateral + ento
    ]
    return customOrder

def customMidCoarseSubdivision(kind:str='id'):
    if kind == 'id':
        midCoarse = [184,500,453,1057,677,247,669,31,972,44,
            714,95,254,22,541,922,895,698,375,726,909,843,1037,1084,
            502,484682470,589508447,703,485,493,275,278,818,835,826,809,864,856,
            157,141,467,290,10671,339,323,348,1132,987,1117,386,370,379,528,519
        ]
    elif kind == 'acronym':
        midCoarse = [
            'FRP', 'MO','SS', 'GU','VISC', 'AUD', 'VIS', 'ACA', 'PL', 'ILA','ORB', 'AI','RSP','PTLp','TEa','PERI','ECT',
            'OLF',
            'CA', 'DG','ENT','PAR','POST','PRE','SUB','ProS','HATA',
            'CTXsp',
            'STRd','STRv','LSX','sAMY',
            'PALd','PALv','PALm','PALc',
            'DORsm','DORpm',
            'PVZ','PVR','MEZ','LZ','ME',
            'MBsen','MBmot','MBsta',
            'P-sen','P-mot','P-sat',
            'MY-sen','MY-mot','MY-sat',
            'CBX','CBN'
        ]
    return midCoarse

def corticalAreasGroups(toDataFrame:bool=False):
    """
    Returns a dict where each cortical area acronym is associated to one of the
    cortical subnetworks defined in (Zingg et al. 2014, Cell). Subnetworks are: 
    motorSS, lateral, audioVisual, mPrefrontal, mAssociation.

        Parameters
        ----------
        None

        Returns
        -------
        areaGroups:dict
            Dict, where keys are area acronyms and values are the corresponding 
            cortical subnetwork
    """
    areaGroups = {
        'FRP':'other',
        'MOp':'motorSS',
        'MOs':'motorSS',
        'SSp-n':'motorSS',
        'SSp-bfd':'motorSS',
        'SSp-ll':'motorSS',
        'SSp-m':'motorSS',
        'SSp-ul':'motorSS',
        'SSp-tr':'motorSS',
        'SSp-un':'motorSS',
        'SSs':'motorSS',
        'GU':'lateral',
        'VISC':'lateral',
        'AUDd':'audioVisual',
        'AUDp':'audioVisual',
        'AUDpo':'audioVisual',
        'AUDv':'audioVisual',
        'VISal':'audioVisual',
        'VISam':'audioVisual',
        'VISl':'audioVisual',
        'VISp':'audioVisual',
        'VISpl':'audioVisual',
        'VISpm':'audioVisual',
        'VISli':'audioVisual',
        'VISpor':'audioVisual',
        'ACAd':'mAssociation',
        'ACAv':'mAssociation',
        'PL':'mPrefrontal',
        'ILA':'mPrefrontal',
        'ORBl':'mAssociation',
        'ORBm':'mPrefrontal',
        'ORBvl':'mAssociation',
        'AId':'lateral',
        'AIp':'lateral',
        'AIv':'lateral',
        'RSPagl':'mAssociation',
        'RSPd':'mAssociation',
        'RSPv':'mAssociation',
        'VISa':'mAssociation',      #'audioVisual',
        'VISrl':'mAssociation',     #'audioVisual',
        'TEa':'lateral',
        'PERI':'lateral',
        'ECT':'lateral',
        'ENTl':'other',
        'ENTm':'other'
    }

    if toDataFrame:
        areaGroups = pd.DataFrame.from_dict(areaGroups,orient='index').rename(columns={0:'function'})            

    return areaGroups

def customSensoryAreaIds(oldAtlasNumbers:bool=False):
    """
    Returns a list of IDs of cortical areas that are specific for sensory-related
    portions of the cortex (Visual, Auditory and Somatosensory)

    Parameters
    ----------
    oldAtlasNumbers:bool = False
        If true it returns the old ID 22 used in the CCFv2 instead of the two visual
        areas 417 and 312782546

    Returns
    -------
    sensoryIds:list[int]
        List of the IDs
    """
    
    if oldAtlasNumbers:
        sensoryIds = [
            353,329,337,345,369,361,182305689,378,              # Somatosensory
            402,394,409,385,425,533,312782574,312782628,22,     # Visual
            1011,1002,1027,1018                                 # Auditory
        ]
    else:
        sensoryIds = [
            353,329,337,345,369,361,182305689,378,                      # Somatosensory
            402,394,409,385,425,533,312782574,312782628,417,312782546,  # Visual
            1011,1002,1027,1018                                         # Auditory
        ]

    return sensoryIds                      

# def loadAllDiffuse(searchPath, channelName):

#     # Select the files to load
#     fileList = [x for x in listdir(searchPath) if ('diffFluo_'+channelName) in x]

#     # List of tuples containing treatment and mouseID for each file
#     keysNames = [('HFD', x.split('_')[0]) if x.startswith('HF') else ('CTR', x.split('_')[0]) for x in fileList]

#     # Assemble all animals in a list of dataFrames
#     dfList = [pd.read_csv(join(searchPath,x), index_col='regionID') for x in fileList]

#     # Concatenate all the dataframes
#     dfMerged = pd.concat(dfList, axis=1,
#                             join='outer',
#                             names=['treat','mouse','params'],
#                             keys=keysNames,
#                             sort=True)
    
#     # Drop the row relative to the background (regionID=0)
#     dfMerged = dfMerged.drop(index=0)
                            
#     return dfMerged

# ------------------------------------------------------------------------------
# FILE LOADING FUNCTIONS
# ------------------------------------------------------------------------------

def loadRegionsDf(searchPath:str, channelName:str, mouseName:str, normCellIntens:bool = False, verbose:bool = True):
    '''
        Loads diffuse and dots data for one mouse and creates a dataframe with 
        raw information (rows = brain regions)

        Parameters
        ----------
        searchPath:str
            folder path for all the raw data
        channelName:str
            either 'wfa' or 'pv'
        mouseName:str
            ID of the mouse to load (e.g., 'CC1A')
        normCellIntens:str
            whether or not to normalize the intensity of cells
        verbose: bool
            additional information about discarded areas

        Returns
        -------
        joined:pd.DataFrame
            output dataframe
        '''

    # Select candidate files based on staining and mouse name
    fileList = listdir(searchPath)
    validFiles = [x for x in fileList if (channelName in x) and (mouseName in x)]

    # Load the CSV files in appropriate variables
    foundDiffuse = False
    foundDots = False
    for fileName in validFiles:
        if '_diffFluo_' in fileName:
            diffuse = pd.read_csv(
                join(searchPath, fileName),
                index_col='regionID')
            foundDiffuse = True
        elif '_dots_' in fileName:
            dots = pd.read_csv(
                join(searchPath, fileName),
                index_col='cellID')
            if normCellIntens:
                # dots['fluoMean'] = dots['fluoMean'] / dots['fluoMean'].mean()  # Old normalization per animal
                dots['fluoMean'] = dots['fluoMean'] / 255  # New normalization of absolute values (0-255) to (0-1)
            foundDots = True
    
    # If verbose, print out information
    if verbose:
        if not foundDiffuse:
            print(f"Diffuse Fluorescence file not found! [staining: {channelName}, mouse:{mouseName}]")
        if not foundDots:
            print(f"Dots file not found! [staining: {channelName}, mouse:{mouseName}]")

    # Group dots by area and aggregate cell number and total intensity
    grouped = dots.groupby(by='regionID')
    dotsGroup = grouped.agg(
            numCells = ('parentImg','count'),
            fluoCellSum = ('fluoMean','sum'),
    )

    # Merge the diffuse and the dots dataframes
    joined = diffuse.join(dotsGroup, how='left')
    # Fill areas without cells (NaN) with 0 numCells
    joined['numCells'] = joined['numCells'].fillna(0)

    # Remove the row with regionID = 0 (background)
    joined = joined.drop(index=0)

    return joined


def allMiceRegions(searchPath:str, channelName:str, normCellIntens:bool = False, verbose:bool=False):
    '''
        Load a unified dataFrame for brain regions data of all the mice in the
        searchPath folder 

        Parameters
        ----------
        searchPath:str
            folder path for all the raw data
        channelName:str
            either 'wfa' or 'pv'
        verbose: bool
            additional information about discarded areas

        Returns
        -------
        dfMerged:pd.DataFrame
            output dataframe with metrics from all mice
        '''
    
    fileList = listdir(searchPath)

    # Create a list of the unique names of all the mice
    miceList = [x.split('_')[0] for x in fileList]
    miceList = list(set(miceList))  # remove duplicate
    miceList.sort()                 # sort alphabetically

    # Assemble all animals in a list of dataFrames
    dfList = [loadRegionsDf(searchPath, channelName, x, verbose=verbose, normCellIntens=normCellIntens) for x in miceList]

    # # List of tuples containing treatment and mouseID for each file
    # keysNames = [('HFD', x) if x.startswith('HF') else ('CTR', x) for x in miceList]

    # Concatenate all the dataframes
    dfMerged = pd.concat(
        dfList,
        axis=1,
        join='outer',
        names=['mouse','params'],
        keys=miceList,
        sort=True
    )
    return dfMerged


def loadDotsDf(searchPath:str,  mouseName:str, verbose:bool=True):
    '''
        Loads processed dots data for a single animal (WFA and PV data after
        colocalization analysis) and creates a dataframe.

        Parameters
        ----------
        searchPath:str
            folder path for the processed dots file (PV and WFA dots after
            colocalization analysis)
        mouseName:str
            ID of the mouse to load (e.g., 'CC1A')
        verbose:bool
            additional information about discarded areas

        Returns
        -------
        dots_df:pd.DataFrame
            output dataframe
        '''

    # Select candidate files based on staining and mouse name
    fileList = listdir(searchPath)
    suffix =  mouseName + '_combinedDots.csv'
    filename = join(searchPath,suffix)

    if not isdir(searchPath):
        raise Exception(f'Wrong path specified![{searchPath} is not a folder]')
    if not isfile(filename) or suffix not in fileList:
        raise FileNotFoundError(f'The selected folder does not contain data referred to mouse :{mouseName}\n')
    
    dots_df = pd.read_csv(filename, index_col='cellID')

    # Normalize intensity values between 0 and 1
    dots_df['fluoMeanWfa'] = dots_df['fluoMeanWfa'] / 255 
    dots_df['fluoMeanPv'] = dots_df['fluoMeanPv'] / 255
    return dots_df


def allMiceDots(searchPath:str, verbose:bool=True):
    '''
        Load a unified dataFrame for colocalized dot data of all the mice in the
        searchPath folder 

        Parameters
        ----------
        searchPath:str
            folder path for all the raw data
        verbose: bool
            additional information about discarded areas

        Returns
        -------
        dfMerged:pd.DataFrame
            output dataframe with dots from all mice
        '''

    fileList = listdir(searchPath)

    # Create a list of the unique names of all the mice
    miceList = [x.split('_')[0] for x in fileList]
    miceList = list(set(miceList))  # remove duplicate
    miceList.sort()                 # sort alphabetically

    # Group all animals in a list of dataFrames
    dfList = [loadDotsDf(searchPath, x, verbose=verbose) for x in miceList]

    dfMerged = pd.concat(dfList, axis=0)
    cellIDs = dfMerged.index.tolist()

    cellList = []
    for c in cellIDs:
        components = c.split('_')
        mouse = components[0]
        cellList.append(tuple([mouse, c]))
        
    mi = pd.MultiIndex.from_tuples(cellList, names = ['mouse', 'cell'])
    dfMerged.index = mi
    return dfMerged


def loadConnectomeFromFile(connectomeFile, atlas):
    """
    Loads the connectome file downloaded from the Supplementary Table 3 in Oh et al., 2014.

    Connection strenght is pre-processed in that ipsi- and contralateral connections 
    are summed. Then a multiindex, including brain areas at the coarse resolution
    is added both to the index and to the columns
    
    Parameters
    ----------
    connectomeFile:str
        Path to the connectome .xlsx file
    atlas: object
        Instance of the atlas from AbaTool

    Returns
    -------
    connectome:pd.DataFrame
        output dataframe with connection strength data
    """
    
    # IPSI connections
    ipsi = pd.read_excel(connectomeFile, sheet_name='W_ipsi', index_col=0)
    # Change acronyms to region IDs
    ipsi.index = atlas.acronyms_to_ids(ipsi.index.tolist())
    ipsi.columns = atlas.acronyms_to_ids(ipsi.columns.tolist())

    # CONTRA connections
    contra = pd.read_excel(connectomeFile, sheet_name='W_contra', index_col=0)
    # Change acronyms to region IDs
    contra.index = atlas.acronyms_to_ids(contra.index.tolist())
    contra.columns = atlas.acronyms_to_ids(contra.columns.tolist())

    # Sum contra and ipsi connectivity to get the final connectome
    connectome = contra.add(ipsi, fill_value=0).copy()

    majdiv = atlas.get_major_divisions_ids()
    # Change index to multiindex, adding coarse level area IDs
    coarse_list, _  = atlas.match_structure_id_lists2(connectome.index.tolist(), majdiv)
    multiindex_list = [(coarse, fine) for coarse, fine in zip(coarse_list, connectome.index.tolist())]
    multiIdx = pd.MultiIndex.from_tuples(multiindex_list, names=['coarse','mid'])
    connectome.index = multiIdx
    # Change columns to multiindex, adding coarse level area IDs
    coarse_list, _  = atlas.match_structure_id_lists2(connectome.columns.tolist(), majdiv)
    multiindex_list = [(coarse, fine) for coarse, fine in zip(coarse_list, connectome.columns.tolist())]
    multiIdx = pd.MultiIndex.from_tuples(multiindex_list)
    connectome.columns = multiIdx

    return connectome

def loadGeneExpressionFromFile(aba_gene_expr_path:str,
                                metric:str = 'energy'):
                                
    if metric.lower() == 'energy':
        expr = pd.read_csv(join(aba_gene_expr_path,"gene_expression_ABA_energy.csv"), index_col=0)
    elif metric.lower() == 'density':
        expr = pd.read_csv(join(aba_gene_expr_path,"gene_expression_ABA_density.csv"), index_col=0)
    elif metric.lower() == 'intensity':
        expr = pd.read_csv(join(aba_gene_expr_path,"gene_expression_ABA_intensity.csv"), index_col=0)
    else:
        raise ValueError('Incorrect metric, only energy, intensity and density supported.')
    
    #convert strings to int
    expr.columns = pd.to_numeric(expr.columns)
    return expr

def loadGOResults(  webgestalt_folder:str,
                    redundancy_reduction = 'affinity_propagation',
                    ):
    folder_name = split(webgestalt_folder)[1]
    res = folder_name.split('_')[2]

    ap = "enriched_geneset_ap_clusters_wg_"+res+".txt"
    ap_file = join(webgestalt_folder,ap)

    topset = "enriched_geneset_wsc_topsets_wg_"+res+".txt"
    topset_file = join(webgestalt_folder,topset)

    data = "enrichment_results_wg_"+res+".txt"
    data_file = join(webgestalt_folder,data)

    go_df = pd.read_csv(data_file, sep='\t', index_col=0 )

    if redundancy_reduction == 'affinity_propagation':
        with open(ap_file, mode = 'r') as f:
            terms = f.read().splitlines()
            ap_terms = []
            for line in terms:
                ap_terms.append(line.split(sep = '\t')[0])
        go_df = go_df.loc[ap_terms]
    elif redundancy_reduction == 'weighted_set_cover':
        with open(topset_file, mode = 'r') as f:
            terms = f.read().splitlines()
            wsc_terms = []
            for line in terms[1:]:
                wsc_terms.append(line)
        go_df = go_df.loc[wsc_terms]
    elif redundancy_reduction == None:
        pass
    else:
        raise ValueError(f'{redundancy_reduction} is not supported as redundancy reduction method. Please choose between affinity_propagation, weighted_set_cover or None')

    go_df.rename(columns = {'description':'set_name',
                            'size':'set_size',
                            'pValue':'p_value',
                            'enrichmentRatio':'enrichment_ratio'},
                            inplace = True)

    go_df = go_df.drop(['link','overlapId'], axis=1)
    percentages = (go_df['overlap'].values/go_df['set_size'].values)*100
    go_df.insert(2,column='percentage_of_term', value =  percentages)
    
    go_df.sort_values(by = 'enrichment_ratio',ascending=False, inplace = True)
    return go_df

def loadMatrisomeTerms(matrisome_file:str) -> dict:
    matrisome = pd.read_excel(matrisome_file)
    matrisome = matrisome.loc[matrisome['Division'] != 'Retired']

    matrisome_terms = {y:matrisome.loc[matrisome['Category'] == y, 'Gene Symbol'].tolist() for y in matrisome['Category'].unique()}
    return matrisome_terms  


if __name__ == "__main__":

    resultsPath = 'D:\PizzorussoLAB\proj_PNN-highFatDiet\RESULTS\diffFluo'
    channelName = 'wfa'
