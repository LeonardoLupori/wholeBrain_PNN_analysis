import pandas as pd
import typing as tp


# ------------------------------------------------------------------------------
# FILE LOADING FUNCTIONS
# ------------------------------------------------------------------------------

def subsetDf_to_mid(rawDf, Atlas, normEnergy, normDiffuse,normalize:bool=True, verbose:bool=False, 
    not_to_drop:tp.Sequence[int] = []) -> pd.DataFrame:
    '''
    Converts the multi index dataframe generated by multiIndexDf_from_fineDf
    in a dataframe of regions at mid ontology levels

    Parameters
    ----------
    rawDf: pd.DataFrame
        Multiindex df generated by multiIndexDf_from_fineDf.
    normalize: bool=True
        Whether to normalize diffuseFluorescence and energy. If True, for each
        animal, energy and diffFluo are divided by energy and diffFluo calculated
        for the entire brain.
    not_to_drop: list = []
        List of area IDs to force the algorithm not to drop.
    verbose: bool=False
        Prints out additional information about areas that are dropped
    Returns
    -------
    df: pd.DataFrame
        multiindex dataframe
    '''
    
    grouped = rawDf.groupby(by=['coarse','mid'], axis=0).sum()

    # Calculate extract dataFrames of single measurements
    areaPx = grouped.xs('areaPx',axis='columns',level='params')
    diffFluo = grouped.xs('diffFluo',axis='columns',level='params')
    areaMm = grouped.xs('areaMm2',axis='columns',level='params')
    numCells = grouped.xs('numCells',axis='columns',level='params')
    fluoCellSum = grouped.xs('fluoCellSum',axis='columns',level='params')

    # Calculate measurements
    diffFluo = diffFluo.divide(areaPx)
    density = numCells.divide(areaMm)
    intensity = fluoCellSum.divide(numCells)
    energy = fluoCellSum.divide(areaMm)

    # Normalize energy and diffuseFluo if requested
    if normalize:
        # Normalize data
        energy = energy.divide(normEnergy)
        diffFluo = diffFluo.divide(normDiffuse)

    # Concatenate the dataframes in a single multi-animal, multi-measurement
    temp = pd.concat(
        [diffFluo,density,intensity,energy],
        axis=1,
        keys=['diffuseFluo','density','intensity','energy'],
        names=['params','mouse']
    )
    # Reorder levels of the columns
    temp = temp.reorder_levels(['mouse','params'], axis=1).sort_index(axis=1,level='mouse')
    
    # Sort rows based on the areas ontology
    temp['sort'] = Atlas.ids_to_graph_order(temp.index.get_level_values('mid').to_list())
    temp = temp.sort_values(by=('sort')).drop('sort', axis=1)

    # Drop non informative areas
    midIndex = temp.index.get_level_values('mid').tolist() 
    midlist = Atlas.get_midontology_structures_ids()
    midlist.extend(not_to_drop)
    areasToDrop = [x for x in midIndex if x not in midlist]
    temp = temp.drop(index=areasToDrop, level='mid')
    if verbose:
        for droppedArea in areasToDrop:
            print(f"- Area ID:{droppedArea} - Name: {Atlas.get_structures_by_id([droppedArea])[0]['name']} dropped.")

    return temp

def subsetDf_to_coarse( rawDf, Atlas, normEnergy, normDiffuse,normalize:bool=True, verbose:bool = False, 
        not_to_drop:tp.Sequence[int] = []) -> pd.DataFrame:
        '''
        Converts the multi index dataframe generated by multiIndexDf_from_fineDf
        in a dataframe of regions at coarse ontology levels

        Parameters
        ----------
        rawDf: pd.DataFrame
            Multiindex df generated by multiIndexDf_from_fineDf.
        normalize: bool=True
            Whether to normalize diffuseFluorescence and energy. If True, for each
            animal, energy and diffFluo are divided by energy and diffFluo calculated
            for the entire brain.
        verbose: bool=False
            Prints out additional information about areas that are dropped
        not_to_drop: list = []
            List of area IDs to force the algorithm not to drop.
        Returns
        -------
        df: pd.DataFrame
            multiindex dataframe
        '''
        
        grouped = rawDf.groupby(by='coarse', axis=0).sum()

        # Calculate extract dataFrames of single measurements
        areaPx = grouped.xs('areaPx',axis=1,level='params')
        diffFluo = grouped.xs('diffFluo',axis=1,level='params')
        areaMm = grouped.xs('areaMm2',axis=1,level='params')
        numCells = grouped.xs('numCells',axis=1,level='params')
        fluoCellSum = grouped.xs('fluoCellSum',axis=1,level='params')

        # Calculate measurements
        diffFluo = diffFluo.divide(areaPx)
        density = numCells.divide(areaMm)
        intensity = fluoCellSum.divide(numCells)
        energy = fluoCellSum.divide(areaMm)

        # Normalize energy and diffuseFluo if requested
        if normalize:

            # Normalize data
            energy = energy.divide(normEnergy)
            diffFluo = diffFluo.divide(normDiffuse)

        # Concatenate the dataframes in a single multi-animal, multi-measurement
        temp = pd.concat(
            [diffFluo,density,intensity,energy],
            axis=1,
            keys=['diffuseFluo','density','intensity','energy'],
            names=['params', 'mouse']
        )
        # Reorder levels of the columns
        temp = temp.reorder_levels(['mouse','params'], axis=1).sort_index(axis=1,level='mouse')
        
        # Sort rows based on the areas ontology
        temp['sort'] = Atlas.ids_to_graph_order(temp.index.get_level_values('coarse').to_list())
        temp = temp.sort_values(by=('sort')).drop('sort', axis=1)

        # Drop non informative areas
        midIndex = temp.index.get_level_values('coarse').tolist() 
        midlist = Atlas.get_major_divisions_ids()
        midlist.extend(not_to_drop)
        areasToDrop = [x for x in midIndex if x not in midlist]
        temp = temp.drop(index=areasToDrop)
        if verbose:
            for droppedArea in areasToDrop:
                print(f"- Area ID:{droppedArea} - Name: {Atlas.get_structures_by_id([droppedArea])[0]['name']} dropped.")

        return temp



