from allensdk.core.structure_tree import StructureTree
from allensdk.core.reference_space_cache import ReferenceSpaceApi, ReferenceSpaceCache
from allensdk.api.queries.ontologies_api import OntologiesApi
import typing as tp
import pandas as pd
import numpy as np
import re
import json
import os
from warnings import warn


class Atlas(StructureTree):
    '''Atlas class provides methods to access and process anatomical information of the Allen Brain Atlas, expanding functionalities of ABA StructureTree.
    Methods listed here are only those specific of this class, for inherited StructureTree methods go to class definition.

    Methods
    -------
    ids_to_graph_order(ids)
        converts list of ids to list of graph orders (i.e., identifiers of the level in the hierarchy)  
    ids_to_acronyms(ids)
        converts list of ids to list of acronyms
    ids_to_names(ids)
        converts list of ids to list of area names
    ids_to_colors(ids)
        converts list of ids to list of colors (different color modes supported)
    acronyms_to_ids(acronyms)
        converts list of colors to list of ids
    names_to_ids(names)
        converts list of colors to list of ids
    remap_area_ids(ids, newmap)
        substitutes ids within a list based on dictionary
    match_structure_id_lists(list_to_match, reference_list)
        convert brain area list so that it is at the same ontology level of a reference list
    get_major_divisions_ids()
        returns the list of ids of  12 major divisions 
    get_midontology_structures_ids()
        returns the list of ids of  316 mid-ontology areas 
    get_cortical_structures_ids()
        returns the list of ids of 43 cortical areas
    get_layer_from_area_id(ids)
        form a list of area ids returns the layer number(1,2/3,4,5,6)
    structure_descends_from_any(area_id, list_of_parents)
        check if structure descends from any of a list of putative parent ids
    '''

    # obsolete_parent_areas ={
    #         934: {
    #             'old_acronym' :'ENTmv',
    #             'new_id' :926,
    #             'new_acronym' :'ENTm',
    #             'new_acronym_list' :['ENTm'],
    #             'parent_id':909,
    #             'parent_acronym':'ENT'},
    #         22: {
    #             'old_acronym' :'PTLp',
    #             'new_id':22,
    #             'new_acronym' :'VISa-VISrl',
    #             'new_acronym_list' :['VISa', 'VISrl'],
    #             'parent_id' : 669,
    #             'parent_acronym':'VIS'}
    #             }
    # obsolete_areas = {
    #         560 : {
    #             'old_acronym': 'CNspg',
    #             'new_id':607,
    #             'new_acronym' :'DCO-VCO',
    #             'new_acronym_list' :['DCO', 'VCO'],
    #             'parent_id' : 607,
    #             'new_parent_acronym':'CN'},
    #         112 : {
    #             'acronym': 'CNlam',
    #             'new_id':607,
    #             'new_acronym' :'DCO-VCO',
    #             'new_acronym_list' :['DCO', 'VCO'],
    #             'parent_id' : 607,
    #             'new_parent_acronym':'CN'}
    #         }
    # CONSTRUCTOR
    def __init__(self, nodes=None, resolution=25, reference_space_key='annotation/ccf_2017'):
        '''Anatomical atlas for navigating ABA structure convention.
        
        Parameters
        ----------
        nodes : list of dict
            Each specifies a structure. Fields are:
            
            'acronym' : str
                Abbreviated name for the structure.
            'rgb_triplet' : str
                Canonical RGB uint8 color assigned to this structure
            'graph_id' : int
                Specifies the structure graph containing this structure.
            'graph_order' : int
                Canonical position in the flattened structure graph.
            'id': int
                Unique structure specifier.
            'name' : str
                Full name of structure.
            'structure_id_path' : list of int
                This structure's ancestors (inclusive) from the root of the 
                tree.
            'structure_set_ids' : list of int
                Unique identifiers of structure sets to which this structure 
                belongs. 

        resolution: int, dimension of atlas voxels (um)
        reference_space_key: str, identifying anatomy reference space
        '''
        
        
        
        if nodes == None:

            nodes = OntologiesApi().get_structures_with_sets(
                        strategy='lazy',
                        pre=self.clean_structures,
                        post=lambda x: self.clean_structures(x), 
                        path=ReferenceSpaceCache(resolution, reference_space_key, manifest='manifest.json').get_cache_path(None, 'STRUCTURE_TREE'),
                        structure_graph_ids=1, 
                        **ReferenceSpaceCache.cache_json())
        
        else:
            if os.path.isfile(nodes):
                with open(nodes) as f:
                    nodes = json.load(f)


        super().__init__(nodes)


    # Converts ids into graph_orders
    def ids_to_graph_order(self, ids:tp.Sequence[int]):
        '''
        converts ids into graph_orders

        Parameters
        ----------
        ids : list
            list of ABA ids

        Returns
        -------
        graph_order: list
            list graph of order numebers        
        '''

        to_fn = lambda x: x['graph_order']
        graph_order = self.nodes_by_property('id', ids, to_fn)
        return graph_order


    # Converts list of acronyms into ids
    def acronyms_to_ids(self, acronyms:tp.Sequence[str]):
        '''
        converts ids into acronyms

        Parameters
        ----------
        acronyms : list
            list of ABA ids

        Returns
        -------
        ids: list
            list of acronyms        
        '''

        to_fn = lambda x: x['id']
        ids = self.nodes_by_property('acronym', acronyms, to_fn)
        return ids


    # Converts list of ids into acronyms
    def ids_to_acronyms(self, ids:tp.Sequence[int]):
        '''
        converts acronyms into ids

        Parameters
        ----------
        ids : list
            list of acronyms

        Returns
        -------
        acronyms: list
            list of ABA ids        
        '''

        to_fn = lambda x: x['acronym']
        ids = self.nodes_by_property('id', ids, to_fn)
        return ids


    # Get colors from ids (rgb triplet or hex)
    def ids_to_colors(self, ids:tp.Sequence[int], color_model = "rgb") -> tp.Sequence[list]:
        '''
        converts ids into graph_orders

        Parameters
        ----------
        ids : list
            list of ABA ids
        color_model: str
            str in ['rgb', 'hex', 'rgb_norm', 'rgb_plotly]

        Returns
        -------
        colors: list
            list of colors        
        '''

        if color_model == "hex":
            to_fn = lambda y: f"{y['rgb_triplet'][0]:02x}{y['rgb_triplet'][1]:02x}{y['rgb_triplet'][2]:02x}" 
        elif color_model == "rgb":
            to_fn = lambda x: x['rgb_triplet']
        elif color_model == "rgb_norm":
            to_fn = lambda x: [c/255 for c in x['rgb_triplet']]
        elif color_model == "rgb_plotly":
            to_fn = lambda x: f"rgb({x['rgb_triplet'][0]},{x['rgb_triplet'][1]},{x['rgb_triplet'][2]})"
        colors = self.nodes_by_property('id', ids, to_fn)
        return colors


    # Get name of areas based on ids
    def ids_to_names(self, ids:tp.Sequence[int]) -> tp.Sequence[str]:
        '''
        converts ids into names

        Parameters
        ----------
        ids : list
            list of ABA ids

        Returns
        -------
        names: list
            list ABA area names        
        '''

        to_fn = lambda x: x['name']
        names = self.nodes_by_property('id', ids, to_fn)
        return names
    

    # Get ids of areas based on names
    def names_to_ids(self, names:tp.Sequence[str]) -> tp.Sequence[int]:
        '''
        converts names into ids

        Parameters
        ----------
        names : list
            list of ABA area names

        Returns
        -------
        ids: list
            list of ABA ids      
        '''

        to_fn = lambda x: x['id']
        ids = self.nodes_by_property('name', names, to_fn)
        return ids


    # Substitute area ids in pre-defined list
    def remap_area_ids(self, ids: tp.Sequence[int], newmap:tp.Mapping) -> tp.Sequence[int]:
        '''
        substitutes area ids in a given list

        Parameters
        ----------
        ids : list
            list of ABA ids
        newmap : dict
            areas to be substituted{old_id: new_id}

        Returns
        -------
        new_ids: list
            list graph order numebers        
        '''

        new_ids = [newmap[i]  if i in newmap.keys() else i for i in ids ]
        return new_ids
    

    #DEPRECATED
    # Align area id lists to the resolution level of another area id lists (supports conversion from a couple of obsolete ids [TO BE CHECKED AGAIN]])
    def match_structure_id_lists(self, list_to_match:tp.Sequence[int], reference_list:tp.Sequence[int], verbose:bool = False) -> tp.List[int] :
        '''
        allign area id lists to the resolution level of another area id lists

        Parameters
        ----------
        list_to_match : list
            list of ABA ids
        list_to_match : list
            list of ABA ids

        Returns
        -------
        matched_list: list
            list of ABA ids, same length of input list
        non_matchable_ids: list
            list of ABA ids for which the ontology matching failed (reported as None in the matched list)
        '''
        warn('the function is deprecated, use match_structure_id_lists2, instead',DeprecationWarning, stacklevel=2)
        matched_list = []
        non_matchable_ids = []
        for area in list_to_match:
            matched  = False
            for ref_area in reference_list:
                if area == ref_area:
                    matched_list.append(area)
                    matched  = True
                    break
                elif self.structure_descends_from(area, ref_area) :
                    matched_list.append(ref_area)
                    matched  = True
                    break
                elif self.structure_descends_from(ref_area, area) :
                    matched_list.append(area)
                    matched  = True
                    break
                elif self.structure_descends_from(area, 73) :
                    matched_list.append(73)
                    matched  = True
                    if verbose:                
                        print('found structures of the ventricular system in the proposed id list')
                    break
                # elif any([self.structure_descends_from(area,obsPar) for obsPar in self.obsolete_parent_areas.keys()]):
                #     # print('ok1')
                #     for obsPar in self.obsolete_parent_areas.keys():
                #         # print('ok2')
                #         if self.structure_descends_from(area,obsPar):
                #             matched_list.append(self.obsolete_parent_areas[obsPar]['new_id'])
                #             matched  = True
                #             # print('ok')
                #             if verbose:    
                #                 old_acr = self.obsolete_parent_areas[obsPar]['old_acronym']
                #                 print(f'Obsolete areas in the source list { old_acr }')
                #                 print(area, obsPar)
                #             break
                #     break
                # elif any([area == obsPar for obsPar in self.obsolete_areas.keys()]):
                #     for obsPar in self.obsolete_areas.keys():
                #         if area == obsPar :
                #             matched_list.append(self.obsolete_areas[obsPar]['new_id'])
                #             matched  = True
                #             if verbose:
                #                 old_acr = self.obsolete_areas[obsPar]['old_acronym']    
                #                 print(f'Obsolete areas in the source list { old_acr }')
                #             break
                #     break
            if matched == False:
                matched_list.append(None)
                non_matchable_ids.append(area)            
        return matched_list, non_matchable_ids


    # Align area id lists to the resolution level of another area id lists (supports conversion from a couple of obsolete ids [TO BE CHECKED AGAIN]])
    def match_structure_id_lists2(self, list_to_match:tp.Sequence[int], reference_list:tp.Sequence[int], verbose:bool = False) -> tp.List[int] :
        '''
        allign area id lists to the resolution level of another area id lists

        Parameters
        ----------
        list_to_match : list
            list of ABA ids
        reference_list : list
            list of ABA ids

        Returns
        -------
        matched_list: list
            list of ABA ids, same length of input list
        non_matchable_ids: list
            list of ABA ids for which the ontology matching failed (reported as None in the matched list)
        '''

        matched_list = []
        non_matchable_ids = []
        for area in list_to_match:
            if  area in reference_list:
                matched_list.append(area)
            elif self.parent_ids([area])[0] in reference_list:
                matched_list.append(self.parent_ids([area])[0])
            elif area in self.parent_ids(reference_list):
                matched_list.append(area)
            elif self.structure_descends_from(area, 73):
                matched_list.append(73)
                if verbose:
                    print('Area of the ventricual system found')
            elif self.structure_descends_from_any(area, reference_list):
                for ancestor in self.ancestor_ids([area])[0]:
                    if ancestor in reference_list:
                        matched_list.append(ancestor)
                        if verbose:
                            print(area)
                            print(f"{self.ids_to_names([area])} not present in reference_list but has ancestor {ancestor} in that list.")
                        break
            elif self.parent_ids([area])[0] in self.parent_ids(reference_list):
                matched_list.append(self.parent_ids([area])[0])
                if verbose:
                    print(f"{self.ids_to_names([area])} not present in reference_list but shares parents with areas in that list.")            
            elif any(self.structures_have_common_ancestors(area, ref_area) for ref_area in reference_list):
                counter = 1
                keep = True
                
                while keep:
                    for ref_area in reference_list:
                        ancestors = self.ancestor_ids([ref_area])[0]
                        # print(f'{counter} and {len(ancestors)}')
                        if (counter<=  len(ancestors)-1)and(ancestors[counter] in self.ancestor_ids([area])[0]):
                            if verbose:
                                print(f'found common ancestor between area {self.ids_to_names([area])} and {self.ids_to_names([ref_area])}')
                            
                            matched_list.append(ancestors[counter])

                            keep = False
                            break
                        if ref_area == reference_list[-1]:
                            matched_list.append(None)
                            non_matchable_ids.append(area)
                            keep= False
                            break
                        
                    counter += 1
            else:
                matched_list.append(None)
                non_matchable_ids.append(area)
                if verbose:
                    print(f"Unable to match area {self.ids_to_names([area])} with reference list: dropped")
        return matched_list, non_matchable_ids


    # Get list of ids of 12 major anatomical structures
    def get_major_divisions_ids(self):
        '''
        Returns list of ids of 12 major anatomical structures

        Returns
        -------
        major_ids: list
            list of ABA ids
        '''
        overlap = lambda x: (set([687527670]) & set(x['structure_set_ids']))
        filtered_dicts =  self.filter_nodes(overlap)
        return [x['id'] for x in filtered_dicts]


    # Get list of ids of mid-ontology areas
    def get_midontology_structures_ids(self):
        '''
        Returns list of ids of 316 mid-ontology structures

        Returns
        -------
        mid_ids: list
            list of ABA ids
        '''
        overlap = lambda x: (set([167587189]) & set(x['structure_set_ids']))
        filtered_dicts =  self.filter_nodes(overlap)
        return [x['id'] for x in filtered_dicts]


    # Get list of ids of cortical areas (43)
    def get_cortical_structures_ids(self):
        '''
        Returns list of ids of 43 cortical structures

        Returns
        -------
        ctrx_ids: list
            list of ABA ids
        '''
        
        overlap = lambda x: (set([688152357]) & set(x['structure_set_ids']))
        filtered_dicts =  self.filter_nodes(overlap)
        ctrx_ids = [x['id'] for x in filtered_dicts]
        return ctrx_ids


    # Return string correspondent to the number of cortical layer (None if area id is not descendant of Isocortex)
    def get_layer_from_area_id(self, ids):
        '''
        Returns list of layers from ids of cortical structures

        Parameters
        ----------
        ids : list
            list of ABA ids

        Returns
        -------
        layer_ids: list
            list of ABA ids (None if structure is not isocortex descendant)
        '''

        isocrtx_ids = self.get_cortical_structures_ids()
        area_names = self.ids_to_names(ids)
        func = lambda x: self.extract_layer_number(x)
        is_isocortex = lambda x: self.structure_descends_from_any(x,isocrtx_ids)
        layer_ids  = list(map(lambda x: func(x) if is_isocortex else None, area_names))
        return layer_ids 


    # Get layer number based on area name
    def extract_layer_number(self,string):
        '''
        Extract layer number from name string (re module)

        Parameters
        ----------
        string : str
            string of area name

        Returns
        -------
        layer: str
            correspondent cortical layer

        '''

        pattern_str = re.compile(r'[Layer]*\w*\s*(\d)\w*\s*(/*)\s*(\d*)\w*', re.IGNORECASE)
        match_obj = re.search(pattern_str , repr(string) )

        layer_tuple = match_obj.groups()


        layer = "".join(list(filter(None, layer_tuple)))
        
        return layer
    

    def find_first_common_ancestor(self, area1, area2):
        '''Finds the first ancestor in area1 ancestor list that is present also in area2 ancestor list.
        
        Parameters
        ----------
        area1 : int
            Id of the putative child structure.
        area2 : int
            Id of the putative parent structure.
            
        Returns
        -------
        acestor 1 : int
            fist common ancestor.
        '''
        for ancestor1 in self.ancestor_ids([area1]):
            if ancestor1 in self.ancestor_ids([area2]):
                return ancestor1
                

    def structures_have_common_ancestors(self, area1, area2):
        '''Tests whether one structure has common ancestors with another.
        
        Parameters
        ----------
        area1 : int
            Id of the first structure.
        area2 : int
            Id of the second structure.
            
        Returns
        -------
        bool : bool
            True if the structures have common ancest. Otherwise False.
        '''

        inters = set(self.ancestor_ids([area1])[0]).intersection(set(self.ancestor_ids([area2])[0]))
        return len(inters) > 0


    def structure_is_child_of(self, child_id, parent_id):
        '''Tests whether one structure descends from another.
        Note: The function returns True only if the child_id structure lies just one step lower in the hierarchy with respect to the parent_id.
        
        Parameters
        ----------
        child_id : int
            Id of the putative child structure.
        parent_id : int
            Id of the putative parent structure.
            
        Returns
        -------
        bool :
            True if the structure specified by child_id is a child of 
            the one specified by parent_id. Otherwise False.
        
        '''
    
        return parent_id in self.parent_ids([child_id])[0]


    def structure_descends_from_any(self, area_id, list_of_parents):
        '''
        Checks if area descends from any area in a given list of candidate ids

        Parameters
        ----------
        area_id: int
            scalar, area id to check
        list_of_parents: list
            parents to be screened

        Returns
        -------
        logical: bool
            True if parent is found in list_of_parents, False otherwise
        '''

        ancestors = self.ancestor_ids([area_id])[0]
        inters = set(ancestors).intersection(set(list_of_parents))
        logical = len(inters)>0
        return logical
    

    def structure_is_child_of_any(self, area_id, list_of_parents):
        '''
        Checks if area is child of any area in a given list of candidate ids. 
        Note: stict descendance is considered. The function returns true only if areas are just one step apart in the hierarchy

        Parameters
        ----------
        area_id: int
            scalar, area id to check
        list_of_parents: list
            parents to be screened

        Returns
        -------
        logical: bool
            True if parent is found in list_of_parents, False otherwise
        '''

        area_dict = self.nodes(area_id)
        inters = set(area_dict[0]['structure_id_path'][-2]).intersection(set(list_of_parents))
        logical = len(inters)>0
        return logical

    def get_mid_ontology_anatomy_sets(self): 
        ''' Returns a dict mapping each major division id to the correspondent descendant in the midontology area id set
        
        Returns
        -------
        anatomy_sets:dict
            keys: 12 major divisions
            values: lists of mid-ontology area ids descending from that major division
        
        '''
        macro_areas = self.get_major_divisions_ids()
        mid_areas = self.get_midontology_structures_ids()
        anatomy_sets = {}
        
        for majdiv in macro_areas:
            majdiv_name = self.ids_to_names([majdiv])[0]
            anatomy_sets[majdiv_name] = list(set(self.descendant_ids([majdiv])[0]).intersection(set(mid_areas)))
        return anatomy_sets

class AnatomyDataFrameManager():
    '''
    Class for managing anatomy-related multiindex dataframes

    Methods
    -------
    multiIndexDf_from_fineDf(fineDf, verbose=False)
        builds up the multiindex dataframe starting from leaf structure nodes
    load_structures_as_df(structuresJsonPath)
        loads structures from json file in the form of dataframe

    DEPRECATED
    multiIndex_to_mid_df(multiindex_df, normalize = True,  verbose = False)
        collapses multiindex dataframe to mid-ontology dataframe
    multiIndex_to_coarse_df(multiindex_df, normalize = True,  verbose = False)
        collapses multiindex dataframe to coarse dataframe
    multiIndex_to_corticalLayer_df(multiIndex_df, normalize = True)
        extracts cortical layers data and organizes them in a dataframe
    '''

    def __init__(self, atlas: Atlas):
        '''
        Dataframe Manager

        Parameters
        ----------
        atlas: AbaTool.Atlas object
            provides anatomical information
        '''
        self._atlas = atlas
    
    def multiIndexDf_from_fineDf(self, fineDf, verbose=False):
        '''
        Creates the multiindex df from fine-ontology df

        Parameters
        ----------
        fineDf: pd.DataFrame
            df containing fluorescence values. Indeces are represented by leaf-node brain areas
        verbose: bool
            additional information about discarded areas
        Returns
        -------
        df: pd.DataFrame
            multiindex dataframe
        '''

        id_list = fineDf.index.tolist()

        midOntologyIDs = self._atlas.get_midontology_structures_ids()
        # Match the fine IDs to the mid IDs
        mid, toDrop_mid = self._atlas.match_structure_id_lists2(id_list, midOntologyIDs)
        if verbose:
            print(f"While matching MID ONTOLOGY structures, {len(toDrop_mid)} structures were dropped:")
            for ID in toDrop_mid:
                print(f"\t- ID: {ID} - Name: {self._atlas.get_structures_by_id([ID])[0]['name']}")
        

        # Adjust dataFrame for COARSE ontology structures
        coarseOntologyIDs = self._atlas.get_major_divisions_ids()
        coarseOntologyIDs.append(1009)     # Add Fiber tracts to the list
        # Match the fine IDs to the coarse IDs
        coarse, toDrop_coarse = self._atlas.match_structure_id_lists2(id_list, coarseOntologyIDs)
        if verbose:
            print(f"While matching COARSE ONTOLOGY structures, {len(toDrop_coarse)} structures were dropped:")
            for ID in toDrop_coarse:
                print(f"\t- ID: {ID} - Name: {self._atlas.get_structures_by_id([ID])[0]['name']}")


        newIndices = [(a,b,c) for a,b,c in zip(coarse,mid,id_list)]
        
        # Mark as regions to drop those that have at least a None in one of the 3 idexes
        toRemove = [not all(x) for x in newIndices]
        
        # Copy of the original dataframe to use as output
        df = fineDf.copy()
        # Put the multiindex
        df.index = pd.MultiIndex.from_tuples(newIndices, names=["coarse","mid","fine"])
        # Remove regions with a None
        if any(toRemove):
            df.drop(df[toRemove].index, inplace=True)
        
        return df
    
    def regionsDf_to_fine(self, rawDf, normalize:bool=True, verbose:bool=False, 
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
        
        # Calculate extract dataFrames of single measurements
        areaPx = rawDf.xs('areaPx',axis=1,level='params')
        diffFluo = rawDf.xs('diffFluo',axis=1,level='params')
        areaMm = rawDf.xs('areaMm2',axis=1,level='params')
        numCells = rawDf.xs('numCells',axis=1,level='params')
        fluoCellSum = rawDf.xs('fluoCellSum',axis=1,level='params')

        # Calculate measurements
        diffFluo = diffFluo.divide(areaPx)
        density = numCells.divide(areaMm)
        intensity = fluoCellSum.divide(numCells)
        energy = fluoCellSum.divide(areaMm)

        # Normalize energy and diffuseFluo if requested
        if normalize:
            # Calculate the energy and diffFluo of the entire brain
            brainEnergy , brainDiffFluo = self.wholeBrainMetrics(rawDf)
            # Normalize data
            energy = energy.divide(brainEnergy)
            diffFluo = diffFluo.divide(brainDiffFluo)

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
        temp['sort'] = self._atlas.ids_to_graph_order(temp.index.get_level_values('fine').to_list())
        temp = temp.sort_values(by=('sort')).drop('sort', axis=1)

        # Drop non informative areas
        midIndex = temp.index.get_level_values('mid').tolist() 
        midlist = self._atlas.get_midontology_structures_ids()
        midlist.extend(not_to_drop)
        areasToDrop = [x for x in midIndex if x not in midlist]
        temp = temp.drop(index=areasToDrop, level='mid')
        if verbose:
            for droppedArea in areasToDrop:
                print(f"- Area ID:{droppedArea} - Name: {self._atlas.get_structures_by_id([droppedArea])[0]['name']} dropped.")

        return temp
    
    def regionsDf_to_mid(self, rawDf, normalize:bool=True, verbose:bool=False, 
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
            # Calculate metric for the total brain by using the raw Df
            brainAreaMm = rawDf.xs('areaMm2', axis=1, level='params').sum()
            brainAreaPx = rawDf.xs('areaPx', axis=1, level='params').sum()
            brainCellIntensity = rawDf.xs('fluoCellSum', axis=1, level='params').sum()
            brainDiffIntensity = rawDf.xs('diffFluo', axis=1, level='params').sum()
            # Calculate the energy and diffFluo of the entire brain
            brainEnergy = brainCellIntensity.divide(brainAreaMm)
            brainDiffFluo = brainDiffIntensity.divide(brainAreaPx)
            # Normalize data
            energy = energy.divide(brainEnergy)
            diffFluo = diffFluo.divide(brainDiffFluo)

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
        temp['sort'] = self._atlas.ids_to_graph_order(temp.index.get_level_values('mid').to_list())
        temp = temp.sort_values(by=('sort')).drop('sort', axis=1)

        # Drop non informative areas
        midIndex = temp.index.get_level_values('mid').tolist() 
        midlist = self._atlas.get_midontology_structures_ids()
        midlist.extend(not_to_drop)
        areasToDrop = [x for x in midIndex if x not in midlist]
        temp = temp.drop(index=areasToDrop, level='mid')
        if verbose:
            for droppedArea in areasToDrop:
                print(f"- Area ID:{droppedArea} - Name: {self._atlas.get_structures_by_id([droppedArea])[0]['name']} dropped.")

        return temp
        
    def regionsDf_to_coarse(self, rawDf, normalize:bool = True, verbose:bool = False, 
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
            # Calculate metric for the total brain by using the raw Df
            brainAreaMm = rawDf.xs('areaMm2', axis=1, level='params').sum()
            brainAreaPx = rawDf.xs('areaPx', axis=1, level='params').sum()
            brainCellIntensity = rawDf.xs('fluoCellSum', axis=1, level='params').sum()
            brainDiffIntensity = rawDf.xs('diffFluo', axis=1, level='params').sum()
            # Calculate the energy and diffFluo of the entire brain
            brainEnergy = brainCellIntensity.divide(brainAreaMm)
            brainDiffFluo = brainDiffIntensity.divide(brainAreaPx)
            # Normalize data
            energy = energy.divide(brainEnergy)
            diffFluo = diffFluo.divide(brainDiffFluo)

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
        temp['sort'] = self._atlas.ids_to_graph_order(temp.index.get_level_values('coarse').to_list())
        temp = temp.sort_values(by=('sort')).drop('sort', axis=1)

        # Drop non informative areas
        midIndex = temp.index.get_level_values('coarse').tolist() 
        midlist = self._atlas.get_major_divisions_ids()
        midlist.extend(not_to_drop)
        areasToDrop = [x for x in midIndex if x not in midlist]
        temp = temp.drop(index=areasToDrop)
        if verbose:
            for droppedArea in areasToDrop:
                print(f"- Area ID:{droppedArea} - Name: {self._atlas.get_structures_by_id([droppedArea])[0]['name']} dropped.")

        return temp

    def regionsDf_to_corticalLayers(self, rawDf, normalize:bool=True, verbose:bool=False) -> pd.DataFrame:
        '''
        Converts the multi index dataframe generated by multiIndexDf_from_fineDf
        in a dataframe of only cortical regions divided by layers

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
        Returns
        -------
        merged_df: pd.DataFrame
            multiindex dataframe
        '''

        latEntId = 918
        medEntId = 926
        isocortexId = 315

        # Create dataframes for cortex, lateral entorhinal and medial entorhinal and
        # merge all the three
        latEnt_df = rawDf.xs(latEntId, level='mid', drop_level=False).droplevel('coarse')
        medEnt_df = rawDf.xs(medEntId, level='mid', drop_level=False).droplevel('coarse')
        cortex_df = rawDf.loc[isocortexId]
        cortex_df = pd.concat([cortex_df, latEnt_df, medEnt_df], axis=0)

        # Add layer as an index level 
        layer = self._atlas.get_layer_from_area_id(cortex_df.index.get_level_values('fine'))
        layer = list(map(lambda x : '2/3' if x == '2' or x== '3' else x, layer))
        cortex_df['layer'] = layer 
        cortex_df = cortex_df.set_index('layer', append=True, drop=True)

        # Sum data from same layers (e.g., both layers 2 and 3 become 2/3; 
        # or layers 6a and 6b become layer 6)
        cortex_df = cortex_df.groupby(by = ['mid', 'layer']).sum()

        # Aggregate prerequisite measurements
        totalIntensity = cortex_df.xs('diffFluo', axis=1, level='params')
        areaPx = cortex_df.xs('areaPx', axis=1, level='params')
        numPnn = cortex_df.xs('numCells', axis=1, level='params')
        areaMm = cortex_df.xs('areaMm2', axis=1, level='params')
        fluoCellSum = cortex_df.xs('fluoCellSum', axis=1, level='params')
        # Compute metrics
        diffuseFluo = totalIntensity.divide(areaPx)     # Diffuse Fluorescence
        density = numPnn.divide(areaMm)                 # Cell Density
        intensity = fluoCellSum.divide(numPnn)          # Average cell intensity
        energy = fluoCellSum.divide(areaMm)             # Cell energy

        # Normalize energy and diffuse fluorescence on the entire brain for each 
        # animal individually
        if normalize:
            # Calculate metric for the total brain by using the raw Df
            brainAreaMm = rawDf.xs('areaMm2', axis=1, level='params').sum()
            brainAreaPx = rawDf.xs('areaPx', axis=1, level='params').sum()
            brainCellIntensity = rawDf.xs('fluoCellSum', axis=1, level='params').sum()
            brainDiffIntensity = rawDf.xs('diffFluo', axis=1, level='params').sum()
            # Calculate the energy and diffFluo of the entire brain
            brainEnergy = brainCellIntensity.divide(brainAreaMm)
            brainDiffFluo = brainDiffIntensity.divide(brainAreaPx)
            # Normalize data
            energy = energy.divide(brainEnergy)
            diffuseFluo = diffuseFluo.divide(brainDiffFluo)

        # Merge the various metrics in a single dataframe
        merged_df = pd.concat(
                    [diffuseFluo,density,intensity,energy],
                    axis=1,
                    keys=['diffuseFluo','density','intensity','energy'],
                    names=['params','mouse']
                )
        merged_df = merged_df.reorder_levels(['mouse','params'], axis=1).sort_index(axis=1,level='mouse')

        # Sort rows so that cortical layers are in ascending order
        merged_df = merged_df.sort_values(by=['mid','layer'], ascending=[True, True])
        # Sort cortical areas in a custom order that makes sense functionally
        if 22 in merged_df.index.get_level_values('mid'):
            customIndex = self.customCortexOrder(oldAtlasNumbers=True)
        else:
            customIndex = self.customCortexOrder()
        merged_df = merged_df.reindex(customIndex, level='mid')

        return merged_df

    def regionsDf_to_sensorySystem(self, rawDf, sensorySystem:str, normalize:bool=True, 
    mergeLayers:bool=True, verbose:bool=False) -> pd.DataFrame:
        """
            Areas
            - Visual:
            - Auditory:
            - Somatosensory:
        """

        assert sensorySystem in ['visual', 'auditory', 'somatosensory'],\
            'sensorySystem can only be one of: "visual", "auditory" or "somatosensory"'


        # Define the IDs of areas to keep as part of the selected sensory system
        allIds = rawDf.index.get_level_values('fine')
        if sensorySystem == 'visual':
            # This solves the bug of not including VISrl and VISa in visual areas
            targetIds = [a for a in allIds if 
                self._atlas.structure_descends_from(a, 669) or  # 669: Visual Areas
                self._atlas.structure_descends_from(a, 22)]     # 22: Posterior parietal association areas  
        elif sensorySystem == 'auditory':
            # 247: Auditory Areas
            targetIds = [a for a in allIds if self._atlas.structure_descends_from(a, 247)]
        elif sensorySystem == 'somatosensory':
            # 453: Somatosensory Areas
            targetIds = [a for a in allIds if self._atlas.structure_descends_from(a, 453)]

        # Filter all areas in rawDf and keep only the ones is the selected sensory system
        boolIdx = [x in targetIds for x in rawDf.index.get_level_values('fine')]
        sensory_df = rawDf.iloc[boolIdx,:]

        # Split areas in primary vs associative
        areatuple = []
        for a in sensory_df.index.get_level_values('fine'):
            layer = self._atlas.get_layer_from_area_id([a])[0]          # Get what layer this area corresponds to
            isPrimaryVis = self._atlas.structure_descends_from(a,385)   # Is a part of primary visual ctx?
            isPrimaryAud = self._atlas.structure_descends_from(a,1002)  # Is a part of primary auditory ctx?
            isPrimarySS = self._atlas.structure_descends_from(a,322)    # Is a part of primary somatosensory ctx?
            isPrimary = isPrimaryVis or isPrimaryAud or isPrimarySS
            if verbose:
                print(f"Area: {self._atlas.ids_to_acronyms([a])[0]}: ", end="")
            if isPrimary:
                areatuple.append(('primary', a, layer))
                if verbose:
                    print("Primary")
            else:
                areatuple.append(('associative', a, layer))
                if verbose:
                    print("Associative")

        multiIdx = pd.MultiIndex.from_tuples(areatuple, names=('sensory', 'fineid', 'layer'))
        sensory_df.index = multiIdx

        # Prepare output dataframe.
        if mergeLayers:
            sensory_df = sensory_df.groupby('sensory').sum()
        else:
            sensory_df = sensory_df.groupby(['sensory','layer']).sum()
        # Calculate extract dataFrames of single measurements
        areaPx = sensory_df.xs('areaPx',axis=1,level='params')
        diffFluo = sensory_df.xs('diffFluo',axis=1,level='params')
        areaMm = sensory_df.xs('areaMm2',axis=1,level='params')
        numCells = sensory_df.xs('numCells',axis=1,level='params')
        fluoCellSum = sensory_df.xs('fluoCellSum',axis=1,level='params')
        # Calculate metrics
        diffFluo = diffFluo.divide(areaPx)
        density = numCells.divide(areaMm)
        intensity = fluoCellSum.divide(numCells)
        energy = fluoCellSum.divide(areaMm)
        # Normalize metrics on the whole brain if requested
        if normalize:
            # Diff Fluo and energy for the whole brain
            brainEnergy, brainDiffFluo = self.wholeBrainMetrics(rawDf)
            # Normalize
            diffFluo = diffFluo.divide(brainDiffFluo)
            energy = energy.divide(brainEnergy)

        # Concatenate all metrics in a single output DataFrame
        output_df = pd.concat(
            [diffFluo,density,intensity,energy],
            axis=1,
            keys=['diffuseFluo','density','intensity','energy'],
            names=['params', 'mouse']
        )
        # Reorder levels of the columns
        output_df = output_df.reorder_levels(['mouse','params'], axis=1).sort_index(axis=1,level='mouse')

        return output_df

    def colocDf_to_sensorySystems(self, rawDf, sensorySystem:str, mergeLayers:bool=True, 
        verbose:bool=False) -> pd.DataFrame:

        # Define the IDs of areas to keep as part of the selected sensory system
        allIds = rawDf.index.get_level_values('fine')
        if sensorySystem == 'visual':
            # This solves the bug of not including VISrl and VISa in visual areas
            targetIds = [a for a in allIds if 
                self._atlas.structure_descends_from(a, 669) or  # 669: Visual Areas
                self._atlas.structure_descends_from(a, 22)]     # 22: Posterior parietal association areas  
        elif sensorySystem == 'auditory':
            # 247: Auditory Areas
            targetIds = [a for a in allIds if self._atlas.structure_descends_from(a, 247)]
        elif sensorySystem == 'somatosensory':
            # 453: Somatosensory Areas
            targetIds = [a for a in allIds if self._atlas.structure_descends_from(a, 453)]

        # Filter all areas in rawDf and keep only the ones is the selected sensory system
        boolIdx = [x in targetIds for x in rawDf.index.get_level_values('fine')]
        sensory_df = rawDf.iloc[boolIdx,:]

        # Split areas in primary vs associative
        areatuple = []
        for a in sensory_df.index.get_level_values('fine'):
            layer = self._atlas.get_layer_from_area_id([a])[0]          # Get what layer this area corresponds to
            isPrimaryVis = self._atlas.structure_descends_from(a,385)   # Is a part of primary visual ctx?
            isPrimaryAud = self._atlas.structure_descends_from(a,1002)  # Is a part of primary auditory ctx?
            isPrimarySS = self._atlas.structure_descends_from(a,322)    # Is a part of primary somatosensory ctx?
            isPrimary = isPrimaryVis or isPrimaryAud or isPrimarySS
            if verbose:
                print(f"Area: {self._atlas.ids_to_acronyms([a])[0]}: ", end="")
            if isPrimary:
                areatuple.append(('primary', a, layer))
                if verbose:
                    print("Primary")
            else:
                areatuple.append(('associative', a, layer))
                if verbose:
                    print("Associative")

        multiIdx = pd.MultiIndex.from_tuples(areatuple, names=('sensory', 'fineid', 'layer'))
        sensory_df.index = multiIdx

        # Prepare output dataframe.
        if mergeLayers:
            sensory_df = sensory_df.groupby('sensory').sum()
        else:
            sensory_df = sensory_df.groupby(['sensory','layer']).sum()

        # Calculate extract dataFrames of single measurements
        npnn = sensory_df.xs('n_pnn',axis=1,level='params')
        npv = sensory_df.xs('n_pv',axis=1,level='params')
        ncol = sensory_df.xs('n_colocalized',axis=1,level='params')
        # Calculate metrics
        pvPositive_pnn = ncol.divide(npnn)
        wfaPositive_pv = ncol.divide(npv)

        # Concatenate all metrics in a single output DataFrame
        output_df = pd.concat(
            [pvPositive_pnn,wfaPositive_pv],
            axis=1,
            keys=['pvPositive_pnn','wfaPositive_pv'],
            names=['params', 'mouse']
        )
        # Reorder levels of the columns
        output_df = output_df.reorder_levels(['mouse','params'], axis=1).sort_index(axis=1,level='mouse')

        return output_df
        
    def midDf_to_avgVector(self,
                            mid_df:pd.DataFrame,
                            metric:str,
                            exclude_last: int = 0,
                            ):
        mid = mid_df.copy()
        # Drop coarse index level
        mid.index = mid.index.droplevel('coarse')
        # Extract metric of interest
        mid = mid.xs(metric, axis = 1, level='params')
        # Average across mice
        avg = mid.mean(axis=1)
        # Exclude areas with largest std
        std = mid.std(axis=1)
        std = std.sort_values(ascending = True)
        if exclude_last:
            avg = avg.drop(std.index[-exclude_last:])
        # Rename the series to match the requested metric
        avg.name = metric

        return avg
    
    def flip_corticalLayers(self, cortex_df:pd.DataFrame) -> pd.DataFrame:
        """
        Flips the cortical dataframe generated by the 'regionsDf_to_corticalLayers()' 
        function.

        Parameters
        ----------
        cortex_df: pd.DataFrame
            Cortical data frame generated by the 'regionsDf_to_corticalLayers()' function.
        Returns
        -------
        horizontalCortex_df: pd.DataFrame
            Flipped dataframe with layers as rows and cortical regions as columns
        """
        # Rearrange the dataframe horizontally
        regions = cortex_df.index.get_level_values('mid').unique()
        df_list =  []
        for id in regions:
            df_list.append(cortex_df.loc[id])

        horizontalCortex_df = pd.concat(df_list, axis=1,
            keys= self._atlas.ids_to_acronyms(regions),
            names=['region','params'])
        horizontalCortex_df.sort_index(inplace=True)

        return horizontalCortex_df

    def load_structures_as_df(self, structuresJsonPath):
        """
        Loads the structures json as a DataFrame and performs some processing to make it more usable.
        
        Parameters
        ----------
        structures_json_path:str
            path to the structures.json file, as might be downloaded by a ReferenceSpaceCache object (allen-sdk)

        Returns
        -------
        structuresDf: pd.DataFrame
            indexes are brain areas ids, columns are the keys of structures.json dictionary. Column 'rgb_plotly' is added (expressing rgb_triplet in plotly format)

        """
        # Load the file
        structuresDf = pd.read_json(structuresJsonPath)
        # Set the region ID as the index
        structuresDf = structuresDf.set_index('id') 
        # Create a column with the RGB color in the plotly format e.g., "rgb(100,200,8)"
        rgb_to_strRgb = lambda x: f"rgb({x[0]},{x[1]},{x[2]})"
        structuresDf['rgb_plotly'] = structuresDf['rgb_triplet'].apply(rgb_to_strRgb)
        return structuresDf

    def wholeBrainMetrics(self, rawDf):
        """
        Calculates energy and diffuse fluorescence metrics for the whole brain
        for each animal.
        Parameters
        ----------
        rawDf: pd.DataFrame
            Multiindex df generated by multiIndexDf_from_fineDf.
        Returns
        -------
        brainEnergy: pd.DataFrame
            DataFrame with whole brain energy for each animal
        brainDiffFluo: pd.DataFrame
            DataFrame with whole brain diffuse fluorescence for each animal
        """
        # Calculate metric for the total brain by using the raw Df
        brainAreaMm = rawDf.xs('areaMm2', axis=1, level='params').sum()
        brainAreaPx = rawDf.xs('areaPx', axis=1, level='params').sum()
        brainCellIntensity = rawDf.xs('fluoCellSum', axis=1, level='params').sum()
        brainDiffIntensity = rawDf.xs('diffFluo', axis=1, level='params').sum()
        # Calculate the energy and diffFluo of the entire brain
        brainEnergy = brainCellIntensity.divide(brainAreaMm)
        brainDiffFluo = brainDiffIntensity.divide(brainAreaPx)

        return brainEnergy, brainDiffFluo

    def customCortexOrder(self, oldAtlasNumbers:bool=False):
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
        if oldAtlasNumbers:
            customOrder = [
                184,985,993,353,329,337,345,369,361,182305689,378,              # SS and motor areas
                402,394,409,385,425,533,312782574,312782628,22,                 # Visual Areas
                1011,1002,1027,1018,                                            # Auditory
                44,972,731,                                                     # Media Prefrontal
                723,746,39,48,894,879,886,                                      # Medial Association
                1057,677,104,111,119,541,895,922,918,926,                       # Lateral + ento
            ]
        else:
            customOrder = [
                184,985,993,353,329,337,345,369,361,182305689,378,              # SS and motor areas
                402,394,409,385,425,533,312782574,312782628,312782546,417,      # Visual Areas
                1011,1002,1027,1018,                                            # Auditory
                44,972,731,                                                     # Media Prefrontal
                723,746,39,48,894,879,886,                                      # Medial Association
                1057,677,104,111,119,541,895,922,918,926,                       # Lateral + ento
            ]



        return customOrder

    def add_resolution_column(self, df:pd.DataFrame, region_column:str, resolution = ['mid', 'coarse']):

        df_out = df.copy()

        ids = df_out[region_column].unique().tolist()
        if 'mid' in resolution:
            mid = self._atlas.get_midontology_structures_ids()
            midpar,_ = self._atlas.match_structure_id_lists2(ids, mid)
            in_mid_dict = {m:p for m, p in zip(ids, midpar)}
            df_out['midontologyID'] = df_out[region_column].map(in_mid_dict)
          
        if 'coarse' in resolution:        
            coarse = self._atlas.get_major_divisions_ids()        
            coarsepar,_ = self._atlas.match_structure_id_lists2(ids, coarse)
            in_coarse_dict = {m:p for m, p in zip(ids, coarsepar)}
            df_out['coarseID'] = df_out[region_column].map(in_coarse_dict)

        if ('mid' not in resolution) and ('coarse' not in resolution):
            warn('only mid and coarse resolutions are supported. Df returned as it is', SyntaxWarning)

        return df_out
    
    def multiIndexDf_from_dotsDf(self, dots_df,verbose = True):

        dots_df = dots_df.copy()

        dots_df = dots_df[dots_df['regionID']!= 0]

        #compute total number of pv cells and pnns
        totcells = dots_df[['pv', 'wfa', 'regionID']]\
                                    .groupby(['mouse', 'regionID'])\
                                    .agg(all_pv = ('pv','sum'),all_pnn = ('wfa','sum'))
        #count double-stained cells per brain areas                          
        coloc = dots_df.loc[dots_df['pv'] == dots_df['wfa'], ['pv', 'regionID']]\
                                    .groupby(['mouse', 'regionID'])\
                                    .agg(colocalized = ('pv','sum'))

        #compute total number of pv cells and pnns
        totcells = dots_df[['pv', 'wfa', 'regionID']]\
                                    .groupby(['mouse', 'regionID'])\
                                    .agg(n_pv = ('pv','sum'),n_pnn = ('wfa','sum'))
        totcells.columns.name = 'params'
        totcells = totcells.unstack(['mouse'])\
                            .reorder_levels(['mouse', 'params'], axis= 1)
        totcells = totcells.sort_index(axis=1)

        #count double-stained cells per brain areas                          
        coloc = dots_df.loc[dots_df['pv'] == dots_df['wfa'], ['pv', 'regionID']]\
                                    .groupby(['mouse', 'regionID'])\
                                    .agg(n_colocalized = ('pv','sum'))
        coloc.columns.name = 'params'
        coloc = coloc.unstack(['mouse'])\
                    .reorder_levels(['mouse', 'params'], axis= 1)
        coloc = coloc.sort_index(axis=1)
        coloc

        dotsDf = pd.concat( [totcells, coloc], axis = 1).sort_index(axis = 1)

        id_list = dotsDf.index.tolist()

        midOntologyIDs = self._atlas.get_midontology_structures_ids()
        # Match the fine IDs to the mid IDs
        mid, toDrop_mid = self._atlas.match_structure_id_lists2(id_list, midOntologyIDs)
        if verbose:
            print(f"While matching MID ONTOLOGY structures, {len(toDrop_mid)} structures were dropped:")
            for ID in toDrop_mid:
                print(f"\t- ID: {ID} - Name: {self._atlas.get_structures_by_id([ID])[0]['name']}")
        

        # Adjust dataFrame for COARSE ontology structures
        coarseOntologyIDs = self._atlas.get_major_divisions_ids()
        coarseOntologyIDs.append(1009)     # Add Fiber tracts to the list
        # Match the fine IDs to the coarse IDs
        coarse, toDrop_coarse = self._atlas.match_structure_id_lists2(id_list, coarseOntologyIDs)
        if verbose:
            print(f"While matching COARSE ONTOLOGY structures, {len(toDrop_coarse)} structures were dropped:")
            for ID in toDrop_coarse:
                print(f"\t- ID: {ID} - Name: {self._atlas.get_structures_by_id([ID])[0]['name']}")


        newIndices = [(a,b,c) for a,b,c in zip(coarse,mid,id_list)]
        
        # Mark as regions to drop those that have at lease a None in one of the 3 idexes
        toRemove = [not all(x) for x in newIndices]
        

        # Put the multiindex
        dotsDf.index = pd.MultiIndex.from_tuples(newIndices, names=["coarse","mid","fine"])
        # Remove regions with a None
        if any(toRemove):
            dotsDf.drop(dotsDf[toRemove].index, inplace=True)

        return dotsDf

    def dotsRegionsDf_to_fine(self, rawDf, verbose:bool=False, numberOfCells:bool=False,
        not_to_drop:tp.Sequence[int] = []) -> pd.DataFrame:
        '''
        Converts the multi index dataframe generated by multiIndexDf_from_dotsDf
        in a dataframe of regions at mid ontology levels

        Parameters
        ----------
        rawDf: pd.DataFrame
            Multiindex df generated by multiIndexDf_from_dotsDf.
        not_to_drop: list = []
            List of area IDs to force the algorithm not to drop.
        verbose: bool=False
            Prints out additional information about areas that are dropped
        Returns
        -------
        df: pd.DataFrame
            multiindex dataframe
        '''
        rawDf = rawDf.copy()
        
        # Extract dataFrames of single measurements
        n_pnn = rawDf.xs('n_pnn',axis=1,level='params')
        n_pv = rawDf.xs('n_pv',axis=1,level='params')
        n_colocalized = rawDf.xs('n_colocalized',axis=1,level='params')


        # Calculate measurements
        wfaPositive_pv = n_colocalized.divide(n_pv)*100
        pvPositive_pnn = n_colocalized.divide(n_pnn)*100


        # Concatenate the dataframes in a single multi-animal, multi-measurement
        if numberOfCells:
            temp = pd.concat(
                [n_colocalized, n_pnn, n_pv, wfaPositive_pv,pvPositive_pnn],
                axis=1,
                keys=['n_colocalized','n_pnn','n_pv','wfaPositive_pv','pvPositive_pnn'],
                names=['params', 'mouse']
            )
        else:
            temp = pd.concat(
                [wfaPositive_pv,pvPositive_pnn],
                axis=1,
                keys=['wfaPositive_pv','pvPositive_pnn'],
                names=['params', 'mouse']
            )
        # Reorder levels of the columns
        temp = temp.reorder_levels(['mouse','params'], axis=1).sort_index(axis=1,level='mouse')
        
        # Sort rows based on the areas ontology
        temp['sort'] = self._atlas.ids_to_graph_order(temp.index.get_level_values('fine').to_list())
        temp = temp.sort_values(by=('sort')).drop('sort', axis=1)

        # Drop non informative areas
        midIndex = temp.index.get_level_values('mid').tolist() 
        midlist = self._atlas.get_midontology_structures_ids()
        midlist.extend(not_to_drop)
        areasToDrop = [x for x in midIndex if x not in midlist]
        temp = temp.drop(index=areasToDrop, level='mid')
        if verbose:
            for droppedArea in areasToDrop:
                print(f"- Area ID:{droppedArea} - Name: {self._atlas.get_structures_by_id([droppedArea])[0]['name']} dropped.")

        return temp

    def dotsRegionsDf_to_mid(self, rawDf, verbose:bool=False, numberOfCells:bool=False,
        not_to_drop:tp.Sequence[int] = []) -> pd.DataFrame:
        '''
        Converts the multi index dataframe generated by multiIndexDf_from_fineDf
        in a dataframe of regions at mid ontology levels

        Parameters
        ----------
        rawDf: pd.DataFrame
            Multiindex df generated by multiIndexDf_from_fineDf.
        not_to_drop: list = []
            List of area IDs to force the algorithm not to drop.
        verbose: bool=False
            Prints out additional information about areas that are dropped
        Returns
        -------
        df: pd.DataFrame
            multiindex dataframe
        '''
        rawDf = rawDf.copy()
        grouped = rawDf.groupby(by=['coarse','mid'], axis=0).sum()
        # Extract dataFrames of single measurements
        n_pnn = grouped.xs('n_pnn',axis=1,level='params')
        n_pv = grouped.xs('n_pv',axis=1,level='params')
        n_colocalized = grouped.xs('n_colocalized',axis=1,level='params')


        # Calculate measurements
        wfaPositive_pv = n_colocalized.divide(n_pv)*100
        pvPositive_pnn = n_colocalized.divide(n_pnn)*100


        # Concatenate the dataframes in a single multi-animal, multi-measurement
        if numberOfCells:
            temp = pd.concat(
                [n_colocalized, n_pnn, n_pv, wfaPositive_pv,pvPositive_pnn],
                axis=1,
                keys=['n_colocalized','n_pnn','n_pv','wfaPositive_pv','pvPositive_pnn'],
                names=['params', 'mouse']
            )
        else:
            temp = pd.concat(
                [wfaPositive_pv,pvPositive_pnn],
                axis=1,
                keys=['wfaPositive_pv','pvPositive_pnn'],
                names=['params', 'mouse']
            )
        # Reorder levels of the columns
        temp = temp.reorder_levels(['mouse','params'], axis=1).sort_index(axis=1,level='mouse')
        
        # Sort rows based on the areas ontology
        temp['sort'] = self._atlas.ids_to_graph_order(temp.index.get_level_values('mid').to_list())
        temp = temp.sort_values(by=('sort')).drop('sort', axis=1)

        # Drop non informative areas
        midIndex = temp.index.get_level_values('mid').tolist() 
        midlist = self._atlas.get_midontology_structures_ids()
        midlist.extend(not_to_drop)
        areasToDrop = [x for x in midIndex if x not in midlist]
        temp = temp.drop(index=areasToDrop, level='mid')
        if verbose:
            for droppedArea in areasToDrop:
                print(f"- Area ID:{droppedArea} - Name: {self._atlas.get_structures_by_id([droppedArea])[0]['name']} dropped.")

        return temp

    def dotsRegionsDf_to_coarse(self, rawDf,  verbose:bool=False, 
        not_to_drop:tp.Sequence[int] = []) -> pd.DataFrame:
        '''
        Converts the multi index dataframe generated by multiIndexDf_from_fineDf
        in a dataframe of regions at mid ontology levels

        Parameters
        ----------
        rawDf: pd.DataFrame
            Multiindex df generated by multiIndexDf_from_fineDf.
        not_to_drop: list = []
            List of area IDs to force the algorithm not to drop.
        verbose: bool=False
            Prints out additional information about areas that are dropped
        Returns
        -------
        df: pd.DataFrame
            multiindex dataframe
        '''
        rawDf = rawDf.copy()
        grouped = rawDf.groupby(by=['coarse'], axis=0).sum()
        # Extract dataFrames of single measurements
        n_pnn = grouped.xs('n_pnn',axis=1,level='params')
        n_pv = grouped.xs('n_pv',axis=1,level='params')
        n_colocalized = grouped.xs('n_colocalized',axis=1,level='params')


        # Calculate measurements
        wfaPositive_pv = n_colocalized.divide(n_pv)*100
        pvPositive_pnn = n_colocalized.divide(n_pnn)*100


        # Concatenate the dataframes in a single multi-animal, multi-measurement
        temp = pd.concat(
            [wfaPositive_pv,pvPositive_pnn],
            axis=1,
            keys=['wfaPositive_pv','pvPositive_pnn'],
            names=['params', 'mouse']
        )
        # Reorder levels of the columns
        temp = temp.reorder_levels(['mouse','params'], axis=1).sort_index(axis=1,level='mouse')
        
        # Sort rows based on the areas ontology
        temp['sort'] = self._atlas.ids_to_graph_order(temp.index.get_level_values('coarse').to_list())
        temp = temp.sort_values(by=('sort')).drop('sort', axis=1)

        # Drop non informative areas
        coarseIndex = temp.index.get_level_values('coarse').tolist() 
        coarselist = self._atlas.get_major_divisions_ids()
        coarselist.extend(not_to_drop)
        areasToDrop = [x for x in coarseIndex if x not in coarselist]
        temp = temp.drop(index=areasToDrop)
        if verbose:
            for droppedArea in areasToDrop:
                print(f"- Area ID:{droppedArea} - Name: {self._atlas.get_structures_by_id([droppedArea])[0]['name']} dropped.")

        return temp

    def dots_leaveOneOut_correction(self, rawDf, verbose:bool=False, min_num_pnn:int=2,
        min_num_pv:int=2, min_num_mice:int=3) -> pd.DataFrame:
        """
        For the colocalization at mid and fine resolution, this function implements
        the leave-one-out analysis of colocalization to solve the problem of erratic
        measurements in small areas with very few PNNs or PV cells
        
        Parameters
        ----------
        rawDf: pd.DataFrame
            Multiindex df generated by dotsRegionsDf_to_mid or dotsRegionsDf_to_fine
            with the option numberOfCells=True.
        verbose: bool=False
            Prints out additional information about areas that are dropped
        min_num_pnn:int=2
            Minimum number of PNNs that is acceptable for a brain region
        min_num_pv:int=2
            Minimum number of PV cells that is acceptable for a brain region
        min_num_mice:int=3
            If a brain area has an acceptable number of PNNs and PV cells in at least
            this number of mice then it is included. Otherwise it is dropped.

        Returns
        -------
        processed_df: pd.DataFrame
            multiindex dataframe
        """

        rawDf = rawDf.copy()
        # First, filter out areas that do not respect the requested criteria of
        # minimum number of PNNs or PVs in a minimum number of animals
        valid_pnn = rawDf.xs('n_pnn', axis=1, level='params') > min_num_pnn
        valid_pnn = valid_pnn.sum(axis=1) > min_num_mice

        valid_pv = rawDf.xs('n_pv', axis=1, level='params') > min_num_pv
        valid_pv = valid_pv.sum(axis=1) > min_num_mice

        validAreas = valid_pnn & valid_pv
        
        print(f'Out of a total of {rawDf.shape[0]} areas:\n\t- {(validAreas==True).sum()} were valid\n\t- {(validAreas==False).sum()} were invalid')
        if verbose:
            finestLevel = 'fine' if 'fine' in validAreas.index.names else 'mid'
            invalidIds = validAreas[validAreas==False].index.get_level_values(finestLevel)
            invalidAcronyms = self._atlas.ids_to_acronyms(invalidIds)
            for id, acro in zip(invalidIds,invalidAcronyms):
                print(f'Area: {id} - ({acro}) was invalid and removed.')

        # Select only valid areas for further processing
        rawDf = rawDf.loc[validAreas]

        # Now perform the leaveOneOut analysis
        # ----------------------------------------------------------------------

        mice = rawDf.columns.get_level_values('mouse').unique()
        colnames = rawDf.columns.get_level_values('mouse')
        # Exclude one mouse at a time and calculate aggregated data on all the
        # remaining mice. Each of these steps is an experimental unit.
        dfList = []
        for excludedMouse in mice:
            temp_df = rawDf.loc[:, colnames!=excludedMouse]
            temp = temp_df.groupby('params',axis=1).sum()   # Sum all cells from all mice except one

            # Calculate measurements
            wfaPositive_pv = temp['n_colocalized'].divide(temp['n_pv'])*100
            wfaPositive_pv.name = 'wfaPositive_pv'
            pvPositive_pnn = temp['n_colocalized'].divide(temp['n_pnn'])*100
            pvPositive_pnn.name = 'pvPositive_pnn'

            expUnit_df = pd.concat([pvPositive_pnn, wfaPositive_pv], axis=1)
            dfList.append(expUnit_df)

        processed_df = pd.concat(dfList,
            axis=1,
            names=['excludedMouse','params'],
            keys=mice
            )

        return processed_df

    def dots_to_colocProbability(self, rawDf, categorized_staining:str='pv', n_bins:int=4, verbose:bool=False):
        
        # Parse which intensity measurement to use
        if categorized_staining=='pv':
            intensityVariable = 'fluoMeanPv'
            otherStaining = 'wfa'
        elif categorized_staining=='wfa':
            intensityVariable = 'fluoMeanWfa'
            otherStaining = 'pv'
        
        edges = np.linspace(0,1,n_bins+1)
        # Select all cells of the proper staining
        data = rawDf.loc[rawDf[categorized_staining]==1].copy()

        # Divide all cells in intensity bins
        intensityClass, bins = pd.cut(
            data[intensityVariable],
            edges,
            labels=[str(x+1) for x in range(n_bins)],
            retbins=True
        )

        # Add intensity class
        data['intClass'] = intensityClass

        # Count colocalized cells
        temp = data.loc[data[otherStaining]==1]
        colocCount = temp.groupby(by=['mouse','intClass']).count()
        colocCount = colocCount[categorized_staining]
        colocCount = colocCount.unstack('intClass')


        # Count all cells of the categorized staining
        allCount = data.groupby(by=['mouse','intClass']).count()
        allCount = allCount[categorized_staining]
        allCount = allCount.unstack('intClass')


        # Frequency
        freqDf = colocCount.divide(allCount)
        return freqDf

    def dots_to_sensorySystems(self, rawDf, sensorySystem:str, n_bins:int=4, 
        verbose:bool=False) -> pd.DataFrame:
        '''
        Takes as an input a dataframe of cells and calculates the number of PV cells
        in a given number of intensity classes for each sensory modality

        Parameters
        ----------
        rawDf: pd.DataFrame
            df of all cells, for example generated by allMiceDots().
        sensorySystem: list = []
            "visual", "auditory" or "somatosensory".
        n_bins:int=4
            number of intensity bins to group cells in
        verbose: bool=False
            Prints out additional information about areas that are dropped
        Returns
        -------
        count_df: pd.DataFrame
            dataframe with counts of PV cells for each intensity class and each 
            sensory hierarchy
        '''

        # Define the IDs of areas to keep as part of the selected sensory system
        allIds = rawDf['regionID'].unique()
        if sensorySystem == 'visual':
            # This solves the bug of not including VISrl and VISa in visual areas
            targetIds = [a for a in allIds if 
                self._atlas.structure_descends_from(a, 669) or  # 669: Visual Areas
                self._atlas.structure_descends_from(a, 22)]     # 22: Posterior parietal association areas  
        elif sensorySystem == 'auditory':
            # 247: Auditory Areas
            targetIds = [a for a in allIds if self._atlas.structure_descends_from(a, 247)]
        elif sensorySystem == 'somatosensory':
            # 453: Somatosensory Areas
            targetIds = [a for a in allIds if self._atlas.structure_descends_from(a, 453)]

        # Filter all areas in rawDf and keep only the ones is the selected sensory system
        boolIdx = [x in targetIds for x in rawDf['regionID']]
        sensory_df = rawDf.iloc[boolIdx,:]
        # Reset the index
        sensory_df = sensory_df.reset_index()

        # Split areas in primary vs associative
        areatuple = []
        for a in sensory_df['regionID']:
            isPrimaryVis = self._atlas.structure_descends_from(a,385)   # Is a part of primary visual ctx?
            isPrimaryAud = self._atlas.structure_descends_from(a,1002)  # Is a part of primary auditory ctx?
            isPrimarySS = self._atlas.structure_descends_from(a,322)    # Is a part of primary somatosensory ctx?
            isPrimary = isPrimaryVis or isPrimaryAud or isPrimarySS
            if verbose:
                print(f"Area: {self._atlas.ids_to_acronyms([a])[0]}: ", end="")
            if isPrimary:
                areatuple.append(('primary', a))
                if verbose:
                    print("Primary")
            else:
                areatuple.append(('associative', a))
                if verbose:
                    print("Associative")

        multiIdx = pd.MultiIndex.from_tuples(areatuple, names=('sensory', 'fineid'))
        sensory_df.index = multiIdx


        edges = np.linspace(0,1,n_bins+1)
        # Split PV cells in intensity classes
        intensityClass, bins = pd.cut(
            sensory_df['fluoMeanPv'],
            edges,
            labels=[str(x+1) for x in range(n_bins)],
            retbins=True
        )
        sensory_df['intClass'] = intensityClass

        # Count the number of PV cells for each mouse, hierarchy and class
        count_df = sensory_df.groupby(['mouse', 'sensory', 'intClass']).count()['pv']
        count_df = count_df.unstack('intClass')

        return count_df

if __name__ == '__main__':

    a = Atlas()
