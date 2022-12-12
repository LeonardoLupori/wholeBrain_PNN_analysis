import pandas as pd
import typing as tp
import os
# from errors import FileError
import numpy as np

import json
import os
import urllib

pv_markers = ['Pvalb', 'Kcnc1', 'Kcna1', 'Scn1a', 'Syt2', 'Acan']
pnn_markers = ['Acan', 'Hapln1', 'Has3', 'Mmp9', 'Adamts5', 'Pvalb']

def get_pnn_markers():
    """
    Returns a list of gene acronyms for genes that are known to be related to PNNs

    Returns
    ----------
    pnn_markers:list[str]

    """
    pnn_markers = ['Acan', 'Hapln1', 'Has3', 'Mmp9','Adamts5','Pvalb']
    return pnn_markers

def check_markers(corr_stat:pd.DataFrame, staining: str = 'wfa'):

    if staining.lower() == 'wfa':
        ans = corr_stat.loc[corr_stat['gene_acronym'].isin(pnn_markers)]
    elif staining.lower() == 'pv':
        ans = corr_stat.loc[corr_stat['gene_acronym'].isin(pv_markers)]
    else:
        raise ValueError('unsupported staining selected')
    return ans

def save_gene_lists(stat_df:pd.DataFrame, 
                    exp_dataset:str, #e.g. 'wfa'
                    gene_dataset:str, #e.g. 'ABA-ISH'
                    foldername = '',
                    gene_identifier = 'gene_acronym', # or 'entrez_id'
                    corr_method = ['spearman', 'pearson'],
                    positive_corr = True,
                    negative_corr = True,
                    bonferroni = True,
                    fdr = True,
                    alpha_bonf = 0.05,
                    alpha_fdr = 0.1,
                    n_max:int = 1000):

    if len(foldername) >0 and os.path.isdir(foldername) == False:
        raise Exception('Invalid forlder path')
    filestr = gene_dataset +  '_background.txt'
    filename = os.path.join(foldername,filestr)
    stat_df[gene_identifier].to_csv(filename, header=None, index=None, mode='w')
    if 'spearman' in corr_method:
        corr_stat =  stat_df.loc[~stat_df['corr_spearman'].isna()].copy()
        if positive_corr == True:
            if bonferroni == True:
                filestr = exp_dataset + '-' + gene_dataset+  '_posCorr_spearman_bonf.txt'
                filename = os.path.join(foldername, filestr)
                to_save = corr_stat.loc[(corr_stat['corr_spearman']>0)&(corr_stat['p_spearman_bonf']<alpha_bonf)][gene_identifier]
                to_save.iloc[:n_max].to_csv(filename, header=None, index=None, mode='w')
            if fdr == True:
                filestr = exp_dataset + '-' + gene_dataset+  '_posCorr_spearman_fdr.txt'
                filename = os.path.join(foldername, filestr)
                to_save = corr_stat.loc[(corr_stat['corr_spearman']>0)&(corr_stat['p_spearman_fdr']<alpha_fdr)][gene_identifier]
                to_save.iloc[:n_max].to_csv(filename, header=None, index=None, mode='w')
        if negative_corr == True:
            if bonferroni == True:
                filestr = exp_dataset + '-' + gene_dataset+  '_negCorr_spearman_bonf.txt'
                filename = os.path.join(foldername, filestr)
                to_save = corr_stat.loc[(corr_stat['corr_spearman']<0)&(corr_stat['p_spearman_bonf']<alpha_bonf)][gene_identifier]
                to_save.iloc[:n_max].to_csv(filename, header=None, index=None, mode='w')
            if fdr == True:
                filestr = exp_dataset + '-' + gene_dataset+  '_negCorr_spearman_fdr.txt'
                filename = os.path.join(foldername, filestr)
                to_save = corr_stat.loc[(corr_stat['corr_spearman']<0)&(corr_stat['p_spearman_fdr']<alpha_fdr)][gene_identifier]
                to_save.iloc[:n_max].to_csv(filename, header=None, index=None, mode='w')
    if 'pearson' in corr_method:
        corr_stat =  stat_df.loc[~stat_df['corr_pearson'].isna()].copy()
        if positive_corr == True:
            if bonferroni == True:
                filestr = exp_dataset + '-' + gene_dataset+  '_posCorr_pearson_bonf.txt'
                filename = os.path.join(foldername, filestr)
                to_save = corr_stat.loc[(corr_stat['corr_pearson']>0)&(corr_stat['p_pearson_bonf']<alpha_bonf)][gene_identifier]
                to_save.iloc[:n_max].to_csv(filename, header=None, index=None, mode='w')
            if fdr == True:
                filestr = exp_dataset + '-' + gene_dataset+  '_posCorr_pearson_fdr.txt'
                filename = os.path.join(foldername, filestr)
                to_save = corr_stat.loc[(corr_stat['corr_pearson']>0)&(corr_stat['p_pearson_fdr']<alpha_fdr)][gene_identifier]
                to_save.iloc[:n_max].to_csv(filename, header=None, index=None, mode='w')
        if negative_corr == True:
            if bonferroni == True:
                filestr = exp_dataset + '-' + gene_dataset+  '_negCorr_pearson_bonf.txt'
                filename = os.path.join(foldername, filestr)
                to_save = corr_stat.loc[(corr_stat['corr_pearson']<0)&(corr_stat['p_pearson_bonf']<alpha_bonf)][gene_identifier]
                to_save.iloc[:n_max].to_csv(filename, header=None, index=None, mode='w')
            if fdr == True:
                filestr = exp_dataset + '-' + gene_dataset+  '_negCorr_pearson_fdr.txt'
                filename = os.path.join(foldername, filestr)
                to_save = corr_stat.loc[(corr_stat['corr_pearson']<0)&(corr_stat['p_pearson_fdr']<alpha_fdr)][gene_identifier]
                to_save.iloc[:n_max].to_csv(filename, header=None, index=None, mode='w')
        

def print_correlation_summary(corr_stat:pd.DataFrame, alpha_bonf = 0.05, alpha_fdr = 0.1, to_file = False, foldername = None, prefix = ''):
    #counts total number of genes analyzed
    genes_analyzed = corr_stat.shape[0]
    #counts genes positively and negatively correlated (spearman)
    poscor_s = corr_stat.loc[corr_stat['corr_spearman']>0].shape[0]
    negcor_s = corr_stat.loc[corr_stat['corr_spearman']<0].shape[0]
    #counts genes significantly and negatively correlated (pearson)
    poscor_p = corr_stat.loc[corr_stat['corr_pearson']>0].shape[0]
    negcor_p = corr_stat.loc[corr_stat['corr_pearson']<0].shape[0]
    ##
    #SPEARMAN's correlation
    ##
    #counts statistically significant tests - Spearman (raw pval, Bonferroni adjusted pval, FDR)
    s_pval = corr_stat.loc[corr_stat['p_spearman']<0.05].shape[0]
    s_bonf = corr_stat.loc[corr_stat['p_spearman_bonf']<alpha_bonf].shape[0]
    s_fdr = corr_stat.loc[corr_stat['p_spearman_fdr']<alpha_fdr].shape[0]
    #counts genes significantly and positively correlated (spearman)
    poscor_s_pval = corr_stat.loc[(corr_stat['corr_spearman']>0)&(corr_stat['p_spearman']<0.05)].shape[0]
    poscor_s_bonf = corr_stat.loc[(corr_stat['corr_spearman']>0)&(corr_stat['p_spearman_bonf']<alpha_bonf)].shape[0]
    poscor_s_fdr = corr_stat.loc[(corr_stat['corr_spearman']>0)&(corr_stat['p_spearman_fdr']<alpha_fdr)].shape[0]
    #counts genes significantly and positively correlated (spearman)
    negcor_s_pval = corr_stat.loc[(corr_stat['corr_spearman']<0)&(corr_stat['p_spearman']<0.05)].shape[0]
    negcor_s_bonf = corr_stat.loc[(corr_stat['corr_spearman']<0)&(corr_stat['p_spearman_bonf']<alpha_bonf)].shape[0]
    negcor_s_fdr = corr_stat.loc[(corr_stat['corr_spearman']<0)&(corr_stat['p_spearman_fdr']<alpha_fdr)].shape[0]
    ##
    #PEARSON's correlation
    ##
    #counts statistically significant tests - Pearson (raw pval, Bonferroni adjusted pval, FDR)    p_pval = corr_stat.loc[corr_stat['p_pearson']<0.05].shape[0]
    p_pval = corr_stat.loc[corr_stat['p_pearson']<0.05].shape[0]
    p_bonf = corr_stat.loc[corr_stat['p_pearson_bonf']<alpha_bonf].shape[0]
    p_fdr = corr_stat.loc[corr_stat['p_pearson_fdr']<alpha_fdr].shape[0]

    poscor_p_pval = corr_stat.loc[(corr_stat['corr_pearson']>0)&(corr_stat['p_pearson']<0.05)].shape[0]
    poscor_p_bonf = corr_stat.loc[(corr_stat['corr_pearson']>0)&(corr_stat['p_pearson_bonf']<alpha_bonf)].shape[0]
    poscor_p_fdr = corr_stat.loc[(corr_stat['corr_pearson']>0)&(corr_stat['p_pearson_fdr']<alpha_fdr)].shape[0]

    negcor_p_pval = corr_stat.loc[(corr_stat['corr_pearson']<0)&(corr_stat['p_pearson']<0.05)].shape[0]
    negcor_p_bonf = corr_stat.loc[(corr_stat['corr_pearson']<0)&(corr_stat['p_pearson_bonf']<alpha_bonf)].shape[0]
    negcor_p_fdr = corr_stat.loc[(corr_stat['corr_pearson']<0)&(corr_stat['p_pearson_fdr']<alpha_fdr)].shape[0]
    
    #print a summary of the correlation analysis to terminal    
    char = 70
    string = f"Analyzed {genes_analyzed} genes"\
    +"\n"+"*"*(char+5)\
    +"\nSPEARMAN CORRELATION"\
    +"\ngenes significantly correlated: "\
    +"\nraw p-value:".ljust(char)+ f"{s_pval}"\
    +"\nfdr-corrected p-value:".ljust(char)+ f"{s_fdr}"\
    +"\nBonferroni-corrected p-value:".ljust(char)+ f"{s_bonf}"\
    +"\n"+"."*(char+5)\
    +"\ngenes positively correlated: ".ljust(char)+f'{poscor_s}'\
    +"\nraw p-value:".ljust(char)+ f"{poscor_s_pval}"\
    +"\nfdr-corrected p-value:".ljust(char)+ f"{poscor_s_fdr}"\
    +"\nBonferroni-corrected p-value:".ljust(char)+ f"{poscor_s_bonf}"\
    +"\n"+"."*(char+5)\
    +"\ngenes negatively correlated: ".ljust(char)+f'{negcor_s}'\
    +"\nraw p-value:".ljust(char)+ f"{negcor_s_pval}"\
    +"\nfdr-corrected p-value:".ljust(char)+ f"{negcor_s_fdr}"\
    +"\nBonferroni-corrected p-value:".ljust(char)+ f"{negcor_s_bonf}"\
    +"\n"+"*"*(char+5)\
    +"\nPEARSON CORRELATION"\
    +"\ngenes significantly correlated: "\
    +"\nraw p-value:".ljust(char)+ f"{p_pval}"\
    +"\nfdr-corrected p-value:".ljust(char)+ f"{p_fdr}"\
    +"\nBonferroni corrected:".ljust(char)+ f"{p_bonf}"\
    +"\n"+"."*(char+5)\
    +"\ngenes positively correlated: ".ljust(char)+f'{poscor_p}'\
    +"\nraw p-value:".ljust(char)+ f"{poscor_p_pval}"\
    +"\nfdr-corrected p-value:".ljust(char)+ f"{poscor_p_fdr}"\
    +"\nBonferroni-corrected p-value:".ljust(char)+ f"{poscor_p_bonf}"\
    +"\n"+"."*(char+5)\
    +"\ngenes negatively correlated: ".ljust(char)+f'{negcor_p}'\
    +"\nraw p-value:".ljust(char)+ f"{negcor_p_pval}"\
    +"\nfdr-corrected p-value:".ljust(char)+ f"{negcor_p_fdr}"\
    +"\nBonferroni-corrected p-value:".ljust(char)+ f"{negcor_p_bonf}"
    print(string)


    if to_file == True:
        try:
            if foldername != None and os.path.isdir(foldername):
                filename = os.path.join(foldername, prefix + "_correlationResult.txt")
            else:
                filename = prefix + "_correlationResult.txt"
                print(f'Foldername ({foldername}) not valid. File saved in current directory')
            with open(filename, "w") as text_file:
                text_file.write(string)
        except PermissionError:
            print("unable to save to file, correlationResult.txt already present!")

# Make a query to the API via a URL.
def QueryAPI(url,verbose=False):
    start_row = 0
    num_rows = 2000
    total_rows = -1
    rows = []
    done = False

    while not done:
        pagedUrl = url + '&start_row=%d&num_rows=%d' % (start_row,num_rows)
        
        if verbose:
            print(pagedUrl)
            
        source = urllib.request.urlopen(pagedUrl).read()

        response = json.loads(source)
        rows += response['msg']
        
        if total_rows < 0:
            total_rows = int(response['total_rows'])

        start_row += len(response['msg'])

        if start_row >= total_rows:
            done = True

    return rows


def query_gene_data():
    # graph_ID = 1
    product_ID = 1

    BASE = "http://api.brain-map.org/api/v2/data"
    target = "/SectionDataSet/query.json?"
    criteria = "criteria=[failed$eq'false'][expression$eq'true'],"
    products = "products[id$eq{}]".format(product_ID)
    inclusions = "&include=genes"
    exclusions = "&except=blue_channel,delegate,expression,failed,failed_facet,green_channel," + \
                    "name,qc_date,red_channel,rnaseq_design_id,sphinx_id,storage_directory,weight"

    url = BASE + target + criteria + products + inclusions + exclusions

    result = QueryAPI(url,verbose=False)
    genesDict = [g['genes'][0] for g in result]
    genesDict = [dict(t) for t in {tuple(d.items()) for d in genesDict}]
    return genesDict


class GeneManager():
    def __init__(self, path:str = None):
        if path !=  None:
            if not os.path.isfile(path):
                raise ValueError(path,msg = 'Specified path does not correspond to any existing file')
            
            if not path.lower().endswith('.json'):
                raise ValueError(path, msg='Specified file must have .json extension')
            try:
                self.genes_df = pd.read_json(path)
                with open(path, "r") as gene_file:
                    self.genes_dict = json.load(gene_file)
                # self.genes_dict = self.genes_df.set_index('gene_id').to_dict()
                # self.genes_dict = [d for d in self.genes_dict if d['id'] != 135754]
            except (ValueError, FileNotFoundError):
                print('Error: specified path not valid!')
        elif path == None:
            try:
                self.genes_df = pd.read_json('genes.json')
                # self.genes_dict = self.genes_df.set_index('gene_id').to_dict()
                with open('genes.json', "r") as gene_file:
                    self.genes_dict = json.load(gene_file)
                # self.genes_dict = [df[1].to_dict() for df in self.genes_df.iterrows()]
                # self.genes_dict = [d for d in self.genes_dict if d['id'] != 135754]
            except (ValueError, FileNotFoundError):
                print('genes.json not found... Dowloading gene data from the ABA server...')
                self.genes_dict = query_gene_data()
                # self.genes_dict = [d for d in self.genes_dict if d['id'] != 135754]
                self.genes_df = pd.DataFrame(self.genes_dict)
                with open('genes.json', "w") as gene_file:
                    json.dump(self.genes_dict, gene_file, indent=4)
                # self.genes_df.to_json('genes.json', index_label=False)
                print('done!\n')
        node_id_cb = lambda s:s['id']
        self._genes = { node_id_cb(n):n for n in self.genes_dict }



    #converts list of acronyms into ids
    def acronyms_to_ids(self, acronyms:tp.Sequence[str]):
        to_fn = lambda x: x['id']
        ids = self.nodes_by_property('acronym', acronyms, to_fn)
        return ids

    #converts list of ids into acronyms
    def ids_to_acronyms(self, ids:tp.Sequence[int]):
        to_fn = lambda x: x['acronym']
        acronyms = self.nodes_by_property('id', ids, to_fn)
        return acronyms

    def ids_to_entrezids(self, ids:tp.Sequence[int]):
        to_fn = lambda x: str(x['entrez_id']) if x['entrez_id'] != None else None
        entrez = self.nodes_by_property('id', ids, to_fn)
        return entrez

    def acronyms_to_names(self, ids:tp.Sequence[str]):
        to_fn = lambda x: x['original_name']
        names = self.nodes_by_property('acronym', ids, to_fn)
        return names

    def ids_to_names(self, ids:tp.Sequence[str]):
        to_fn = lambda x: x['original_name']
        names = self.nodes_by_property('id', ids, to_fn)
        return names

    def get_gene_by_id(self, gene_ids:tp.Sequence[int]):
        return self.nodes_by_property('id', gene_ids)

    def get_gene_by_acronym(self, gene_acro:tp.Sequence[str]):
        return self.nodes_by_property('acronym', gene_acro)


    def nodes_by_property(self, key, values, to_fn=None):
        '''Get nodes by a specified property

        Parameters
        ----------
        key : hashable or function
            The property used for lookup. Should be unique. If a function, will 
            be invoked on each node.
        values : list
            Select matching elements from the lookup.
        to_fn : function, optional
            Defines the outputs, on a per-node basis. Defaults to returning 
            the whole node.
  
        Returns
        -------
        list : 
            outputs, 1 for each input value.

        '''

        if to_fn is None:
            to_fn = lambda x: x

        if not callable( key ):
            from_fn = lambda x: x[key]
        else:
            from_fn = key

        value_map = self.value_map( from_fn, to_fn )#        names = self.nodes_by_property('id', ids, to_fn)
        return [ value_map[vv] for vv in values ]

    def value_map(self, from_fn, to_fn):
        '''Obtain a look-up table relating a pair of node properties across 
        nodes
        
        Parameters
        ----------
        from_fn : function | node dict => hashable value
            The keys of the output dictionary will be obtained by calling 
            from_fn on each node. Should be unique.
        to_fn : function | node_dict => value
            The values of the output function will be obtained by calling 
            to_fn on each node.
            
        Returns
        -------
        dict :
            Maps the node property defined by from_fn to the node property 
            defined by to_fn across nodes.

        '''
        
        vm = {}
        
        for node in list(self._genes.values()):
            key = from_fn(node)
            value = to_fn(node)
            
            # if key in vm:
            #     print('from_fn is not unique across nodes. \nCollision between {0} and {1}.'.format(value, vm[key]))
            #     vm2[key] = value
            #     continue
            #     #                    'Collision between {0} and {1}.'.format(value, vm[key])))
            #     # # raise RuntimeError('from_fn is not unique across nodes. '
            #     #                    'Collision between {0} and {1}.'.format(value, vm[key]))    
            vm[key] = value
  
        return vm



    def save_csv(self, folder:str = ''):
        if isinstance(folder, str) == False:
            raise TypeError('Folder path must be a string')
        if os.path.isdir(folder):
            filename = os.path.join(folder, 'genesDataFrame.csv')
            # self.genes_dict = query_gene_data()
            # self.genes_df = pd.DataFrame(self.genes_dict)
            try:
                self.genes_df.to_csv(filename,index=False)
            except PermissionError:
                print(f'Error: no wtriting permission for {filename}')
        elif len(folder) == 0:
            filename =  'genesDataFrame.csv'
            # self.genes_dict = query_gene_data()
            # self.genes_df = pd.DataFrame(self.genes_dict)
            try:
                self.genes_df.to_csv(filename,index=False)
            except PermissionError:
                print(f'Error: no wtriting permission for {filename}')        
        elif os.path.isfile(folder):
            filename = folder
            # self.genes_dict = query_gene_data()
            # self.genes_df = pd.DataFrame(self.genes_dict)
            try:
                self.genes_df.to_csv(filename, index=False)
            except PermissionError:
                print(f'Error: no wtriting permission for {filename}')
            



if __name__ == '__main__':

    g = GeneManager()