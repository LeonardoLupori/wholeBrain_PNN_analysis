import numpy as np
import pandas as pd

from scipy.stats import zscore, hypergeom, ttest_ind, ttest_rel, spearmanr, pearsonr

from statsmodels.stats.multitest import fdrcorrection, multipletests
from statsmodels.stats.anova import AnovaRM
from sklearn.decomposition import PCA
import typing as tp
import math

def multiple_custom_set_overrepresentation_analysis(exp_list:list,
                                                    custom_sets_dict:dict,
                                                    background:list,
                                                    multiple_test:str = 'fdr_bh', 
                                                    alpha:float = 0.05
                                                    ):
    pvals = []
    enr_ratio =[]
    df_indexes = []
    overlaps = []
    sizes = []
    percentages =[]
    for macroarea in custom_sets_dict.keys():
        if macroarea:
            macroarea_list = custom_sets_dict[macroarea]
            # print(macroarea)
            bg_brain = background
            p, e, o = overrepresentation_analysis(exp_list,macroarea_list, bg_brain)
            sz = len(macroarea_list)

            pvals.append(p)
            enr_ratio.append(e)
            sizes.append(sz)
            overlaps.append(o)
            percentages.append((o/float(sz))*100)
            df_indexes.append(macroarea)

    ora_df = pd.DataFrame(
                                {'set_name':df_indexes,
                                'set_size':sizes,
                                'overlap':overlaps,
                                'percentage_of_term':percentages,
                                'p_value':pvals,
                                'enrichment_ratio':enr_ratio})
    _, fdr,_,_ = multipletests(ora_df['p_value'], method=multiple_test, alpha=alpha)
    # _, fdr = fdrcorrection(area_ora['pvals'])
    ora_df['FDR'] = fdr
    ora_df.sort_values('enrichment_ratio',  inplace=True, ascending=False)
    ora_df = ora_df.reset_index(drop = True)
    return ora_df

def pairwise_comparison(data:pd.DataFrame,
                        dfm:object,
                        sex_dict:dict,
                        metric:str = 'energy',
                        resolution:str = 'mid',
                        variable:str = 'treat',
                        group_1:str = 'CTR',
                        group_2:str = 'HFD', 
                        isolate_group:str = None,
                        multiple_test:str = 'fdr_bh',
                        alpha:float = 0.05):
    '''Performs statistical comparisons between two groups.
        Parameters
        ----------
        data: pd.Dataframe
            multiindex dataframe as returned by AnatomyDataFrameManager
            columns must be multiindex with 3 levels: treat, mouse, params
        metric:str
            selected metric for statistical comparison
        resolution:str
            resolution at which statistical testing is performed. Can be fine, mid, coarse
        variable:str
            variable defining experimental groups (sex, treat)
        group1:str
            name of the control group, used for indexing in the dataframe
        group2:str
            name of the treatment group, used for indexing in the dataframe
        pval_correction:str
            method chosen for correction of multiple testing (FDR or Bonferroni FWER correction)
        multiple_test:str
            method for multiple test correction, see statsmodels.stats.multitest.multipletests
        alpha:float
            alpha significance level set for multiple comparisons (FWER or FDR)
        Returns
        -------
        stat_df:pd.DataFrame
            with index: area_id
            with columns:
                tStat: Student's t
                p: pval resulting from t-test'
                -log10p: negative log of pvalue
                FC: fold change group2 vs. group1
                log2FC: log2 of fold change
                pvals_corr: corrected pvals or q values as returned by statsmodels.stats.multitest.multipletests
                cohen_d: effect size
        '''

    if resolution.lower() == 'fine':
        data = dfm.regionsDf_to_fine(data)
        data.index = data.index.get_level_values(resolution.lower())
    elif resolution.lower() == 'mid':
        data = dfm.regionsDf_to_mid(data)
        data.index = data.index.get_level_values(resolution.lower())
    elif resolution.lower() == 'coarse':
        data = dfm.regionsDf_to_coarse(data)
        data.index = data.index.get_level_values(resolution.lower())
    else:
        raise ValueError('Incorrect resolution. Allowed resolution levels: fine, mid, coarse')
    data = dfm.add_sex_index(data, sex_dict)
    

    if isolate_group:
        try:
            if variable == 'treat':
                data = data.xs(isolate_group, level = 'sex', axis = 1)
            elif variable == 'sex':
                data = data.xs(isolate_group, level = 'treat', axis = 1)
        except KeyError:
            print('isolate_group should refer to elements of the level NOT selected in variable')
    
    data = data.xs(metric, level='params', axis = 1)
    
    stat = list()
    for area, series in data.iterrows():
        group1 = series.xs(group_1, level = variable).dropna()
        group2 = series.xs(group_2, level = variable).dropna()
        # Perform the t-test
        if (len(group2)>2) and (len(group1)>2):
            result = ttest_ind(group1, group2, nan_policy='omit')
            try:
                d = cohend_ind(group1.values, group2.values)
            except RuntimeWarning:
                d = np.nan
            if (group1.mean()>0) and ((group2.mean()/group1.mean())!=np.nan)and ((group2.mean()/group1.mean())!=0):
                fc = group2.mean()/group1.mean()
                logfc = np.log2(float(fc))
            else:
                fc = np.nan
                logfc = np.nan
                
                
        else:
            result = [np.nan, np.nan]
            fc = np.nan
            logfc = np.nan
            d =np.nan

        # Save the results
        temp = {'area_id':area,
            'tStat': result[0],
            'p': result[1], 
            '-log10p':-np.log10(float(result[1])),
            'FC':fc,
            'log2FC':logfc, 
            'cohen_d':d,
                }
        stat.append(temp)
    stat = pd.DataFrame(stat).set_index('area_id')
    stat.dropna(axis=0, how='any', inplace=True)

    _ , fdr,_,_ = multipletests(stat['p'], method = multiple_test, alpha=alpha)
    # _ , fdr  = fdrcorrection(stat['p'])

    stat['FDR'] = fdr
    return stat

def compute_pca(exp_data:pd.DataFrame,
        miceinfo:pd.DataFrame,
        dfm:object,
        metric:str = 'energy',
        resolution:str  ='mid',
        ):

    data = exp_data.copy()
    micedata = miceinfo.copy()
    if resolution.lower() == 'fine':
        data = dfm.regionsDf_to_fine(data)
        data.index = data.index.get_level_values(resolution.lower())
    elif resolution.lower() == 'mid':
        data = dfm.regionsDf_to_mid(data)
        data.index = data.index.get_level_values(resolution.lower())
    elif resolution.lower() == 'coarse':
        data = dfm.regionsDf_to_coarse(data)
        data.index = data.index.get_level_values(resolution.lower())
    else:
        raise ValueError('Incorrect resolution. Allowed resolution levels: fine, mid, coarse')
    data = data.xs(metric, axis = 1, level = 'params')
    data = zscore(data, axis = 1, ddof=1, nan_policy = 'omit')

    pca = PCA(n_components=2)
    comp = pca.fit_transform(data.dropna(axis=0, how =  'any').values.T)
    mice = data.columns.get_level_values('mouse').tolist()
    treat = data.columns.get_level_values('treat').tolist()
    comp = pd.DataFrame(data = {'PC1':comp[:,0],'PC2':comp[:,1], 'treat':treat}, index=mice)
    micedata = micedata.merge(comp, left_index=True,right_index=True, copy = False)
    
    var = pca.explained_variance_ratio_
    return micedata, var

# function to calculate Cohen's d for independent samples
def cohend_ind(d1, d2):
    # calculate the size of samples
    n1=  len(d1)
    n2 = len(d2)
    # calculate the variance of the samples
    s1 = np.var(d1, ddof=1)
    s2 =  np.var(d2, ddof=1)
    # calculate the pooled standard deviation
    s = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    # calculate the means of the samples
    u1 = np.mean(d1)
    u2 =  np.mean(d2)

    # calculate the effect size
    if s != 0:
        d = (u1 - u2) / s
    else:
        d = np.nan

    return d

def cohend_rel(d1,d2):
    
    # calculate the size of samples
    n1=  len(d1)
    n2 = len(d2)
    if n1 != n2:
        raise SyntaxError('This function is supposed to operate in a paired design')

        # calculate the variance of the samples

    
    s = np.std(d1-d2, ddof=1)
    # calculate the means of the samples
    m = np.mean((d1-d2))
    u1 = np.mean(d1)
    u2 =  np.mean(d2)

    # calculate the effect size
    if s != 0:
        dav = m / s
    else:
        dav = np.nan

    return dav

def correlation_with_abagenes(avgseries:pd.Series, genes_df:pd.DataFrame, gene_manager:object):
    gene_id =[]
    gene_acro = []
    gene_name = []
    # Spearmann Correlation
    corr_s = []
    pVal_s =[]
    # Pearson Correlation between log-transformed data
    corr_p = []
    pVal_p =[]
    #coefficient of the log-log regression
    beta = []
    # dropped_rows = []

   
    for x in genes_df.iterrows():
        #concatenate series to perform filtering
        temp = pd.concat([x[1],avgseries], axis=1)

        # replace inf values (if present) with nans 
        temp.replace([np.inf, -np.inf], np.nan, inplace=True)
        # remove areas with unmatched data for either datasets
        temp.dropna(axis=1,how='all', inplace=True)
        temp.dropna(axis=0,how='any', inplace=True)

        if not temp.empty and temp.shape[1] == 2:
            # Spearmann
            correlation_s, p_s = spearmanr(temp.iloc[:,0],temp.iloc[:,1], axis=0)
        else:
            correlation_s = np.nan
            p_s = 1
            print(f'gene {x[0]}, {gene_manager.ids_to_acronyms([int(x[0])])[0]},  impossible to perform Spearman\'s correlation')

        # Pearson
        log_data= temp.loc[~(temp==0).any(axis=1)].copy()
        if (not temp.empty) and (temp.shape[1] == 2) and (log_data.shape[0]>=2):
            log_data.iloc[:,0] = np.log10(log_data.iloc[:,0])
            log_data.iloc[:,1] = np.log10(log_data.iloc[:,1])
            correlation_p, p_p = pearsonr(log_data.iloc[:,0], log_data.iloc[:,1])
            b = np.polyfit(log_data.iloc[:,0], log_data.iloc[:,1], 1)[0]
        else:
            correlation_p = np.nan
            p_p = 1
            b = np.nan
            print(f'gene {x[0]}, {gene_manager.ids_to_acronyms([int(x[0])])[0]},  impossible to perform Pearson\'s correlation')
            
        gene_id.append(x[0])
        corr_s.append(correlation_s)
        pVal_s.append(p_s)
        corr_p.append(correlation_p)
        pVal_p.append(p_p)
        beta.append(b)

    #fill additional gene information
    gene_acro = gene_manager.ids_to_acronyms(gene_id)
    gene_name = gene_manager.ids_to_names(gene_id)
    gene_entrez = gene_manager.ids_to_entrezids(gene_id)


    # correct for multiple testing: Bonferroni
    bonferroni_s = [x*len(pVal_s) for x in pVal_s]
    bonferroni_p = [x*len(pVal_p) for x in pVal_s]
    # fdr Benjamini-Hochberg
    _,bh_p,_,_ = multipletests(pvals=pVal_p, method="fdr_bh")
    _,bh_s,_,_ = multipletests(pvals=pVal_s, method="fdr_bh")


    # Create a unique dataframe to store the statistics
    corr_stat = pd.DataFrame({'gene_id':gene_id,
                        'gene_acronym': gene_acro,
                        'gene_entrez_id':gene_entrez,
                        'gene_name': gene_name, 
                        'corr_spearman':corr_s,
                        'p_spearman':pVal_s,
                        'p_spearman_bonf':bonferroni_s,
                        'p_spearman_fdr': bh_s,
                        'corr_pearson':corr_p,
                        'p_pearson':pVal_p,
                        'p_pearson_bonf':bonferroni_p,
                        'p_pearson_fdr': bh_p,
                        'loglog_reg':beta})

    corr_stat.sort_values(
                            by  = ['corr_spearman', 'p_spearman'] ,
                            axis = 0,
                            ascending = [False,True],
                            inplace = True,
                            ignore_index = True)
  

    return corr_stat

def overrepresentation_analysis(exp_set:tp.Sequence, target_set:tp.Sequence, reference_set:tp.Sequence = None):
    if reference_set != None:

        n_matched = len(set(exp_set).intersection(set(target_set)))
        n_experimental = len(set(exp_set))
        n_target = len(set(target_set).intersection(set(reference_set)))
        n_bg = len(set(reference_set))

        pval = hypergeom.sf(n_matched,n_bg,n_target,n_experimental)

        n_expected = (n_experimental/float(n_bg))*n_target
        enrichment_ratio = n_matched/n_expected
    else:
        raise NotImplementedError('Fisher test not implemented yet!')

    # fraction_matched = (n_matched/float(n_target))*100
    overlap = n_matched
    return pval, enrichment_ratio, overlap

def sensoryCortexByLayers(data:pd.DataFrame, metric:str='energy', printResults:bool=True):
    
    # Melt the input dataframe
    melted = data.xs(metric,axis=1,level='params').melt(ignore_index=False).reset_index()
    # 2-Way Repeated measurements ANOVA
    aovrm2way = AnovaRM(melted, 'value', 'mouse', within=['sensory', 'layer'])
    res2way = aovrm2way.fit()
    # Statistics - post-hoc
    resPostHoc = []
    for layer in melted['layer'].unique():
        temp = melted.loc[melted['layer']==layer]
        # Paired T-test
        stat, pval = ttest_rel(
            temp.loc[temp['sensory']=='primary','value'],
            temp.loc[temp['sensory']=='associative','value'])
        resPostHoc.append((layer, stat, pval))
    # Sidak method for correction
    pAdj = multipletests([x[2] for x in resPostHoc], method='sidak')
    # Store all results in a dataframe
    resPostHoc = pd.DataFrame(resPostHoc,columns=['layer','t','p'])
    resPostHoc['adjP'] = pAdj[1]

    if printResults:
        print(res2way)
        print("Post-hoc (Paired T-test, corrected with Sidak method)")
        for _, x in resPostHoc.iterrows():
            print(f"-Layer {x['layer']:10} t={x['t']:8.4f}  p={x['p']:.4f}  pAdj={x['adjP']:.4f}")

    return res2way, resPostHoc


def sensoryCortexIntClass(data:pd.DataFrame, printResults:bool=True):

    melted = data.melt(ignore_index=False).reset_index()

    # 2-Way Repeated measurements ANOVA
    aovrm2way = AnovaRM(melted, 'value', 'mouse', within=['sensory', 'intClass'])
    res2way = aovrm2way.fit()

    # Statistics - post-hoc
    resPostHoc = []
    for intClass in melted['intClass'].unique():
        temp = melted.loc[melted['intClass']==intClass]
        # Paired T-test
        stat, pval = ttest_rel(
            temp.loc[temp['sensory']=='primary','value'],
            temp.loc[temp['sensory']=='associative','value'])
        resPostHoc.append((intClass, stat, pval))
    # Sidak method for correction
    pAdj = multipletests([x[2] for x in resPostHoc], method='sidak')
    # Store all results in a dataframe
    resPostHoc = pd.DataFrame(resPostHoc,columns=['intClass','t','p'])
    resPostHoc['adjP'] = pAdj[1]

    if printResults:
        print(res2way)
        print("Post-hoc (Paired T-test, corrected with Sidak method)")
        for _, x in resPostHoc.iterrows():
            print(f"-Intensity Class {x['intClass']:10} t={x['t']:8.4f}  p={x['p']:.4f}  pAdj={x['adjP']:.4f}")

    return res2way, resPostHoc


def hiLowWfa_IntClass(data:pd.DataFrame, printResults:bool=True):

    melted = data.melt(ignore_index=False).reset_index()

    # 2-Way Repeated measurements ANOVA
    aovrm2way = AnovaRM(melted, 'value', 'mouse', within=['group', 'intClass'])
    res2way = aovrm2way.fit()

    # Statistics - post-hoc
    resPostHoc = []
    for intClass in melted['intClass'].unique():
        temp = melted.loc[melted['intClass']==intClass]
        # Paired T-test
        stat, pval = ttest_rel(
            temp.loc[temp['group']=='high','value'],
            temp.loc[temp['group']=='low','value'])
        resPostHoc.append((intClass, stat, pval))
    # Sidak method for correction
    pAdj = multipletests([x[2] for x in resPostHoc], method='sidak')
    # Store all results in a dataframe
    resPostHoc = pd.DataFrame(resPostHoc,columns=['intClass','t','p'])
    resPostHoc['adjP'] = pAdj[1]

    if printResults:
        print(res2way)
        print("Post-hoc (Paired T-test, corrected with Sidak method)")
        for _, x in resPostHoc.iterrows():
            print(f"-Intensity Class {x['intClass']:10} t={x['t']:8.4f}  p={x['p']:.4f}  pAdj={x['adjP']:.4f}")

    return res2way, resPostHoc