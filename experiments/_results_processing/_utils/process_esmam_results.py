'''
SCRIPT TO PROCESS THE ESMAM(-DS) ALGORITHM _LOG.JSON FILE OF STATISTICAL (REPEATED) EXECUTIONS OVER ALL 14 DATASETS
- SAVES PROCESSED LOG FILES TO BE USED IN THE GENERATION OF THE FINAL RESULTS FOR ANALYSIS (THROUGH THE SCRIPT <generate_results>)

'''

import json
import os
import pandas as pd
import pathlib
import statsmodels.api as sm
import zipfile

from collections import Counter
from itertools import chain, combinations_with_replacement
from lifelines import KaplanMeierFitter
from scipy.stats import ttest_ind as t_test

# Global variables
LOG_FILE_NAME = '{}_{}_exp{}_log.json'

ROOT = "EsmamDS"
ROOT_PATH = str(pathlib.Path(__file__).parent.absolute()).split(ROOT)[0] + ROOT + '/'
DATA_PATH = ROOT_PATH + 'data sets/final data sets/{}_disc.xz'

METRICS = ['num_rules','length','rule_coverage','set_coverage']
DATASETS = ['actg320','breast-cancer','cancer','carcinoma','gbsg2','lung','melanoma','mgus2','mgus','pbc','ptc','uis','veteran','whas500']

VAR_SURV_TIME = 'survival_time'
VAR_SURV_EVENT = 'survival_status'
RUNS = 30
ALPHA = 0.05


def __read_data(db_name):

    data_file = DATA_PATH.format(db_name)
    json_file = data_file.split('.')[0] + '_dtypes.json'
    with open(json_file, 'r') as f:
        dtypes = json.load(f)
    return pd.read_csv(data_file, delimiter=',', header=0, index_col=False, compression='xz', dtype=dtypes)

def __get_terms(rule):
    terms = []
    for attr in rule['antecedent'].keys():
        if isinstance(rule['antecedent'][attr], list):
            for v in rule['antecedent'][attr]:
                terms.append((attr, v))
        else:
            terms.append((attr,rule['antecedent'][attr]))
    return set(terms)


def __calc_CR(cases_count, size):
    hat = sum(cases_count.values())/size
    _cr = list(map(lambda c: abs(c - hat)/hat , cases_count.values()))
    return sum(_cr)/size


def __get_survmodels_matrix(version_log, db_name, exp):
    '''
    Returns a table with the survival models [population and rules].
    The returned table is related to a single rule-model
    [the survival models are saved as pd.DataFrame.to_dict() and should be read back into pandas obj]
    '''
    db = __read_data(db_name)
    rules = version_log['model']
    
    # kapplan meier fitter to Population
    kmf_pop = KaplanMeierFitter()
    kmf_pop.fit(db[VAR_SURV_TIME], db[VAR_SURV_EVENT], label='pop.', alpha=ALPHA)
    
    # generating dataframe setting
    index = kmf_pop.survival_function_.index.copy()
    columns = ['times', 'population'] + list(rules.keys())
    df = pd.DataFrame(columns=columns)
    df.times = index.values
    df.population = kmf_pop.survival_function_.values

    for rule in rules:
        
        # fitting KM model
        kmf = KaplanMeierFitter()
        kmf.fit(db[VAR_SURV_TIME][rules[rule]['cases']], db[VAR_SURV_EVENT][rules[rule]['cases']], label=rule, alpha=ALPHA)
        
        survival_fnc = kmf.survival_function_.reindex(index)
        survival_fnc.fillna(method='ffill', inplace=True)
        df[rule] = survival_fnc.values

    return df.astype(dict.fromkeys(df.columns, 'float64')).to_dict()

def __get_pval_matrix(version_log, db_name, exp):
    '''
    Returns a triangular matrix [Rules x Rules] of logrank p-values between combinations of rules.
    The returned matrix is related to a single rule-model (a single algorithm run - experiment)
    return: pd.DataFrame
    '''
    
    db = __read_data(db_name)
    rules = version_log['model']
    rules_idx = list(rules.keys())
    
    matrix = pd.DataFrame(data=None, index=rules_idx, columns=rules_idx)
    
    for r1,r2 in combinations_with_replacement(rules_idx, 2):

        times = db[VAR_SURV_TIME][rules[r1]['cases']].to_list() + db[VAR_SURV_TIME][rules[r2]['cases']].to_list()
        events = db[VAR_SURV_EVENT][rules[r1]['cases']].to_list() + db[VAR_SURV_EVENT][rules[r2]['cases']].to_list()
        group_id = ['r1'] * db[VAR_SURV_TIME][rules[r1]['cases']].shape[0] + ['r2'] * db[VAR_SURV_TIME][rules[r2]['cases']].shape[0]
        try:
            _, p_value = sm.duration.survdiff(time=times, status=events, group=group_id)
        except: # apparently, get error when rules are equal
            p_value = 1.0
        
        if r1==r2: matrix.loc[r1,r1] = p_value
        else:
            matrix.loc[r1,r2] = p_value
            matrix.loc[r2,r1] = p_value
    return matrix.astype(dict.fromkeys(matrix.columns, 'float64')).to_dict()

def __get_pval_rules(version_log, db_name, exp, baseline):
    '''
    Return a dictionary with keys as the representatives of each rule in the specific [db-exp] model and the values 
    are the logrank p-value of the rules models in comparison to population
    return: dictionary
    '''
    
    db = __read_data(db_name)
    rules = version_log['model']
    rules_idx = list(rules.keys())
    
    rules_pval = {}.fromkeys(rules_idx)
    
    for rule in rules:
        
        rules_pval[rule] = {}.fromkeys(['pval_logrank', 'pval_ttest'])
        
        times_sg = db[VAR_SURV_TIME][rules[rule]['cases']]
        events_sg = db[VAR_SURV_EVENT][rules[rule]['cases']]
        group_sg = ['rule']*times_sg.shape[0]
        
        if baseline=='complement':
            cases_cpm = set(rules[rule]['cases'])^set(db.index)
            times_comp = db[VAR_SURV_TIME][cases_cpm]
            events_comp = db[VAR_SURV_EVENT][cases_cpm]
            group_comp = ['cpm'] * times_comp.shape[0]            
        elif baseline=='population':
            times_comp = db[VAR_SURV_TIME]
            events_comp = db[VAR_SURV_EVENT]
            group_comp = ['pop'] * times_comp.shape[0]
        
        # logrank
        times = times_sg.to_list() + times_comp.to_list()
        events = events_sg.to_list() + events_comp.to_list()
        group_id = group_sg + group_comp
        _, pval_logrank = sm.duration.survdiff(time=times, status=events, group=group_id)
        
        # t-test
        if times_comp.shape[0]==0: pval_ttest = 1
        else: pval_ttest = t_test(times_sg,times_comp).pvalue
        
        rules_pval[rule]['pval_logrank'] = pval_logrank
        rules_pval[rule]['pval_ttest'] = pval_ttest
        
    return rules_pval

def __get_jaccard(version_log, db, exp, baseline, small=True):
    '''
    Returns a triangular matrix [Rules x Rules] of jaccard index based on the coverage (baseline=cover) 
    or based on the rule's description items (baseline=descript).
    A modified version of the index (that uses the lenght of the smallest set instead of the union of both) can be computed 
    through the parameter small=True.
    The returned matrix is related to a single rule-model (a single algorithm run - experiment)
    return: pd.DataFrame
    '''

    def jaccard_index(set1, set2, small):
        minor = min([len(set1),len(set2)])
        union = set1.union(set2)
        intersection = set1.intersection(set2)
        if small:
            return len(intersection)/minor
        else:
            return len(intersection)/len(union) 
    
    rules = version_log['model']
    rules_idx = list(rules.keys())
    matrix = pd.DataFrame(data=None, index=rules_idx, columns=rules_idx)
    
    for r1,r2 in combinations_with_replacement(rules_idx, 2):
        
        if baseline == 'cover':
            set1 = set(rules[r1]['cases'])
            set2 = set(rules[r2]['cases'])
        else: 
            set1 = __get_terms(rules[r1])
            set2 = __get_terms(rules[r2])    
            
        j_index = jaccard_index(set1, set2, small)

        if r1==r2: matrix.loc[r1,r1] = j_index
        else:
            matrix.loc[r1,r2] = j_index
            matrix.loc[r2,r1] = j_index

    return matrix.astype(dict.fromkeys(matrix.columns, 'float64')).to_dict()

def __get_metrics_exp(version_log, db, exp):
    '''
    Returns a dictionary of the defined analysis metrics: {metric_name: metric_value}.
    The metrics are related to a single algorithm run (experiment).
    returns: dictionary
    '''
    
    metrics_exp = version_log['metrics'].copy()  # read a dictionary
    for k in metrics_exp.keys():
        if k not in METRICS:
            del metrics_exp[k]
    metrics_exp['cr'] = None

    # calc CR
    ruleset = version_log['model']
    rules_cov = list(chain.from_iterable(rule['cases'] for rule in ruleset.values()))
    db_size = int(version_log['run']['data_shape'][0])
    metrics_exp['cr'] = __calc_CR(Counter(rules_cov), db_size)
    
    return metrics_exp
    

def __generate_log_results(version, version_folder, metric, func, save_path, **kwargs):
    '''
    Generates processed log-results to be used for generating the final results for analysis.
    The saved log-results are related to the hole statistical proceadure, 
    i.e. they refer to all algorithms runs (14dbs - 30 experiments each).
    >> Generated files:
    1. file: <algorithm-baseline>_metrics.csv 
        A table <dataset(14x30) x metrics> with the performance of the specific <algorithm-baseline> on the N-experiments.
    2. file: <algorithm-baseline>_jaccard-cover.json
        A nested dictionary containing the jaccard-matrix [based on rule coverage] (Rules x Rules) of each algorithm final model.
        The jaccard-matrix is saved as pd.DataFrame.to_dict() and should be read back into pandas obj.
        .. file structure: {db_name: exp_N: jaccard_matrix}
    3. file: <algorithm-baseline>_jaccard-descript.json
        A nested dictionary containing the jaccard-matrix [based on rule description] (Rules x Rules) of each algorithm final model.
        The jaccard-matrix is saved as pd.DataFrame.to_dict() and should be read back into pandas obj.
        .. file structure: {db_name: exp_N: jaccard_matrix}
    4. file: <algorithm-baseline>_rules-pval.json
        A nested dictionary containing the logrank p-value [in comparison to population] of each rule in the final rule-models
        .. file structure: {db_name: exp_N: {Rule_ID: p-value}}
    5. file: <algorithm-baseline>_pval-matrix.json
        A nested dictionary containing the logrank p-value matrix (Rules x Rules) of each combination of rules in the final rule-models
        The pval-matrix is saved as pd.DataFrame.to_dict() and should be read back into pandas obj.
        .. file structure: {db_name: exp_N: pvalue_matrix}
    6. file: <algorithm-baseline>_survival-models.json
        A nested dictionary containing a table with the survival models for population and each rule in the final rule-models
        The survival table is saved as pd.DataFrame.to_dict() and should be read back into pandas obj.
        .. file structure: {db_name: exp_N: survival-table}
    '''

    dic_db_exp = {}
    
    for db in DATASETS:
        dic_exp = {}
        
        for exp in range(RUNS):
            
            log_file = version_folder + LOG_FILE_NAME.format(version, db, exp)
            if os.path.exists(log_file): # exists a .json file
                with open(log_file, 'r') as f:
                    version_log = json.load(f)
            else: # dont exist a .json file >> read (.json).zip file
                zip_path = log_file.split('.')[0] + '.zip'
                with zipfile.ZipFile(zip_path) as z:
                    for filename in z.namelist():
                        with z.open(filename) as f:  
                            data = f.read()  
                            dic_exp = json.loads(data)
            
            dic_exp[exp] = func(version_log, db, exp, **kwargs)

        dic_db_exp[db] = dic_exp.copy()
    
    ## SAVES METRICS TABLE
    if metric == 'metrics':
        save_name = save_path+'{}_metrics.csv'.format(version)
        df = pd.concat([pd.DataFrame(dic_db_exp[db].values(), index=dic_db_exp[db].keys())
                        for db in dic_db_exp.keys()],
                       axis=0, keys=list(dic_db_exp.keys()), names=['Dataset','expID']).reset_index()
        df.to_csv(save_name, index=False)
        print('..saved: {}'.format(save_name))
    
    ## SAVES JACCARD-COVER FILE
    if metric == 'jaccard_c':
        save_name = save_path+'{}_jaccard-cover.json'.format(version)
        with open(save_name, 'w') as f:
            json.dump(dic_db_exp, f)
        print('..saved: {}'.format(save_name))
    
    ## SAVES JACCARD-DESCRIPTION FILE
    if metric == 'jaccard_d':
        save_name = save_path+'{}_jaccard-descript.json'.format(version)
        with open(save_name, 'w') as f:
            json.dump(dic_db_exp, f)
        print('..saved: {}'.format(save_name))
    
    ## SAVES PVALUE-RULES FILE
    if metric == 'pval_r':
        save_name = save_path+'{}_rules-pval.json'.format(version)
        with open(save_name, 'w') as f:
            json.dump(dic_db_exp, f)
        print('..saved: {}'.format(save_name))
     
    ## SAVES PVALUE-MATRIX FILE
    if metric == 'pval_m':
        save_name = save_path+'{}_pval-matrix.json'.format(version)
        with open(save_name, 'w') as f:
            json.dump(dic_db_exp, f)
        print('..saved: {}'.format(save_name))
    
    ## SAVES PVALUE-MATRIX FILE
    if metric == 'surv_model':
        save_name = save_path+'{}_survival-models.json'.format(version)
        with open(save_name, 'w') as f:
            json.dump(dic_db_exp, f)
        print('..saved: {}'.format(save_name))
    
    return

def run(algorithm, folder, save_path):
    '''
    Script to extract processed results log-files for Esmam class algorithms

    Parameters
    ----------
    algorithm : string
        Name of the Esmam-class algorithm to process the results.
        [EsmamDS-pop, EsmamDS-cpm, Esmam-pop, Esmam-cpm]
    folder : string
        Folder containing the (statistical evaluation) < {algorithm}_{dataset-name}_exp{N}_log.json > 
        output files from the algorithm.
        This folder should contain the Esmam-class algorithm output files, including the .json file 
        for all 14 data sets and 30 runs.
    save_path : string
        Path to save the processed results.

    Returns
    -------
    None.
    (This script saves log files on the provided < save_path >)

    '''
        
    if algorithm[-3:] == 'pop':
        baseline = 'population'
    if algorithm[-3:] == 'cpm':
        baseline = 'complement'
    
    __generate_log_results(algorithm, folder, metric='metrics', func=__get_metrics_exp)
    __generate_log_results(algorithm, folder, metric='jaccard_c', func=__get_jaccard, **{'baseline':'cover'})
    __generate_log_results(algorithm, folder, metric='jaccard_d', func=__get_jaccard, **{'baseline':'descrpt'})
    __generate_log_results(algorithm, folder, metric='pval_r', func=__get_pval_rules, **{'baseline':baseline})
    __generate_log_results(algorithm, folder, metric='pval_m', func=__get_pval_matrix)
    __generate_log_results(algorithm, folder, metric='surv_model', func=__get_survmodels_matrix)
    return


