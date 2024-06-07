'''
SCRIPT FOR PROCESSING THE RESULTS FROM < cortana1782.jar > EXECUTION INTO LOG-FILES 
TO BE USED IN THE GENERATION OF FINAL RESULTS ANALYSIS
'''

import json
import pandas as pd
import pathlib
import statsmodels.api as sm

from collections import Counter
from itertools import chain, combinations_with_replacement
from lifelines import KaplanMeierFitter
from scipy.stats import ttest_ind as t_test


FILE = "cortana_{}_results.txt"

ROOT = "EsmamDS"
ROOT_PATH = str(pathlib.Path(__file__).parent.absolute()).split(ROOT)[0] + ROOT + '/'
DATA_PATH = ROOT_PATH + 'data sets/final data sets/{}_disc.xz'

DATASETS = ['actg320','breast-cancer','cancer','carcinoma','gbsg2','lung','melanoma','mgus2','mgus','pbc','ptc','uis','veteran','whas500']
METRICS = ['num_rules','length','rule_coverage','set_coverage','cr']

VAR_SURV_TIME = 'survival_time'
VAR_SURV_EVENT = 'survival_status'
ALPHA = 0.05

NO_YES = [{'no','yes'}, {'No','Yes'}, {'NO','YES'}]


def __read_data(db_name):
    data_file = DATA_PATH.format(db_name)
    json_file = data_file.split('.')[0] + '_dtypes.json'
    with open(json_file, 'r') as f:
        dtypes = json.load(f)
    db = pd.read_csv(data_file, delimiter=',', header=0, index_col=False, compression='xz', dtype=dtypes)

    for attr in list(db.columns):
        if set(db[attr].unique()) in NO_YES:
            db[attr] = db[attr].replace(regex=[r'no', r'No', r'NO'], value='0')
            db[attr] = db[attr].replace(regex=[r'yes', r'Yes', r'YES'], value='1')
    return db

def __calc_CR(cases_count, size):
    hat = sum(cases_count.values())/size
    _cr = list(map(lambda c: abs(c - hat)/hat , cases_count.values()))
    return sum(_cr)/size

def __get_antecedent(series_val):
    list_of_items = [x.split(' = ') for x in series_val]
    return dict(list_of_items)

def __get_jaccard_matrix(db_results, db, baseline, small=True):
    '''
    Returns a triangular matrix [Rules x Rules] of jaccard index based on the coverage (baseline=cover) 
    or based on the rule's description items (baseline=descript).
    A modified version of the index (that uses the lenght of the smallest set instead of the union of both) can be computed 
    through the parameter small=True.
    The returned matrix is related to a single rule-model
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
    
    antecedents = db_results['Conditions'].apply(lambda x: x.replace("'","")).apply(lambda x: x.split(' AND ')).apply(lambda x: __get_antecedent(x))   
    queries = antecedents.apply(lambda x: ' and '.join([f'`{k}` == {repr(v)}' for k, v in x.items()]))
    cases = queries.apply(lambda x: set(db.query(x).index))
    rules_idx = ['R{}'.format(idx) for idx in range(antecedents.shape[0])]
    matrix = pd.DataFrame(data=None, index=rules_idx, columns=rules_idx)
    
    for r1,r2 in combinations_with_replacement(list(range(len(rules_idx))), 2):
        
        if baseline == 'cover':
            set1 = cases.iloc[r1]
            set2 = cases.iloc[r2]
        else: 
            set1 = set(antecedents.iloc[r1].items())
            set2 = set(antecedents.iloc[r2].items())   
            
        j_index = jaccard_index(set1, set2, small)
        
        idx_r1 = 'R{}'.format(r1)
        idx_r2 = 'R{}'.format(r2)
        if r1==r2: matrix.loc[idx_r1,idx_r1] = j_index
        else:
            matrix.loc[idx_r1,idx_r2] = j_index
            matrix.loc[idx_r2,idx_r1] = j_index

    return matrix.astype(dict.fromkeys(matrix.columns, 'float64')).to_dict()

def __get_rules_pval(db_results, db):
    
    antecedents = db_results['Conditions'].apply(lambda x: x.replace("'","")).apply(lambda x: x.split(' AND ')).apply(lambda x: __get_antecedent(x))   
    queries = antecedents.apply(lambda x: ' and '.join([f'`{k}` == {repr(v)}' for k, v in x.items()]))
    cases = queries.apply(lambda x: set(db.query(x).index))
    rules_idx = ['R{}'.format(idx) for idx in range(antecedents.shape[0])]
    
    rules_pval = {}.fromkeys(rules_idx)
    
    for rule in range(len(rules_idx)):
        idx_rule = 'R{}'.format(rule)
        rules_pval[idx_rule] = {}.fromkeys(['pval_logrank', 'pval_ttest'])
        
        times_sg = db[VAR_SURV_TIME][cases.iloc[rule]]
        events_sg = db[VAR_SURV_EVENT][cases.iloc[rule]]
        times_comp = db[VAR_SURV_TIME]
        events_comp = db[VAR_SURV_EVENT]
        
        # logrank
        times = times_sg.to_list() + times_comp.to_list()
        events = events_sg.to_list() + events_comp.to_list()
        group_id = ['rule'] * times_sg.shape[0] + ['pop'] * times_comp.shape[0]
        _, pval_logrank = sm.duration.survdiff(time=times, status=events, group=group_id)
        
        # t-test
        if times_comp.shape[0]==0: pval_ttest = 1
        else: pval_ttest = t_test(times_sg,times_comp).pvalue
        
        rules_pval[idx_rule]['pval_logrank'] = pval_logrank
        rules_pval[idx_rule]['pval_ttest'] = pval_ttest
    
    return rules_pval

def __get_pval_matrix(db_results, db):
    
    antecedents = db_results['Conditions'].apply(lambda x: x.replace("'","")).apply(lambda x: x.split(' AND ')).apply(lambda x: __get_antecedent(x))   
    queries = antecedents.apply(lambda x: ' and '.join([f'`{k}` == {repr(v)}' for k, v in x.items()]))
    cases = queries.apply(lambda x: set(db.query(x).index))
    rules_idx = ['R{}'.format(idx) for idx in range(antecedents.shape[0])]
    matrix = pd.DataFrame(data=None, index=rules_idx, columns=rules_idx)
    
    for r1,r2 in combinations_with_replacement(list(range(len(rules_idx))), 2):

        times = db[VAR_SURV_TIME][cases.iloc[r1]].to_list() + db[VAR_SURV_TIME][cases.iloc[r2]].to_list()
        events = db[VAR_SURV_EVENT][cases.iloc[r1]].to_list() + db[VAR_SURV_EVENT][cases.iloc[r2]].to_list()
        group_id = ['r1'] * db[VAR_SURV_TIME][cases.iloc[r1]].shape[0] + ['r2'] * db[VAR_SURV_TIME][cases.iloc[r2]].shape[0]
        try:
            _, p_value = sm.duration.survdiff(time=times, status=events, group=group_id)
        except: # apparently, get error when rules are equal
            p_value = 1.0
        
        idx_r1 = 'R{}'.format(r1)
        idx_r2 = 'R{}'.format(r2)
        if r1==r2: matrix.loc[idx_r1,idx_r1] = p_value
        else:
            matrix.loc[idx_r1,idx_r2] = p_value
            matrix.loc[idx_r2,idx_r1] = p_value
    return matrix.astype(dict.fromkeys(matrix.columns, 'float64')).to_dict()

def __get_models_matrix(db_results, db):
    
    # process rules information
    rules = list(range(db_results['Conditions'].shape[0]))
    antecedents = db_results['Conditions'].apply(lambda x: x.replace("'","")).apply(lambda x: x.split(' AND ')).apply(lambda x: __get_antecedent(x))
    queries = antecedents.apply(lambda x: ' and '.join([f'`{k}` == {repr(v)}' for k, v in x.items()]))
    cases = queries.apply(lambda x: set(db.query(x).index))
    
    # kapplan meier fitter to Population
    kmf_pop = KaplanMeierFitter()
    kmf_pop.fit(db[VAR_SURV_TIME], db[VAR_SURV_EVENT], label='pop.', alpha=ALPHA)
    
    # generating dataframe setting
    index = kmf_pop.survival_function_.index.copy()
    columns = ['times', 'population'] + ['R{}'.format(rule) for rule in rules]
    df = pd.DataFrame(columns=columns)
    df.times = index.values
    df.population = kmf_pop.survival_function_.values

    for rule in rules:
        
        # fitting KM model
        kmf = KaplanMeierFitter()
        kmf.fit(db[VAR_SURV_TIME][cases.iloc[rule]], db[VAR_SURV_EVENT][cases.iloc[rule]], label='R{}'.format(rule), alpha=ALPHA)
        
        survival_fnc = kmf.survival_function_.reindex(index)
        survival_fnc.fillna(method='ffill', inplace=True)
        df['R{}'.format(rule)] = survival_fnc.values

    return df.astype(dict.fromkeys(df.columns, 'float64')).to_dict()

def __get_survival_models(folder, save_path):
    '''
    Computes the [Rules X Rules] matrix of jaccard index for both coverage and description.
    The saved log-results refer to Cortana runs over all 14dbs.
    Generated files:
    >> file: cortana_survival-models.json
        A dictionary containing the tables for the survival models [population and rules] of the Cortana results over each dataset. 
        .. file structure: {db_name: survival-models}
        [the survival models are saved as pd.DataFrame.to_dict() and should be read back into pandas obj]
    '''
    models = {}.fromkeys(DATASETS)

    for db_name in DATASETS:
        
        # read dataset
        db = __read_data(db_name)
        # read cortana results
        file = folder + FILE.format(db_name)
        db_results = pd.read_csv(file, sep='\t', index_col=0)

        models[db_name] = __get_models_matrix(db_results, db)
    
    # saving processed files
    save_name = save_path+'cortana_survival-models.json'
    with open(save_name, 'w') as f:
        json.dump(models, f)
    print('..saved: {}'.format(save_name))
    return

def __get_logrank_results(folder, save_path):
    '''
    Computes the [Rules X Rules] matrix of jaccard index for both coverage and description.
    The saved log-results refer to Cortana runs over all 14dbs.
    >> Generated files:
    1. file: cortana_rules-pval.json
        A nested dictionary containing the logrank p-value for each rule in comparison to the population, for each dataset final rule-model 
        .. file structure: {db_name: {rule_idx: p-value}}
    2. file: cortana_pval-matrix.json
        A dictionary containing the [Rules X Rules] matrix of logrank p-values
        .. file structure: {db_name: pval_matrix}
        [saved as pd.DataFrame.to_dict() and should be read back into pandas obj]
    '''
    rules_pval = {}.fromkeys(DATASETS)
    pval_matrix = {}.fromkeys(DATASETS)

    for db_name in DATASETS:
        
        # read dataset
        db = __read_data(db_name)
        # read cortana results
        file = folder + FILE.format(db_name)
        db_results = pd.read_csv(file, sep='\t', index_col=0)

        rules_pval[db_name] = __get_rules_pval(db_results, db)
        pval_matrix[db_name] = __get_pval_matrix(db_results, db)
    
    # saving processed files
    save_name = save_path+'cortana_rules-pval.json'
    with open(save_name, 'w') as f:
        json.dump(rules_pval, f)
    print('..saved: {}'.format(save_name))

    save_name = save_path+'cortana_pval-matrix.json'
    with open(save_name, 'w') as f:
        json.dump(pval_matrix, f)
    print('..saved: {}'.format(save_name))
    return

def __get_jaccard_results(folder, save_path):
    '''
    Computes the [Rules X Rules] matrix of jaccard index for both coverage and description.
    The saved log-results refer to Cortana runs over all 14dbs.
    >> Generated files:
    1. file: cortana_jaccard-cover.json
        A nested dictionary containing the jaccard-matrix based on rule coverage [saved as pd.DataFrame.to_dict() and should be read back into pandas obj]
        .. file structure: {db_name: jaccard_matrix}
    2. file: cortana_jaccard-descript.json
        A nested dictionary containing the jaccard-matrix based on rule description [saved as pd.DataFrame.to_dict() and should be read back into pandas obj]
        .. file structure: {db_name: jaccard_matrix}
    '''
    jac_cover = {}.fromkeys(DATASETS)
    jac_descript = {}.fromkeys(DATASETS)

    for db_name in DATASETS:
        
        # read dataset
        db = __read_data(db_name)
        # read cortana results
        file = folder + FILE.format(db_name)
        db_results = pd.read_csv(file, sep='\t', index_col=0)

        jac_cover[db_name] = __get_jaccard_matrix(db_results, db, 'cover')
        jac_descript[db_name] = __get_jaccard_matrix(db_results, db, 'descript')
    
    # saving processed files
    save_name = save_path+'cortana_jaccard-cover.json'
    with open(save_name, 'w') as f:
        json.dump(jac_cover, f)
    print('..saved: {}'.format(save_name))

    save_name = save_path+'cortana_jaccard-descript.json'
    with open(save_name, 'w') as f:
        json.dump(jac_descript, f)
    print('..saved: {}'.format(save_name))
    return

def __get_metrics(folder, save_path):
    '''
    Compute defined metrics for Cortana results over all 14 datasets and saves a .csv table (dataset X metrics)
    Generated file:
    >> cortana_metrics.csv
    '''
    metrics = {}.fromkeys(DATASETS)

    for db_name in DATASETS:

        db = __read_data(db_name)
        metrics[db_name] = {}.fromkeys(METRICS,0)

        # read cortana results
        file = folder + FILE.format(db_name)
        db_results = pd.read_csv(file, sep='\t', index_col=0)

        # process rules information 
        antecedents = db_results['Conditions'].apply(lambda x: x.replace("'","")).apply(lambda x: x.split(' AND ')).apply(lambda x: __get_antecedent(x))
        queries = antecedents.apply(lambda x: ' and '.join([f'`{k}` == {repr(v)}' for k, v in x.items()]))
        cases = queries.apply(lambda x: set(db.query(x).index))

        all_cases = cases.iloc[0].union(*cases.iloc[1:])                                 # union of all covered cases
        rules_cov = list(chain.from_iterable(cov_list for cov_list in cases.to_list()))  # list combining all covered cases (with repetition)

        # compute metrics
        metrics[db_name]['num_rules'] = db_results.shape[0]
        metrics[db_name]['length'] = db_results['Depth'].mean()
        metrics[db_name]['rule_coverage'] = db_results['Coverage'].sum()/db_results.shape[0]/db.shape[0]
        metrics[db_name]['set_coverage'] = len(all_cases)/db.shape[0]
        metrics[db_name]['cr'] = __calc_CR(Counter(rules_cov), db.shape[0])
    
    save_name = save_path+'cortana_metrics.csv'
    pd.DataFrame(metrics).T.to_csv(save_name)
    print('..saved: {}'.format(save_name))
    return


def run(folder, save_path):
    '''
    Script to extract processed results log-files for DSSD-CBSS algorithm

    Parameters
    ----------
    folder : string
        Folder containing the < cortana_{dataset-name}_results.txt > output files from the algorithm.
    save_path : string
        Path to save the processed results.

    Returns
    -------
    None.
    (This script saves log files on the provided < save_path >)

    '''
    # pipeline
    __get_metrics(folder, save_path)
    __get_jaccard_results(folder, save_path)
    __get_logrank_results(folder, save_path)
    __get_survival_models(folder, save_path)
    return

