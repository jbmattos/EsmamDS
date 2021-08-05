"""
SCRIPT TO GENERATE THE RESULTS OF RULE-SET COMPARISON BETWEEN APPROACHED ALGORITHMS. 
"""

import json
import os
import pandas as pd
import pathlib
import re
import zipfile

from lifelines.statistics import logrank_test

# Global variables
ROOT = "EsmamDS"
ROOT_PATH = str(pathlib.Path(__file__).parent.absolute()).split(ROOT)[0] + ROOT + '/'
DATA_PATH = ROOT_PATH + 'data sets/final data sets/{}_disc.xz'
READ_PATH = ROOT_PATH+'experiments/_results_processing/_algorithms_output_files/'
SAVE_PATH = ROOT_PATH+'experiments/_results_processing/_processed_output_files/_inter-set_similarity/'
SAVE_FILE = '{}_{}'

ALG_BASE = {'population': 'EsmamDS-pop',
            'complement': 'EsmamDS-cpm'
           }
ALGORITHMS = {'population': ['Esmam-pop', 'BS-EMM-pop', 'BS-SD-pop', 'DSSD-CBSS'],
              'complement': ['Esmam-cpm', 'BS-EMM-cpm', 'BS-SD-cpm','LR-Rules']
             }
ALG_FILE = {'Esmam-pop': READ_PATH + 'esmam-pop/Esmam-pop_{}_exp{}_log.json',
            'Esmam-cpm': READ_PATH + 'esmam-cpm/Esmam-cpm_{}_exp{}_log.json',
            'EsmamDS-pop': READ_PATH + 'esmamds-pop/EsmamDS-pop_{}_exp{}_log.json',
            'EsmamDS-cpm': READ_PATH + 'esmamds-cpm/EsmamDS-cpm_{}_exp{}_log.json',
            'BS-EMM-pop': READ_PATH + 'bs-emm-pop/{}_result.csv',
            'BS-EMM-cpm': READ_PATH + 'bs-emm-cpm/{}_result.csv',
            'BS-SD-pop': READ_PATH + 'bs-sd-pop/{}_result.csv',
            'BS-SD-cpm': READ_PATH + 'bs-sd-cpm/{}_result.csv',
            'DSSD-CBSS': READ_PATH + 'dssd-cbss/cortana_{}_results.txt',
            'LR-Rules': READ_PATH + 'lr-rules/{}_ruleSet.txt'}

LOG_FILE = {'descr': 'interset_descrSimilarity.json',
            'cover': 'interset_coverSimilarity.json',
            'model': 'interset_modelSimilarity.json',
           }

BASELINES = ['population', 'complement']
DATASETS = ['actg320','breast-cancer','cancer','carcinoma','gbsg2','lung','melanoma','mgus2','mgus','pbc','ptc','uis','veteran','whas500']
VAR_SURV_TIME = 'survival_time'
VAR_SURV_EVENT = 'survival_status'
RUNS = 30
ALPHA = 0.05

VAL_CONVERT = {'0': r'^[nN][oO]$',
               '1': r'^[yY][eE][sS]$'}

class RuleSet():
    
    def __init__(self, pdseries_rules_set_items, db=None):
        self.name = pdseries_rules_set_items.name
        self.__ruleset = self.__set(pdseries_rules_set_items)
        self.set_cases(db)
    
    def __set(self, s):
        if s.shape[0]==0:
            return pd.Series([])
        
        rules = []
        rules.append(s.iloc[0])
        for idx in range(s.shape[0]):
            if s.iloc[idx] not in rules:
                rules.append(s.iloc[idx])
        return pd.Series(rules)
    
    def __adjust_values(self, db):
        
        for idx in range(len(self.__ruleset.index)): # iterates over all rules
            terms = self.__ruleset.iloc[idx]
            new_terms = []
            
            # iterates over all terms adjusting values
            for (attr, val) in terms:
                uniques = list(db[attr].unique())
                if val not in uniques:
                    convert_val = VAL_CONVERT[val]
                    new_val = [x for x in uniques if re.match(convert_val, x)][0]
                    new_terms.append((attr,new_val))
                else:
                    new_terms.append((attr,val))
            
            self.__ruleset.iloc[idx] = set(new_terms)
        return
    
    def set_cases(self, db):
        
        self.__adjust_values(db)
        
        def get_antecedent(set_items):
            dic = {}
            for (attr, v) in set_items:
                if attr in dic:
                    dic[attr].append(v)
                else:
                    dic[attr] = [v]
            return dic
        def get_queries(dic):
            queries = []
            for attr, values in dic.items():
                if len(values)==1:
                    queries.append(f'`{attr}` == {repr(values[0])}')
                else:
                    queries.append("(" + ' or '.join([f'`{attr}` == {repr(v)}' for v in values]) + ")")
            return ' and '.join(queries)     
        
        if db is None:
            self.__cases = pd.Series([])
        elif self.name in ALG_BASE.values(): # EsmamDS algorithm: treat disjunctions of terms
            antecedent = self.__ruleset.apply(lambda x: get_antecedent(x))
            queries = antecedent.apply(lambda x: get_queries(x))
            cases = queries.apply(lambda x: set(db.query(x).index))
            self.__cases = cases.rename('{}-cases'.format(self.name))
        else:
            queries = self.__ruleset.apply(lambda x: ' and '.join([f'`{k}` == {repr(v)}' for (k, v) in x]))
            cases = queries.apply(lambda x: set(db.query(x).index))
            self.__cases = cases.rename('{}-cases'.format(self.name))
        return
    
    @property
    def ruleset(self):
        return self.__ruleset
    
    @property
    def size(self):
        return self.__ruleset.shape[0]
    
    @property
    def cases(self):
        return self.__cases
    
    def intersection(self, ruleset):
        rules = []
        self_rules = self.__ruleset.to_list()
        
        for rule in range(ruleset.size):
            if ruleset.ruleset.iloc[rule] in self_rules:
                rules.append(ruleset.ruleset.iloc[rule])
        return pd.Series(rules)
    
    def union(self, ruleset):
        rules = self.__ruleset.to_list()
        
        for rule in range(ruleset.size):
            if ruleset.ruleset.iloc[rule] not in rules:
                rules.append(ruleset.ruleset.iloc[rule])
        return pd.Series(rules)
    
    def not_in(self, ruleset): # returns rules in self but not in ruleset
        rules = []
        other_rules = ruleset.ruleset.to_list()
        
        for rule in range(self.size):
            if self.__ruleset.iloc[rule] not in other_rules:
                rules.append(self.__ruleset.iloc[rule])
        return RuleSet(pd.Series(rules).rename('{}_not-in_{}'.format(self.name, ruleset.name)))


#### SCRIPT FUNCTIONS
def __read_data(db_name):

    data_file = DATA_PATH.format(db_name)
    json_file = data_file.split('.')[0] + '_dtypes.json'
    with open(json_file, 'r') as f:
        dtypes = json.load(f)
    return pd.read_csv(data_file, delimiter=',', header=0, index_col=False, compression='xz', dtype=dtypes)

def __load_model(db_name, alg, exp=None):
    def __get_terms(dic):
        terms = []
        for attr, val in dic.items():
            if isinstance(val, list):
                for v in val:
                    terms.append((attr, v))
            else:
                terms.append((attr,val))
        return set(terms)
    def __get_antecedent(series_val, alg):
        if alg in ['BS-EMM-pop','BS-SD-pop','BS-EMM-cpm','BS-SD-cpm']:
            list_of_items = [x.split('==') for x in series_val]
            return dict(list_of_items)
        else:
            list_of_items = [x.split(' = ') for x in series_val]
            return dict(list_of_items)

    if alg in ['BS-EMM-pop','BS-SD-pop','BS-EMM-cpm','BS-SD-cpm']:
        file = ALG_FILE[alg].format(db_name)
        db_results = pd.read_csv(file)
        model = db_results['subgroup'].apply(lambda x: x.replace("'","")).apply(lambda x: x.split(' AND ')).apply(lambda x: __get_antecedent(x, alg))
        model_terms = model.apply(lambda x: __get_terms(x))
        return model_terms.rename(alg)

    elif alg=='LR-Rules':
        file = ALG_FILE['LR-Rules'].format(db_name)
        db_results = pd.read_table(file, squeeze=True)
        final_idx = db_results[db_results == 'General information:'].index[0]
        db_results = db_results.iloc[:final_idx]
        model = db_results.apply(lambda x: x.split(' THEN ')[0]).apply(lambda x: x.split(': IF ')[1]).apply(lambda x: x.replace("{","")).apply(lambda x: x.replace("}","")).apply(lambda x: x.split(' AND ')).apply(lambda x: __get_antecedent(x, alg))
        model_terms = model.apply(lambda x: __get_terms(x))
        return model_terms.rename(alg)

    elif alg=='DSSD-CBSS':
        file = ALG_FILE['DSSD-CB'].format(db_name)
        db_results = pd.read_csv(file, sep='\t', index_col=0)
        model = db_results['Conditions'].apply(lambda x: x.replace("'","")).apply(lambda x: x.split(' AND ')).apply(lambda x: __get_antecedent(x, alg))
        model_terms = model.apply(lambda x: __get_terms(x))
        return model_terms.rename(alg)
    
    elif alg in ["Esmam-cpm","Esmam-pop","EsmamDS-cpm","EsmamDS-pop"]:
        if exp is None:
            raise ValueError('Exp number missing')
        
        file_name = ALG_FILE[alg].format(db_name, exp)
        
        if os.path.exists(file_name): # exists a .json file
            with open(file_name, 'r') as f:
                log = json.load(f)
        else: # dont exist a .json file >> read (.json).zip file
            zip_path = file_name.split('.')[0] + '.zip'
            with zipfile.ZipFile(zip_path) as z:
                for filename in z.namelist():
                    with z.open(filename) as f:  
                        data = f.read()
                        log = json.loads(data)
            
        model = pd.Series(map(lambda dic: dic['antecedent'], log['model'].values()))
        model_terms = model.apply(lambda x: __get_terms(x))    
        return model_terms.rename(alg)
    
    else:
        raise ValueError('Algorithm not found')

def __get_interset_matrix(my_model, other_model, db_name, index=None):
    '''
    Returns a table of either a jaccard index (for cover or description) or the logrank p-value for Rule-to-Rule comparison (rows vs columns)
    Row locations are ruleset (all rules) of the other model
    Column locations are all rules in the my model
    return: pd.DataFrame.to_dict()
    '''
    if index is None:
        raise ValueError('Define a metric for computing <__get_jaccard_matrix_rulesets>')
    
    def jaccard(set1, set2, small=True):
        union = set1.union(set2)
        if small:
            minor = min([len(set1),len(set2)])
            return len(intersection)/minor
        else:
            intersection = set1.intersection(set2)
            return len(intersection)/len(union) 
    
    def pval(cases1, cases2):
        time1 = db[VAR_SURV_TIME][cases1].to_numpy()
        status1 =  db[VAR_SURV_EVENT][cases1].to_numpy()
        time2 = db[VAR_SURV_TIME][cases2].to_numpy()
        status2 = db[VAR_SURV_EVENT][cases2].to_numpy()
        p_value = logrank_test(time1, time2, event_observed_A=status1, event_observed_B=status2).p_value
        return p_value
    
    # read dataset
    db = __read_data(db_name)
    
    # create RuleSet objects
    my_ruleset = RuleSet(my_model, db)
    other_ruleset = RuleSet(other_model, db)

    matrix = pd.DataFrame(data=None, index=list(my_ruleset.ruleset.index), columns=list(other_ruleset.ruleset.index))
    if matrix.shape[0]==0:
        return matrix.to_dict()
    
    ### CONSTRUCT MATRIX
    for row in my_ruleset.ruleset.index:
        for col in other_ruleset.ruleset.index:            
            if index=='cover':
                matrix.loc[row,col] = jaccard(other_ruleset.cases.loc[row], my_ruleset.cases.loc[col])
            elif index=='descr':
                matrix.loc[row,col] = jaccard(other_ruleset.ruleset.loc[row], my_ruleset.ruleset.loc[col])
            elif index=='model':
                matrix.loc[row,col] = pval(other_ruleset.cases.loc[row], my_ruleset.cases.loc[col])
    
    return matrix.astype(dict.fromkeys(matrix.columns, 'float64')).to_dict()


def generate_logs(metric, baseline):
    
    print('\n>> generating logs for baseline={}, metric={}'.format(baseline, metric))
    log_matrix = {}.fromkeys(DATASETS)
    
    base_alg = ALG_BASE[baseline]
    
    for db_name in DATASETS:
        print('..data set: {}'.format(db_name.upper()))
        log_matrix[db_name] = {}.fromkeys(ALGORITHMS[baseline])
        
        for alg in ALGORITHMS[baseline]:
            print('{} vs {}'.format(base_alg, alg))
                
            log_matrix[db_name][alg] = {}.fromkeys(list(range(RUNS)))
            
            for exp in range(RUNS):
                comp_model = __load_model(db_name, alg, exp)
                base_model = __load_model(db_name, base_alg, exp)
                log_matrix[db_name][alg][exp] = __get_interset_matrix(base_model, comp_model, db_name, metric)
                
    save_name = SAVE_PATH + SAVE_FILE.format(baseline, LOG_FILE[metric])
    with open(save_name, 'w') as f:
        json.dump(log_matrix, f)
    print('..saved: {}'.format(save_name))
    return

def run(metric='all', baseline=['population', 'complement']):
    
    for base in baseline:
        
        if metric in ['descr', 'all']:
            generate_logs('descr', base)
        if metric in ['cover', 'all']:
            generate_logs('cover', base)
        if metric in ['model', 'all']:
            generate_logs('model', base)
        
    return

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Script to generate the file logs for inter-set similarity analysis')
    
    parser.add_argument("--base", type=str,
                        choices=["all", "complement", "population"],
                        default="all",
                        help="Baseline for generating the results.")
    parser.add_argument("--metric", type=str,
                        choices=["all", "descr", "cover", "model"],
                        default="all",
                        help="Metric for generating the results.")
    args = parser.parse_args()
    
    if args.base == 'all':
        baselines = ['population', 'complement']
    else:
        baselines = [args.base]
    
    run(args.metric, baselines)