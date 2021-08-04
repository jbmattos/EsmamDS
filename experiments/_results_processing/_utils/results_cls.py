# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 10:30:23 2021

@author: jubma
"""

import json
import numpy as np
import pandas as pd
import pathlib


class Results():
    def __init__(self):
        self.ROOT = "EsmamDS"
        self.ROOT_PATH = str(pathlib.Path(__file__).parent.absolute()).split(self.ROOT)[0].replace('\\','/') + self.ROOT + '/'
        self.DATA_PATH = self.ROOT_PATH + 'data sets/final data sets/{}_disc.xz'
        self.PROC_PATH = self.ROOT_PATH + 'experiments/_results_processing/_processed_output_files/'
        self.SAVE_FILE = '{}_{}'
        
        self.ALGORITHMS = {'population': ['EsmamDS-pop', 'Esmam-pop', 'BS-EMM-pop', 'BS-SD-pop', 'DSSD-CBSS'],
                           'complement': ['Esmam-cpm', 'EsmamDS-cpm', 'BS-EMM-cpm', 'BS-SD-cpm','LR-Rules']}
        self.ESMAM_VARS = ["Esmam-cpm","Esmam-pop","EsmamDS-cpm","EsmamDS-pop"]
        self.ALG_FILE = {'Esmam-pop': 'esmam-pop/esmam-pop_{}',
                        'Esmam-cpm': 'esmam-cpm/esmam-cpm_{}',
                        'EsmamDS-pop': 'esmamds-pop/esmamds-pop_{}', 
                        'EsmamDS-cpm': 'esmamds-cpm/esmamds-cpm_{}',
                        'BS-EMM-pop': 'bs-emm-pop/pysg_{}',
                        'BS-EMM-cpm': 'bs-emm-cpm/pysg_{}',
                        'BS-SD-pop': 'bs-sd-pop/pysg_{}',
                        'BS-SD-cpm': 'bs-sd-cpm/pysg_{}',
                        'DSSD-CBSS': 'dssd-cbss/cortana_{}',
                        'LR-Rules': 'lr-rules/lrrules_{}'}
        self.LOG_FILE = {'metrics': 'metrics.csv',
                        'jaccard_c': 'jaccard-cover.json',
                        'jaccard_d': 'jaccard-descript.json',
                        'pval_r': 'rules-pval.json',
                        'pval_m': 'pval-matrix.json',
                        'surv_models': 'survival-models.json'
                       }
        self.METRICS = {'num_rules': '#sg',
                       'length': 'length',
                       'rule_coverage': 'sgCov',
                       'set_coverage': 'setCov',
                       'cr': 'CR',
                      }
        self.DATASETS = ['actg320','breast-cancer','cancer','carcinoma','gbsg2','lung','melanoma','mgus2','mgus','pbc','ptc','uis','veteran','whas500']
        self.RUNS = 30
        self.ALPHA = 0.05
        

class Table(Results):
    def __init__(self):
        super().__init__()
        self.SAVE_PATH = self.ROOT_PATH + 'experiments/metrics results (tables and statistics)//'
        self.__tbl_metrics = None
        self.__tbl_redundancy = None
        self.__tbl_exceptionality = None
        self.__final_tbl = None
        self.__baseline = None
    
    def __adjust_tbls(self):
        
        # adjust redundancy
        red_names = ['description redundancy', 'cover redundancy', 'CR', 'model redundancy']
        df_red = pd.concat([self.__tbl_redundancy['description redundancy'],
                            self.__tbl_redundancy['cover redundancy'],
                            self.__tbl_metrics['CR'],
                            self.__tbl_redundancy['model redundancy']], 
                           axis=1, keys=red_names, names=['Metrics', 'Algorithms'])
        self.__tbl_redundancy = df_red.copy()
        
        # adjust metrics
        met_names = ['#sg', 'length', 'sgCov', 'setCov']
        self.__tbl_metrics = self.__tbl_metrics[met_names].copy()
        return
    
    def __set_tbl_metrics(self, baseline):
        '''
        Generates a dataframe with the performance of all analysed algorithms in the following metrics:
            >> #sg, length, sgCov, setCov, CR            
        '''
        log_file = 'metrics'        
        def __read_log(alg):
            if alg in self.ESMAM_VARS:
                return pd.read_csv(self.PROC_PATH + self.ALG_FILE[alg].format(self.LOG_FILE[log_file])).set_index('Dataset')
            else:
                df = pd.read_csv(self.PROC_PATH + self.ALG_FILE[alg].format(self.LOG_FILE[log_file]), index_col=0)
                df = df.reindex(pd.Series(index_expand, name='Dataset'))
                return df
        
        index_expand = [db for db in self.DATASETS for i in range(self.RUNS)]
        df_metrics = {}
        
        for alg in self.ALGORITHMS[baseline]:
            df_metrics[alg] = __read_log(alg)

        df_unif = pd.concat([pd.concat([df[metric].rename(k) for k, df in df_metrics.items()], axis=1) for metric in self.METRICS.keys()],
                            axis=1, keys=self.METRICS.values(), names=['Metrics','Algorithms'])
        
        df_unif.index.name = 'Data sets'
        df_unif['Exp'] = list(range(30))*14
        df_unif.set_index('Exp', append=True, inplace=True)
        self.__tbl_metrics = df_unif.copy()
        if isinstance(self.__tbl_redundancy, pd.DataFrame):
            self.__adjust_tbls()
        return
    
    def __set_tbl_redundancy(self, baseline):
        '''
        RE WRITE
        Generates a dataframe with the average jaccard (cover/descript) and pvalue performance of each algorithm [models]
        > metrics: jaccard-cover, jaccard-descript, pval-matrix
        [for a proper csv-importing use >> pd.read_csv('metrics_table.csv', header=0, index_col=[0,1])]
        '''
        
        def __read_matrix(alg, dataset, metric, exp):
            with open(self.PROC_PATH + self.ALG_FILE[alg].format(self.LOG_FILE[metric]), 'r') as f:
                log = json.load(f)
            if alg in self.ESMAM_VARS:
                matrix = pd.DataFrame(log[dataset][str(exp)])
            else:
                matrix = pd.DataFrame(log[dataset])
            if metric == "pval_m":
                matrix = matrix > self.ALPHA
                matrix = matrix.astype(int)
            return matrix
        
        def __compute_metric(matrix): # normalized sum of lower triangular (w/o diagnonal)
            m_sum = pd.DataFrame(np.tril(matrix, k=-1)).sum().sum()
            idx = np.tril_indices(len(matrix.index), k=-1)
            return m_sum/len(list(zip(idx[0], idx[1])))
            
        dfs = {}.fromkeys(['jaccard_d','jaccard_c','pval_m'])
        
        for metric in ['jaccard_d','jaccard_c','pval_m']:
            
            df_metric = {}.fromkeys(self.DATASETS)
            for db in self.DATASETS:
                df_metric[db] = {}.fromkeys(self.ALGORITHMS[baseline]).copy()
                
                for alg in self.ALGORITHMS[baseline]:
                    df_metric[db][alg] = {}.fromkeys([str(i) for i in range(self.RUNS)]).copy()

                    for exp in range(self.RUNS):
                        matrix = __read_matrix(alg, db, metric, exp)
                        
                        df_metric[db][alg][str(exp)] = __compute_metric(matrix)
            
            dfs[metric] = pd.concat([pd.DataFrame(df_metric[db]) for db in df_metric.keys()], axis=0, keys=df_metric.keys(), names=['Datasets','Exp'])
        
        df_unif = pd.concat(dfs.values(), axis=1, keys=['description redundancy', 'cover redundancy', 'model redundancy'], names=['Metrics', 'Algorithms'])
        df_unif.index.names = ('Data sets', 'Exp')
        df_unif.index.set_levels(df_unif.index.levels[1].astype(int), level=1, inplace=True)
        self.__tbl_redundancy = df_unif.copy()
        if isinstance(self.__tbl_metrics, pd.DataFrame):
            self.__adjust_tbls()
        return
    
    def __set_tbl_exceptionality(self, baseline):
        
        res_except = {}.fromkeys(self.ALGORITHMS[baseline], None)
        res = {}.fromkeys(self.DATASETS, None)

        for db in self.DATASETS:

            for alg in self.ALGORITHMS[baseline]:

                res_except[alg] = {}.fromkeys(range(self.RUNS), None)

                file = self.PROC_PATH + self.ALG_FILE[alg].format('rules-pval.json')
                with open(file, 'r') as f:
                    log = json.load(f)

                result_db = None
                if 'Esmam' in alg:
                    for idx in range(self.RUNS):
                        df = pd.DataFrame(log[db][str(idx)]).T < self.ALPHA
                        result_db = df.sum().T / df.shape[0] # number of not exceptional by number of rules
                        res_except[alg][idx] = result_db['pval_logrank']

                else:
                    df = pd.DataFrame(log[db]).T < self.ALPHA
                    result_db = df.sum() / df.shape[0] # number of not exceptional by number of rules
                    for idx in range(self.RUNS):
                        res_except[alg][idx] = result_db['pval_logrank']

            res[db] = pd.DataFrame(res_except)
        df_unif = pd.concat([df for df in res.values()], axis=0, keys=res.keys(), names=['Data sets', 'Exp'])
        df_unif.columns = pd.MultiIndex.from_product([['exceptionality'], df_unif.columns])
        self.__tbl_exceptionality = df_unif.copy()
        return 
    
    def set_table(self, baseline):
        
        self.__baseline = baseline
        self.__set_tbl_metrics(baseline)
        self.__set_tbl_redundancy(baseline)
        self.__set_tbl_exceptionality(baseline)
        
        tbl = self.__tbl_exceptionality.join(self.__tbl_metrics, how='outer').join(self.__tbl_redundancy, how='outer')
        tbl.columns.names = ('Metrics', 'Algorithms')
        self.__final_tbl = tbl
        return
    
    def get_table(self):
        return self.__final_tbl
    
    def save(self):
        file_name = self.SAVE_PATH+'metrics_baseline-{}.csv'.format(self.__baseline)
        self.__final_tbl.to_csv(file_name)
        print(".saved: {}".format(file_name))
        return
    
    
if __name__ == '__main__':
    tbl = Table()
    print(tbl._LOG_PATH)
        
        
        
        
        