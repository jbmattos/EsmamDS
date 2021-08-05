"""
EMPIRICAL EVALUATION RESULTS CLASS

    This file contains the classes for generating the final results of 
    EsmamDS empirical evaluation
"""

import itertools
import json
import matplotlib.backends.backend_pdf
import numpy as np
import pandas as pd
import pathlib

from matplotlib import pyplot as plt
from matplotlib.pylab import rcParams



class Results():
    '''
    Parent class for defining the general structure of the repository 
    and the empirical evaluation procedure
    '''
    def __init__(self):
        self.ROOT = "EsmamDS"
        self.ROOT_PATH = str(pathlib.Path(__file__).parent.absolute()).split(self.ROOT)[0].replace('\\','/') + self.ROOT + '/'
        self.DATA_PATH = self.ROOT_PATH + 'data sets/final data sets/{}_disc.xz'
        self.PROC_PATH = self.ROOT_PATH + 'experiments/_results_processing/_processed_output_files/'
        self.SAVE_FILE = '{}_{}'
        
        self.ALGORITHMS = {'population': ['EsmamDS-pop', 'Esmam-pop', 'BS-EMM-pop', 'BS-SD-pop', 'DSSD-CBSS'],
                           'complement': ['EsmamDS-cpm', 'Esmam-cpm', 'BS-EMM-cpm', 'BS-SD-cpm', 'LR-Rules']}
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
    '''
    Class for generating the Table results of empirical evaluation
    '''
    def __init__(self):
        super().__init__()
        self.SAVE_PATH = self.ROOT_PATH + 'experiments/metrics results (tables and statistics)/'
        self.__tbl_metrics = None
        self.__tbl_redundancy = None
        self.__tbl_exceptionality = None
        self.__final_tbl = None
        self.__baseline = None
    
    def __adjust_tbls(self):
        '''
        Adjusts the metrics and redundancy tables, realocating the CR metric
        from the first to the second.

        Returns
        -------
        None.

        '''
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
    
    def __set_tbl_metrics(self):        
        '''
        Generates a dataframe with the performance of all analysed algorithms in the following metrics:
            - #sg, length, sgCov, setCov, CR  

        Returns
        -------
        None.

        '''
        
        baseline = self.__baseline
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
    
    def __set_tbl_redundancy(self):
        '''
        Generates a dataframe with the performance of all analysed algorithms in the following metrics:
            - description redundancy, coverage redundancy and model redundancy

        Returns
        -------
        None.

        '''
        
        baseline = self.__baseline
        
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
    
    def __set_tbl_exceptionality(self):
        '''
        Generates a dataframe with the performance of all analysed algorithms in the following metrics:
            - exceptionality

        Returns
        -------
        None.

        '''
        
        baseline = self.__baseline
        
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
        '''
        Generates the results' dataframe for the empirical evaluation performed 
        over 14 data sets (with 30 runs each).
        The dataframe contains the performance of all analysed algorithms in the following metrics:
            - exceptionality
            - #sg
            - length
            - sgCov
            - setCov  
            - description redundancy
            - coverage redundancy
            - CR 
            - model redundancy

        Parameters
        ----------
        baseline : string
            options=[population, baseline]
            The baseline family of algorithms to generate the results

        Returns
        -------
        None.

        '''
        print('> Call to Table().set_table()')
        
        self.__baseline = baseline
        self.__set_tbl_metrics()
        self.__set_tbl_redundancy()
        self.__set_tbl_exceptionality()
        
        tbl = self.__tbl_exceptionality.join(self.__tbl_metrics, how='outer').join(self.__tbl_redundancy, how='outer')
        tbl.columns.names = ('Metrics', 'Algorithms')
        self.__final_tbl = tbl
        return
    
    def get_table(self):
        '''
        Access to the final results' dataframe.

        Returns
        -------
        pandas.DataFrame

        '''
        return self.__final_tbl
    
    def save(self):
        '''
        To save the final results' dataframe in its specif location inside 
        EsmamDS repository.

        Returns
        -------
        None.

        '''
        
        file_name = self.SAVE_PATH+'metrics_baseline-{}.csv'.format(self.__baseline)
        self.__final_tbl.to_csv(file_name)
        print(".saved: {}\n".format(file_name))
        return


class SurvivalPlots(Results):
    '''
    Class for generating the Survival Plots results of the empirical evaluation
    '''
    def __init__(self):
        super().__init__()
        self.SAVE_PATH = self.ROOT_PATH + 'experiments/survival models/'
        self.COLORS = ['#67001f','#b2182b','#d6604d','#92c5de','#4393c3','#2166ac']
        rcParams["font.family"] = "Times New Roman"
    
    
    def __plot_single_run(self, db, exp, dic_matrix, baseline, ALGS, save):
        '''
        This function is responsible for generating the plot result 
        of a single run.

        Parameters
        ----------
        db : string
            Name of the data set related to the results to be plotted
        exp : int
            Index of the data set experiment (run) to plot the results
        dic_matrix : dictionary {alg: info}
            For each algorithm, it contains the info to generate a dataframe of 
            the KM estimates to be plotted.
        baseline : string
            options=[population, baseline]
            The baseline family of algorithms to generate the results
        ALGS : list
            List of compared algorithms (for a given baseline).
        save : bool
            Whether to save the plot as pdf file or not (just show).

        Returns
        -------
        None.

        '''
        
        # set figure
        rcParams["font.size"] = 32
        rows = 1
        cols = len(ALGS)
        rcParams['figure.figsize'] = 10*cols, 6*rows
        fig, axes = plt.subplots(nrows=rows, ncols=cols, num=db, clear=True, sharex=True, sharey=True)
        
    
        for col_id, alg in enumerate(ALGS):
            
            color_it = itertools.cycle(self.COLORS)
    
            #load km-estimates file
            if alg in self.ESMAM_VARS:
                kmModels = pd.DataFrame(dic_matrix[alg][db][str(exp)])
            else:
                kmModels = pd.DataFrame(dic_matrix[alg][db])
            
            x = kmModels.times.values
            columns = list(kmModels.columns)
            columns.remove('times')
            
            ax = axes[col_id]
            for column in columns:
                if column == "population":
                    ax.plot(x, kmModels[column], color='k', label='{}'.format(column), linestyle=':')
                else:
                    ax.plot(x, kmModels[column], label='{}'.format(column), color=next(color_it))
    
            ax.set_title('{}'.format(alg))
            if col_id==0:
                ax.set_xlabel('Time')
                ax.set_ylabel('Survival probability')
            else:
                ax.set_xlabel('')
                ax.set_ylabel('')
            ax.set_ylim([0, 1])
            ax.set_xlim([0, x[-1]+3])
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.legend().set_visible(False)
            plt.xticks([])
        
        if save:
            save_name = self.SAVE_PATH +'survivalModels_baseline-{}_{}-exp{}.pdf'.format(baseline, db, exp)
            plt.savefig(save_name, bbox_inches='tight')  
            print('..saved: {}'.format(save_name))
        else: 
            plt.show()
        return
    
    def __plot_full_report(self, dic_matrix, baseline, ALGS):
        '''
        This function is responsible for generating the pdf file of the 
        full report.

        Parameters
        ----------
        dic_matrix : dictionary {alg: info}
            For each algorithm, it contains the info to generate a dataframe of 
            the KM estimates to be plotted.
        baseline : string
            options=[population, baseline]
            The baseline family of algorithms to generate the results
        ALGS : list
            List of compared algorithms (for a given baseline).

        Returns
        -------
        None.

        '''
          
        # set figure
        rcParams.update({'font.size': 16})
        rows = self.RUNS
        cols = len(ALGS)
        rcParams['figure.figsize'] = 12*cols, 8*rows
        
        # initialize pdf
        file_name = 'TsurvivalModels_baseline-{}.pdf'.format(baseline)
        save_name = self.SAVE_PATH + file_name
        pdf = matplotlib.backends.backend_pdf.PdfPages(save_name)
        
        for db in self.DATASETS:
            print('...generating page results for {}'.format(db))
            fig, axes = plt.subplots(nrows=rows, ncols=cols, num=db, clear=True, sharex=True, sharey=True)
                        
            for col_id, alg in enumerate(ALGS):
    
                for exp in range(self.RUNS):
                    
                    #load km-estimates file
                    if alg in self.ESMAM_VARS:
                        kmModels = pd.DataFrame(dic_matrix[alg][db][str(exp)])
                    else:
                        kmModels = pd.DataFrame(dic_matrix[alg][db])
                    
                    color_it = itertools.cycle(self.COLORS)
                    x = kmModels.times.values
                    columns = list(kmModels.columns)
                    columns.remove('times')
                    
                    ax = axes[exp, col_id]
                    for column in columns:
                        if column == "population":
                            ax.plot(x, kmModels[column], color='k', label='{}'.format(column), linestyle=':')
                        else:
                            ax.plot(x, kmModels[column], label='{}'.format(column), color=next(color_it))
                    
                    ax.set_title('{}'.format(alg))
                    if col_id==0:
                        ax.set_xlabel('Time')
                        ax.set_ylabel('Survival probability')
                    else:
                        ax.set_xlabel('')
                        ax.set_ylabel('')
                    ax.set_ylim([0, 1])
                    ax.set_xlim([0, x[-1]+3])
                    ax.spines['top'].set_visible(False)
                    ax.spines['right'].set_visible(False)
                    ax.legend().set_visible(False)
                    plt.xticks([])
                
            fig.suptitle('{} Discovered Subgroups Survival Models'.format(db.upper()), fontsize=32)
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
        pdf.close()
        print('..saved: {}'.format(save_name))
        return
        
    def full_report(self, baseline):
        '''
        Generates the full empirical evaluation report of the Survival Models 
        associated with the subgroups discovered by each compared algorithm.
        This function generates a pdf file where each page presents the results
        regarding one data set.
        In the page, the rows are the experiment runs, and the columns are the
        results of different algorithms.

        Parameters
        ----------
        baseline : string
            options=[population, baseline]
            The baseline family of algorithms to generate the results

        Returns
        -------
        None.

        '''
        
        print('> Call to Survival().full_report()')
        metric = 'surv_models'
        
        def __read_logs(alg, log_file):
            with open(self.PROC_PATH + self.ALG_FILE[alg].format(self.LOG_FILE[log_file]), 'r') as f:
                log = json.load(f)
            return log
        
        algs = self.ALGORITHMS[baseline]
        dic_matrix = {}.fromkeys(algs)
        for alg in algs:
            dic_matrix[alg] = __read_logs(alg, metric)
        
        self.__plot_full_report(dic_matrix, baseline, algs)
        return
    
    def single_run_plot(self, baseline, db_name, exp,  save=True):
        '''
        Generates a plot with the survival model of the discovered subgroups 
        delivered by each compared algorithm.

        Parameters
        ----------
        baseline : string
            options=[population, baseline]
            The baseline family of algorithms to generate the results
        db_name : string
            Name of the data set related to the results to be plotted
        exp : int
            Index of the data set experiment (run) to plot the results
        save : bool, optional
            Whether to save the plot as pdf file or not (just show). 
            The default is True.

        Returns
        -------
        None.

        '''
        print('> Call to Survival().single_run_plot(baseline={}, db_name={}, exp={}, save={})'.format(baseline, db_name, exp, save))
        metric = 'surv_models'
        
        def __read_logs(alg, log_file):
            with open(self.PROC_PATH + self.ALG_FILE[alg].format(self.LOG_FILE[log_file]), 'r') as f:
                log = json.load(f)
            return log
        
        # read matrixes for all algorithms
        dic_matrix = {}.fromkeys(self.ALGORITHMS[baseline])
        for alg in self.ALGORITHMS[baseline]:
            dic_matrix[alg] = __read_logs(alg, metric)
                
        self.__plot_single_run(db_name, exp, dic_matrix, baseline, self.ALGORITHMS[baseline], save)
        
        return  