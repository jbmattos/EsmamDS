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
import seaborn as sns

from matplotlib import pyplot as plt
from matplotlib.gridspec import SubplotSpec
from matplotlib.pylab import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable


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
        rcParams["font.size"] = 40
        rows = 1
        cols = len(ALGS)
        rcParams['figure.figsize'] = 9*cols, 6*rows
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
            ax.set_yticks([0.5,1])
            ax.set_xlim([0, x[-1]+3])
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.legend().set_visible(False)
            plt.xticks([])
            plt.subplots_adjust(wspace=0.05)
        
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
        file_name = 'survivalModels_baseline-{}.pdf'.format(baseline)
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
        Analysis of the Survival Models associated with the subgroups 
        discovered by each compared algorithm.
         
        This function generates a report of the full empirical evaluation in a 
        pdf file where each page presents the results regarding one data set.
        In each page, the rows are the experiment runs, and the columns are the
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


class Heatmaps(Results):
    '''
    Class for generating the Heatmap results of the empirical evaluation
    '''
    def __init__(self):
        super().__init__()
        self.SAVE_PATH = {'sim': self.ROOT_PATH + 'experiments/sets similarities (heatmap matrix)/',
                          'red': self.ROOT_PATH + 'experiments/set redundancy (heatmap matrix)/'}
        self.__LOG_PATH = self.ROOT_PATH + 'experiments/_results_processing/_processed_output_files/_inter-set_similarity/'
        self.__color_map_singlerun = sns.diverging_palette(240, 10, n=100, l=30)
        self.__color_map = 'RdYlGn_r'
        self.__ALGORITHMS = {'population': ['Esmam-pop', 'BS-EMM-pop', 'BS-SD-pop', 'DSSD-CBSS'],
                              'complement': ['Esmam-cpm', 'BS-EMM-cpm', 'BS-SD-cpm','LR-Rules']
                             }
        rcParams["font.family"] = "Times New Roman"
 
    
    def __generate_single_interset(self, baseline, dic_matrix, max_x, max_y, db, exp, save, _metrics):
        '''
        Generate the pdf file that comprises the < single_run_interset_similarity >.
        
        Parameters
        ----------
        baseline : string
            options=['population', 'baseline']
            The baseline family of algorithms to generate the results.
        dic_matrix : dictionary
            Is the '{baseline}_interset_*Similarity.json' files.
        metric : string
            options=['jaccard_c', 'jaccard_d', 'pval_m']
            (Log file) Information considered to be plotted
        max_x : int
            Maximum value of the x-axis for all considered plots. 
            Used for adjustment on the aspect ratio of the plots.
        max_y : int
            Maximum value of the y-axis for all considered plots. 
            Used for adjustment on the aspect ratio of the plots.
        db : string
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
        
        # set the EsmamDS variant as baseline algorithm
        if baseline=='population': alg_base = 'EsmamDS-pop'
        else: alg_base = 'EsmamDS-cpm'
        
        # set figure
        rcParams["font.family"] = "Times New Roman"
        rcParams["font.size"] = 22
        rows = 1
        cols = len(self.__ALGORITHMS[baseline])
        rcParams['figure.figsize'] = 5*cols, 4*rows
        sns.set_palette("deep")
        
        
        for metric in _metrics:                                     # rows: iterates over metrics (descript, cover, model)
        
            fig, axes = plt.subplots(nrows=rows, ncols=cols, num=metric, clear=True, sharey=True, sharex=False)
            plt.subplots_adjust(wspace=0.3)
        
            for col, alg in enumerate(dic_matrix[metric].keys()):   # cols: iterates over algorithms results
                            
                matrix = dic_matrix[metric][alg]
                # adjust for boolean plot
                if metric=='modelSimilarity' and not matrix.shape[0]==0:      
                    matrix = matrix >= self.ALPHA
                    matrix = matrix.astype(int)
                # square-ratio: insert Nan rows at the right of dataframe for matrix with less that max_y
                # if matrix.shape[1]<max_y:
                #     delta = max_y - matrix.shape[1]
                #     delta_m = pd.DataFrame(data=np.nan, index=matrix.index, columns=range(delta))
                #     matrix = pd.concat([matrix,delta_m], axis=1).reset_index(drop = True)
                
                ### plot
                ax = axes[col]
                with sns.axes_style("white"):
                    ax = sns.heatmap(matrix, cmap=self.__color_map, ax=ax, vmin=0, vmax=1, cbar=False)
                # insert colorbar on the right plot
                if col==cols-1: 
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes('right', size='5%', pad=0.3)
                    plt.colorbar(ax.collections[0], cax=cax)
                    
                ax.set_yticks([])
                ax.set_xticks([])
                # title
                if metric=='descrSimilarity': ax.set_title(alg, {'fontweight':'bold'})
                # x-labels (on bottom plots)
                if metric=='modelSimilarity' and col==0: ax.set_xlabel('Subgroups')
                else: ax.set_xlabel('')
                # y-labels (on left plots)
                if col==0: ax.set_ylabel('{}'.format(alg_base), {'fontweight':'bold'})
                else: ax.set_ylabel('')
                
            if save:
                save_name = self.SAVE_PATH['sim'] +'_intersetSimilarity_baseline-{}_{}-exp{}_{}.pdf'.format(baseline, db, exp, metric)
                plt.savefig(save_name, bbox_inches='tight') 
                print('..saved: {}'.format(save_name))
            else:
                plt.show()
        return
    
    def __generate_single_intraset(self, dic_matrix, metric, baseline, db, exp, save, boolean=False):
        '''
        Generate the results' file that comprises the < single_run_intraset_redundancy >.
        
        Parameters
        ----------
        dic_matrix : dictionary
            For each compared algorithm(keys), the dictionary value is its 
            respective < _results_processing/_processed_output_files > log file
            related to the metric attribute
        metric : string
            options=['jaccard_c', 'jaccard_d', 'pval_m']
            (Log file) Information considered to be plotted
        baseline : string
            options=['population', 'baseline']
            The baseline family of algorithms to generate the results.
        db : string
            Name of the data set related to the results to be plotted
        exp : int
            Index of the data set experiment (run) to plot the results
        save : bool, optional
            Whether to save the plot as pdf file or not (just show). 
            The default is True.
        boolean : bool, optional
            Whether to consider a boolean plot (transform a float matrix into 
                                                boolean according to the
                                                self.ALPHA threshold).
            The default is False.

        Returns
        -------
        None.
        '''

        # figure settings
        plt.rcParams.update({'font.size': 22})
        rows = 1
        cols = len(self.ALGORITHMS[baseline])     # number of algorithms
        rcParams['figure.figsize'] = 5*cols, 4*rows
        
        fig, axes = plt.subplots(nrows=rows, ncols=cols, num='plt_{}'.format(db), clear=True)
        plt.subplots_adjust(wspace=0.005)
        
        for col, alg in enumerate(self.ALGORITHMS[baseline]):    # iterates over subplots COLUMNS
            
            if alg in self.ESMAM_VARS:
                matrix = pd.DataFrame(dic_matrix[alg][db][str(exp)])
            else:
                matrix = pd.DataFrame(dic_matrix[alg][db])
            
            if boolean:
                matrix = matrix >= self.ALPHA
                matrix = matrix.astype(int)
            
            plt.figure(num='plt_{}'.format(db))
            ax = axes[col]  
            
            with sns.axes_style("white"):
                mask = np.zeros_like(matrix)
                mask[np.triu_indices_from(matrix)] = True
                ax = sns.heatmap(matrix, mask=mask, square=True, cmap=self.__color_map, ax=ax, vmin=0, vmax=1, cbar=False)
            if col==cols-1: 
                divider = make_axes_locatable(ax)
                cax = divider.append_axes('right', size='5%', pad=0.05)
                plt.colorbar(ax.collections[0], ax=ax, cax=cax)
            
            ax.set_yticks([])
            ax.set_xticks([])
            if metric=='modelRedundancy' and col==0: 
                ax.set_ylabel('Subgroups')
                ax.set_xlabel('Subgroups')
            else: 
                ax.set_ylabel('')
                ax.set_xlabel('')
            if metric=='descrRedundancy': ax.set_title(alg, {'fontweight':'bold'})
        
        if save:
            save_name = self.SAVE_PATH['red'] +'_intrasetRedundancy_baseline-{}_{}-exp{}_{}.pdf'.format(baseline, db, exp, metric)
            plt.savefig(save_name, bbox_inches='tight') 
            print('..saved: {}'.format(save_name))
        else:
            plt.show()
        return
    
    def __generate_redundancy_heatmap(self, dic_matrix, ALGS, metric, baseline, boolean=False):
        '''
        Generate the pdf file that comprises the < redundancy_full_report >.
        
        Parameters
        ----------
        dic_matrix : dictionary
            For each compared algorithm(keys), the dictionary value is its 
            respective < _results_processing/_processed_output_files > log file
            related to the metric attribute
        metric : string
            options=['jaccard_c', 'jaccard_d', 'pval_m']
            (Log file) Information considered to be plotted
        baseline : string
            options=['population', 'baseline']
            The baseline family of algorithms to generate the results.
        boolean : bool, optional
            Whether to consider a boolean plot (transform a float matrix into 
                                                boolean according to the
                                                self.ALPHA threshold).
            The default is False.

        Returns
        -------
        None.
        '''

        # figure settings
        plt.rcParams.update({'font.size': 10})
        rows = self.RUNS    # number of experiments
        cols = len(ALGS)    # number of algorithms
        rcParams['figure.figsize'] = 5*cols, 4*rows
        
        # pdf file
        _m_dic = {'jaccard_d': 'descrRedundacy', 
                  'jaccard_c': 'coverRedundacy', 
                  'pval_m': 'modelRedundacy'}
        save_name = self.SAVE_PATH['red'] + '{}_{}.pdf'.format(baseline, _m_dic[metric])
        pdf = matplotlib.backends.backend_pdf.PdfPages(save_name)
        
        for db in self.DATASETS:
            print('...generating page results for {}'.format(db))
            
            fig, axes = plt.subplots(nrows=rows, ncols=cols, num='plt_{}'.format(db), clear=True)
            
            for col, alg in enumerate(ALGS):    # iterates over subplots COLUMNS
                
                for row in range(self.RUNS):         # iterates over subplots ROWS
                    
                    if alg in self.ESMAM_VARS:
                        matrix = pd.DataFrame(dic_matrix[alg][db][str(row)])
                    else:
                        matrix = pd.DataFrame(dic_matrix[alg][db])
                    
                    if boolean:
                        matrix = matrix >= self.ALPHA
                        matrix = matrix.astype(int)
                    
                    # jaccard plot
                    plt.figure(num='plt_{}'.format(db))
                    ax = axes[row, col]
                    
                    if col==cols-1: cbar = True
                    else: cbar = False
                    
                    with sns.axes_style("white"):
                        mask = np.zeros_like(matrix)
                        mask[np.triu_indices_from(matrix)] = True
                        ax = sns.heatmap(matrix, mask=mask, square=True, cmap=self.__color_map, ax=ax, vmin=0, vmax=1, cbar=cbar)
                    ax.title.set_text(alg)
                    ax.set_yticks([])
                    ax.set_xticks([])
            
            # page titles
            fig.suptitle('{} DATASET (Experiments x Algorithms)'.format(db.upper()), fontsize=22) 
            # saving on pdf
            pdf.savefig(fig,bbox_inches='tight')
            plt.close(fig)
            
        pdf.close()
        print('..saved: {}'.format(save_name))
        return
       
    def __generate_similarity_heatmap(self, dic_matrix, metric, baseline, boolean=False):
        '''
        Generate the pdf file that comprises the < similarity_full_report >.
        
        Parameters
        ----------
        dic_matrix : dictionary
            Is the '{baseline}_interset_*Similarity.json' files.
        metric : string
            options=['descrSimilarity', 'coverSimilarity', 'modelSimilarity']
            The type of similarity measures to compute the results.
        baseline : string
            options=['population', 'baseline']
            The baseline family of algorithms to generate the results.
        boolean : bool, optional
            Whether to consider a boolean plot (transform a float matrix into 
                                                boolean according to the
                                                self.ALPHA threshold).
            The default is False.

        Returns
        -------
        None.

        '''
        
        save_name = self.SAVE_PATH['sim'] + '{}_{}.pdf'.format(baseline, metric)
        pdf_jaccard = matplotlib.backends.backend_pdf.PdfPages(save_name)
        
        # figure settings
        plt.rcParams.update({'font.size': 10})
        rows = self.RUNS                            # number of experiments
        cols = len(self.__ALGORITHMS[baseline])     # number of algorithms
        rcParams['figure.figsize'] = 5*cols, 4*rows
        
        for db in self.DATASETS:
            print('...generating page results for {}'.format(db))
            
            fig_jaccard, axes_jaccard = plt.subplots(nrows=rows, ncols=cols, num='jaccard_{}'.format(db), clear=True)
            
            for col, alg in enumerate(self.__ALGORITHMS[baseline]):    # iterates over subplots COLUMNS
                
                for row in range(self.RUNS):                         # iterates over subplots ROWS
                    
                    matrix = pd.DataFrame(dic_matrix[db][alg][str(row)])
                    if boolean and not matrix.shape[0]==0:
                        matrix = matrix >= self.ALPHA
                        matrix = matrix.astype(int)
                    if matrix.shape[0]==0:
                        matrix = pd.DataFrame(data=np.nan, index=matrix.columns, columns=matrix.columns)
                    
                    # jaccard plot
                    plt.figure(num='jaccard_{}'.format(db))
                    ax = axes_jaccard[row, col]
                    
                    if col==cols-1: cbar = True
                    else: cbar = False
                    
                    with sns.axes_style("white"):
                        ax = sns.heatmap(matrix, square=True, cmap=self.__color_map, ax=ax, vmin=0, vmax=1, cbar=cbar)
                    
                    if baseline=='population': alg_base = 'EsmamDS-pop'
                    else: alg_base = 'EsmamDS-cpm'
                    if col==0: ax.set_ylabel('{}'.format(alg_base))
                    else: ax.set_ylabel('')
                    ax.set_xlabel(alg)
                    ax.set_yticks([])
                    ax.set_xticks([])
                    
            # page titles
            fig_jaccard.suptitle('{} DATASET (Experiments x Algorithms)'.format(db.upper()), fontsize=22) 
            # saving on pdf
            pdf_jaccard.savefig(fig_jaccard,bbox_inches='tight')
            plt.close(fig_jaccard)
            
        pdf_jaccard.close()
        print('..saved: {}'.format(save_name))
        return
    
    def similarity_full_report(self, baseline):
        '''
        Analysis of the similarity between the EsmamDS sets of subgroups and 
        the set of subgroups discovered by the compared algorithm.
        
        This function generates a report of the full empirical evaluation in 
        three pdf files, one for each similarity measure: descrition, coverage
        and model.
        In the pdf files, each page presents the results regarding one data set.
        In each page, the rows are the experiment runs, and the columns are the
        comparison of EsmamDS with the different algorithms.
        
        In the heatmaps, each row represents a subgroup discovered by the EsmamDS,
        and each column is a subgroup discovered by the compared algorithm.
        
        Parameters
        ----------
        baseline : string
            options=[population, baseline]
            The baseline family of algorithms to generate the results

        Returns
        -------
        None.

        '''
        print('> Call to Heatmaps().similarity_full_report(baseline={})'.format(baseline))
        metrics = ['descrSimilarity', 'coverSimilarity', 'modelSimilarity']
        
        for metric in metrics:
            print(".processing {} result pdf file".format(metric))
            boolean = False                
            if metric=='modelSimilarity': boolean = True
            
            file = self.__LOG_PATH +'{}_interset_{}.json'.format(baseline, metric)
            with open(file, 'r') as f:
                log_dic = json.load(f)             
            
            self.__generate_similarity_heatmap(log_dic, metric, baseline, boolean=boolean)
        return 

    def redundancy_full_report(self, baseline):
        '''
        Analysis of the similarity between the subgroups comprising the final 
        sets of subgroups discovered by each compared approach. 
        In other words, analysis of (set) redundancy.
        
        This function generates a report of the full empirical evaluation in a 
        pdf file where each page presents the results regarding one data set.
        In each page, the rows are the experiment runs, and the columns are the
        results of different algorithms.
        
        
        In the heatmaps, the rows and columns are the subgroups in a set, 
        and the values plotted (in color) are the similarity between the pairs
        of subgroups.

        Parameters
        ----------
        baseline : string
            options=[population, baseline]
            The baseline family of algorithms to generate the results

        Returns
        -------
        None.

        '''
        print('> Call to Heatmaps().redundancy_full_report(baseline={})'.format(baseline))
        metrics = ['jaccard_c', 'jaccard_d', 'pval_m']
        _metric_dic = {'jaccard_d': 'descrRedundancy',
                       'jaccard_c': 'coverRedundancy',
                       'pval_m': 'modelRedundancy'}
        
        def __read_logs(alg, log_file):
            with open(self.PROC_PATH + self.ALG_FILE[alg].format(self.LOG_FILE[log_file]), 'r') as f:
                log = json.load(f)
            return log
        
        # read matrixes for all algorithms
        algs = self.ALGORITHMS[baseline]
        dic_matrix = {}.fromkeys(algs)
        for metric in metrics:
            print(".processing {} result pdf file".format(_metric_dic[metric]))
            
            for alg in algs:
                dic_matrix[alg] = __read_logs(alg, metric)
            
            if metric == 'pval_m':
                self.__generate_redundancy_heatmap(dic_matrix, algs, metric, baseline, boolean=True)
            else:
                self.__generate_redundancy_heatmap(dic_matrix, algs, metric, baseline, boolean=False)
        
        return

    def single_run_interset_similarity(self, baseline, db_name, exp, save=True):
        '''
        Generates a plot for inter-set (description, coverage and model) similarity
        between the EsmamDS sets of subgroups and the set of subgroups discovered 
        by the compared algorithm.
        
        Each type of similarity between pairs of subgroups is displayed in a row, 
        where the plots in a row are the comparison between EsmamDS and the different algorithms.

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
        print('> Call to Heatmaps().single_run_interset_similarity(baseline={}, db_name={}, exp={}, save={})'.format(baseline, db_name, exp, save))
        _metrics = ['descrSimilarity', 'coverSimilarity', 'modelSimilarity']
        
        # read logs for db_name/exp
        log_dic = {}.fromkeys(_metrics)
        max_x, max_y = 0,0
        
        for metric in _metrics:  
    
            file = self.__LOG_PATH +'{}_interset_{}.json'.format(baseline, metric)
            with open(file, 'r') as f:
                log = json.load(f)
            
            log_dic[metric] = {}.fromkeys(log[db_name].keys())
            for alg in log_dic[metric].keys():
                log_dic[metric][alg] = pd.DataFrame(log[db_name][alg][str(exp)])
                if log_dic[metric][alg].shape[0] > max_x: max_x = log_dic[metric][alg].shape[0]
                if log_dic[metric][alg].shape[1] > max_y: max_y = log_dic[metric][alg].shape[1]
        
        # generate plots
        self.__generate_single_interset(baseline, log_dic, max_x, max_y, db_name, exp, save, _metrics)
        return 

    def single_run_intraset_redundancy(self, baseline, db_name, exp, save=True):
        '''
        Generates a plot for intra-set (description, coverage and model) redundancy
        between the subgroups discovered by the EsmamDS and all the compared approaches.
        
        The plots are triangular matrixes < subgroup x subgroup >, 
        i.e. both x and y-axis are the discovered subgroups.
        The plots display the metrics of similarity between pairs of subgroups
        for each subgroup compared to all others in the set.
        
        Each similarity metric is plotted in a different pdf file. 
        The 3 pdf files (for description, cover and model similarities) are 
        formatted so they are displayed vertically in such an order. 

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
        print('> Call to Heatmaps().single_run_intraset_redundancy(baseline={}, db_name={}, exp={}, save={})'.format(baseline, db_name, exp, save))
        metrics = ['jaccard_c', 'jaccard_d', 'pval_m']
        _metric_dic = {'jaccard_d': 'descrRedundancy',
                       'jaccard_c': 'coverRedundancy',
                       'pval_m': 'modelRedundancy'}
        
        def __read_logs(alg, log_file):
            with open(self.PROC_PATH + self.ALG_FILE[alg].format(self.LOG_FILE[log_file]), 'r') as f:
                log = json.load(f)
            return log
        
        # dictionary of matrix (for all algorithms)
        algs = self.ALGORITHMS[baseline]
        dic_matrix = {}.fromkeys(algs)
        
        for metric in metrics:
            
            # read matrixes for all algorithms
            for alg in algs:
                dic_matrix[alg] = __read_logs(alg, metric)
            
            # call function to plot 
            if metric == 'pval_m':
                self.__generate_single_intraset(dic_matrix, _metric_dic[metric], baseline, db_name, exp, save, boolean=True)
            else:
                self.__generate_single_intraset(dic_matrix, _metric_dic[metric], baseline, db_name, exp, save, boolean=False)
        
        return

if __name__ == '__main__':
    
    #m = SurvivalPlots()
    #m.single_run_plot('population', 'pbc', 0)
    
    #########################################################
    ## Plot inter-set similarities for single run comparison
    #########################################################
    #
    m = Heatmaps()
    m.single_run_interset_similarity('population', 'pbc', 0)
    
    #########################################################
    ## Plot intra-set redundancy for single run comparison
    #########################################################
    #
    m = Heatmaps()
    m.single_run_intraset_redundancy('population', 'pbc', 0)