# -*- coding: utf-8 -*-
"""
THIS SCRIPT GENERATES THE PROCESSED LOG-FILES FOR ALL ALGORITHMS COMPARED IN 
THE EMPIRICAL EVALUATION

It runs the algorithm-specific processing script from the _utils folder
"""

import argparse
import errno
import os
import pathlib

from datetime import datetime
from _utils import process_cortana_results as dssd_cbss
from _utils import process_esmam_results as esmam_class
from _utils import process_lrrules_results as lr_rules
from _utils import process_pysg_results as bs_class


'''
ADJUST PATHS
'''
ROOT = "EsmamDS"
ROOT_PATH = str(pathlib.Path(__file__).parent.absolute()).split(ROOT)[0] + ROOT + '/'
SAVE_FOLDER = ROOT_PATH + 'experiments/_results_processing/_processed_output_files/'
SAVE_PATH = ROOT_PATH + 'experiments/_results_processing/_processed_output_files/{}/'
FOLDER_PATH = ROOT_PATH + 'experiments/_results_processing/_algorithms_output_files/{}/'

DB_NAMES = ['actg320','breast-cancer','cancer','carcinoma','gbsg2','lung','melanoma','mgus2','mgus','pbc','ptc','uis','veteran','whas500']
ALGORITHMS = ['esmamds-pop', 'esmamds-cpm', 'esmam-pop', 'esmam-cpm',
              'bs-emm-pop', 'bs-emm-cpm', 'bs-sd-pop', 'bs-sd-cpm',
              'lr-rules', 'dssd-cbss']


def __get_save_path(alg):
    save_path = SAVE_PATH.format(alg)
    # creates directory for saving results and logs
    if not os.path.exists(os.path.dirname(save_path)):
        try:
            os.makedirs(os.path.dirname(save_path))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    return save_path

def __write_info(file, _direc, _base):
    
    with open(file, 'w') as f:
        f.write("This file was automatically created by the execution of")
        f.write("\n< {}\{} > {}".format(_direc,_base, datetime.now().strftime('%Y%m%d %H:%M:%S')))
        f.write("\n.. path to read the output files: ")
        f.write("\n< {} >".format(FOLDER_PATH))
        
        f.write("\n\nFile's contents description:")
        f.write('\n>> 1. file: <algorithm-baseline>_metrics.csv')
        f.write('\n   A table <dataset(14x30) x metrics> with the performance of the specific <algorithm-baseline> on the N-experiments.')
        f.write('\n>> 2. file: <algorithm-baseline>_jaccard-cover.json')
        f.write('\n   A nested dictionary containing the jaccard-matrix [based on rule coverage] (Rules x Rules) of each algorithm final model.')
        f.write('\n   The jaccard-matrix is saved as pd.DataFrame.to_dict() and should be read back into pandas obj.')
        f.write('\n   .. file structure: {db_name: exp_N: jaccard_matrix}')
        f.write('\n>> 3. file: <algorithm-baseline>_jaccard-descript.json')
        f.write('\n   A nested dictionary containing the jaccard-matrix [based on rule description] (Rules x Rules) of each algorithm final model.')
        f.write('\n   The jaccard-matrix is saved as pd.DataFrame.to_dict() and should be read back into pandas obj.')
        f.write('\n   .. file structure: {db_name: exp_N: jaccard_matrix}')
        f.write('\n>> 4. file: <algorithm-baseline>_rules-pval.json')
        f.write('\n   A nested dictionary containing the logrank p-value [in comparison to population] of each rule in the final rule-models')
        f.write('\n   .. file structure: {db_name: exp_N: {Rule_ID: p-value}}')
        f.write('\n>> 5. file: <algorithm-baseline>_pval-matrix.json')
        f.write('\n   A nested dictionary containing the logrank p-value matrix (Rules x Rules) of each combination of rules in the final rule-models')
        f.write('\n   The pval-matrix is saved as pd.DataFrame.to_dict() and should be read back into pandas obj.')
        f.write('\n   .. file structure: {db_name: exp_N: pvalue_matrix}')
        f.write('\n>> 6. file: <algorithm-baseline>_survival-models.json')
        f.write('\n   A nested dictionary containing a table with the survival models for population and each rule in the final rule-models')
        f.write('\n   The survival table is saved as pd.DataFrame.to_dict() and should be read back into pandas obj.')
        f.write('\n   .. file structure: {db_name: exp_N: survival-table}')
    return

def __append_info(file, algs, _direc, _base):
    with open(file, 'a') as f: # check a+ type
        f.write("\n\n================================")
        f.write("\n APPENDED INFO {}".format(datetime.now().strftime('%Y%m%d %H:%M:%S')))
        f.write("\n================================")
        f.write("\n.. info generated by <{}\{}>".format(_direc,_base))
        f.write("\n.. path to read the output files: {}".format(FOLDER_PATH))
        f.write("\n.. new generated files for:")
        f.write("\n   (folders) {}".format(algs))
    return

def __save_log_folder(algs):
    '''
    ADJUST FUNCTION:
        If the log file exists in the folder:
            - [false] print all information on the containing files
            - [true] append information on new generated files 
    '''
    base = str(os.path.basename(__file__))
    direc = str(os.path.dirname(__file__))
    
    log_file = SAVE_FOLDER+'_info.txt'
    exists = os.path.exists(log_file)
    
    if not exists:
        __write_info(log_file, direc, base)
        __append_info(log_file, algs, direc, base)
    else:
        __append_info(log_file, algs, direc, base)
    return

def __run(alg):
    
    folder = FOLDER_PATH.format(alg)
    save_path = __get_save_path(alg)
    
    if alg in ['esmamds-pop', 'esmamds-cpm', 'esmam-pop', 'esmam-cpm']:
        esmam_class.run(alg, folder, save_path)
    elif alg in ['bs-emm-pop', 'bs-emm-cpm', 'bs-sd-pop', 'bs-sd-cpm']:
        bs_class.run(alg, folder, save_path)
    elif alg == 'lr-rules':
        lr_rules.run(folder, save_path)
    elif alg == 'dssd-cbss':
        dssd_cbss.run(folder, save_path)
    return

def run(algs):
    '''
    Process algorithms output files and saves the processed log-files in 
    'EsmamDS/experiments/_results_processing/_processed_output_files/'

    Parameters
    ----------
    algs : list
        List of algorithm's names to process the results.
        (all options) ['esmamds-pop', 'esmamds-cpm', 'esmam-pop', 'esmam-cpm',
                       'bs-emm-pop', 'bs-emm-cpm', 'bs-sd-pop', 'bs-sd-cpm',
                       'lr-rules', 'dssd-cbss']

    Returns
    -------
    None.

    '''
    
    if len(algs)==1:
        __run(algs[0])
    else:
        for alg in algs:
            __run(alg)  
    return

if __name__ == '__main__':

    # parse args setting
    parser = argparse.ArgumentParser(description='Script to generate the processed log-files from the original algorithms output files.')
    
    parser.add_argument("--alg", type=str,
                        choices=['esmamds-pop', 'esmamds-cpm', 'esmam-pop', 'esmam-cpm',
                                 'bs-emm-pop', 'bs-emm-cpm', 'bs-sd-pop', 'bs-sd-cpm',
                                 'lr-rules', 'dssd-cbss'],
                        help="Algorithm name to single process the results.")
    args = parser.parse_args()
    
    if args.alg:
        algs = [args.alg]
    else:
        algs = ALGORITHMS
    
    run(algs)
    __save_log_folder(algs)