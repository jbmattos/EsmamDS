'''
SCRIPT FOR EXTRACTING THE SEARCH PARAMETERS <MAX_DEPTH, NUM_RULES, MIN_COV> 
FROM ESMAM-DS STATISTICAL EXECUTIONS:
    - max_depth: the average maximum depth rounded to ceil, for each dataset (on the 30 executions)
    - max_depth: the average number of rules rounded to ceil, for each dataset (on the 30 executions)
    - min_coverage: the minimum rule coverage is a percentage (EsmamDS parameter) of the dataset size, rounded to ceil
    
SCRIPT INPUT PARAMETER:
    The folder containing the statistical executions of EsmamDS for all datasets
'''

import argparse
from datetime import datetime
import json
import math
import numpy as np
import os
import pandas as pd
import pathlib
import zipfile

ROOT = "EsmamDS"
ROOT_PATH = str(pathlib.Path(__file__).parent.absolute()).split(ROOT)[0] + ROOT + '/'
MAIN_FOLDER = ROOT_PATH + 'experiments/'

LOG_FILE_NAME = "/EsmamDS-{}_{}_exp{}_log.json"
DATASETS = ['actg320','breast-cancer','cancer','carcinoma','gbsg2','lung','melanoma','mgus2','mgus','pbc','ptc','uis','veteran','whas500']
N_RUNS = 30

def __write_info(file, _direc, _base):
    
    with open(file, 'w') as f:
        f.write("This file was automatically created by the execution of")
        f.write("\n{}/{} {}".format(_direc,_base, datetime.now().strftime('%Y%m%d %H:%M:%S')))
    return

def __append_info(file, call, _direc, _base):
    with open(file, 'a') as f: # check a+ type
        f.write("\n\n=================================")
        f.write("\n APPENDED INFO {}".format(datetime.now().strftime('%Y%m%d %H:%M:%S')))
        f.write("\n=================================")
        f.write("\n.. exe {}/{}".format(_direc,_base))
        f.write("\n.. call to: baseline={}".format(call[0]))
        f.write("\n.. esmamds results read from: {}".format(call[1]))
    return

def __save_log_folder(call):
    base = str(os.path.basename(__file__))
    direc = str(os.path.dirname(__file__))
    
    log_file = MAIN_FOLDER +'_info.txt'
    exists = os.path.exists(log_file)
    if not exists:
        __write_info(log_file, direc, base)
        __append_info(log_file, call, direc, base)
    else:
        __append_info(log_file, call, direc, base)
    return


def pipeline(baseline, logs_path):
    
    depth = {}.fromkeys(DATASETS)
    rules = {}.fromkeys(DATASETS)
    min_cov = {}.fromkeys(DATASETS)
    
    # iterates over datasets 
    for db in DATASETS:

        max_depth = []
        num_rules = []
        
        # iterates over each data set experiment (0-29)
        for exp in range(N_RUNS):
            
            # read EsmamDS log_results
            log_file = logs_path + LOG_FILE_NAME.format(baseline,db,exp)
            if os.path.exists(log_file): # exists a .json file
                with open(log_file, 'r') as f:
                    log = json.load(f)
            else: # dont exist a .json file >> read (.json).zip file
                zip_path = log_file.split('.')[0] + '.zip'
                with zipfile.ZipFile(zip_path) as z:
                    for filename in z.namelist():
                        with z.open(filename) as f:  
                            data = f.read()  
                            log = json.loads(data)
            
            # compute maximum depth of all outter-ant-iterations in a single experiment
            max_exp = 0            
            for outter_it, inner_its_dic in log['run']['log_iterations'].items():
                depths = list(map(lambda inner_dic_val: len(inner_dic_val['r_const']), inner_its_dic.values()))     # list of length of all constructed rules in a colony
                maximum = max(depths)                                                                               # size of the large constructed rule
                if maximum > max_exp:
                    max_exp = maximum

            max_depth.append(max_exp)                                                                               # maximum rule length in a db/experiment
            num_rules.append(log['metrics']['num_rules'])                                                           # number of discovered rules in a db/experiment
            min_cov[db] = math.ceil(log['params']['sg_min_percent'] * log['run']['data_shape'][0])                  # minimum rule-coverage for the db (all experiments have same min_cov)

        depth[db] = max_depth[:]    # maximum depths for all db-experiments
        rules[db] = num_rules[:]    # number of rules for all db-experiments

    dic_params = {'max_depth': pd.DataFrame(depth).mean().apply(np.ceil).apply(int),
                  'num_rules': pd.DataFrame(rules).mean().apply(np.ceil).apply(int),
                  'min_cov':   pd.Series(min_cov)}
    
    # saves table
    if baseline == 'pop':
        file_name = 'population_paramsConfig.csv'
    elif baseline == 'cpm':
        file_name = 'complement_paramsConfig.csv'
    save_file = MAIN_FOLDER + '{}'.format(file_name)
    pd.DataFrame(dic_params).to_csv(save_file)
    print('.. saved: {}'.format(save_file))
    
    return

if __name__ == '__main__':
    
    # parse args setting
    parser = argparse.ArgumentParser(description="Scritp for extracting < depth, num_rules > parameters from EsmamDS statistical executions")
    
    parser.add_argument("--p", type=str,
                        help="Path of <esmamds> folder with the logs (json files) of the statistical runs. If not provided, use the repository default path.")
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pop", action='store_true',
                        help="<Population> baseline results.")
    group.add_argument("--cpm", action='store_true',
                        help="<Complement> baseline results.")
    args = parser.parse_args()
    
    if args.pop:
        baseline = 'pop'
    else:
        baseline = 'cpm'
    
    if args.p:
        read_path = args.p
    else:
        read_path = ROOT_PATH + 'experiments/_results_processing/_algorithms_output_files/esmamds-{}/'.format(baseline)
    
    print('\nExtracting paramenters from \n>> {}'.format(read_path))
    pipeline(baseline, read_path)
    __save_log_folder((baseline,read_path))