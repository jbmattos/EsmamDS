"""

"""
import argparse
import warnings
import os
import pathlib

from datetime import datetime
from _utils.results_cls import Table, SurvivalPlots, Heatmaps

warnings.filterwarnings("ignore")

ROOT = "EsmamDS"
ROOT_PATH = str(pathlib.Path(__file__).parent.absolute()).split(ROOT)[0].replace('\\','/') + ROOT + '/'
MAIN_FOLDER= ROOT_PATH + 'experiments/'

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
        f.write("\n.. call to: baseline={} and results={}".format(call[0], call[1]))
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

def __table(baseline):
    tbl = Table()
    tbl.set_table(baseline)
    tbl.save()
    return

def __survival_models(baseline):
    surv = SurvivalPlots()
    surv.full_report(baseline)
    return

def __set_similarity(baseline):
    sim = Heatmaps()
    sim.similarity_full_report(baseline)
    pass

def __set_redundancy(baseline):
    sim = Heatmaps()
    sim.redundancy_full_report(baseline)
    pass

def __pipeline(baselines, results):
    
    for base in baselines:
        print('\n...processing baseline={}\n'.format(base))
        
        if results in ['all', 'tbl']:
            __table(base)
        
        if results in ['all', 'surv']:
            __survival_models(base)
        
        if results in ['all', 'setsim']:
            __set_similarity(base)
        
        if results in ['all', 'setred']:
            __set_redundancy(base)
        
    return


if __name__ == '__main__':

    # parse args setting
    parser = argparse.ArgumentParser(description='Script to generate the final results of EsmamDS empirical evaluation.')
    
    parser.add_argument("--base", type=str,
                        choices=["all", "complement", "population"],
                        default="all",
                        help="Baseline for generating the results.")
    parser.add_argument("--res", type=str,
                        choices=["all", "tbl", "surv", "setsim", "setred"],
                        default="all",
                        help="The results to be generated.")
    args = parser.parse_args()
    
    if args.base == 'all':
        baselines = ['population', 'complement']
    else:
        baselines = [args.base]
    
    if args.res == 'all': res=["tbl", "surv", "setsim", "setred"]
    else: res=[args.res] 
    print('\n>> Generate results={} for baselines={}'.format(res, baselines))
    
    __pipeline(baselines, args.res)
    __save_log_folder((baselines,res))