# -*- coding: utf-8 -*-
"""

"""
import argparse
import warnings
#import errno
#import json
#import os
#import pandas as pd
#import pathlib

#from datetime import datetime
from _utils.results_cls import Table

warnings.filterwarnings("ignore")

def __table(baseline):
    # table
    tbl = Table()
    tbl.set_table(baseline)
    tbl.save()
    return

def __survival_models(baseline):
    pass

def __set_similarity(baseline):
    pass

def __set_redundancy(baseline):
    pass

def __pipeline(baselines, results):
    
    for base in baselines:
        print('\n...processing baseline={}'.format(base))
        
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
    
    '''
    ADJUST
    receive baseline options [population, complement, default=both]
    receive results options  [table, survival, setSimilarity, setRedundancy, default=all]
    '''
    # parse args setting
    parser = argparse.ArgumentParser(description='Script to generate the final results of EsmamDS empirical evaluation.')
    
    parser.add_argument("--base", type=str,
                        choices=["all", "complement", "population"],
                        default="all",
                        help="Baseline for generating the results.")
    parser.add_argument("--res", type=str,
                        choices=["all", "tbl", "surv", "setsim", "setred"],
                        default="all",
                        help="Baseline for generating the results.")
    args = parser.parse_args()
    
    if args.base == 'all':
        baselines = ['population', 'complement']
    else:
        baselines = [args.base]
    
    if args.res == 'all': res=["tbl", "surv", "setsim", "setred"]
    else: res=[args.res] 
    print('\n>> Generate resuls={} for baselines={}'.format(res, baselines))
    
    __pipeline(baselines, args.res)