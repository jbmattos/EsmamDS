"""
This script executes the EsmamDS algorithm over a data set

.. this file should be run from EsmamDS/code/ folder
"""

import argparse
import errno
import os
import pathlib

from algorithm import EsmamDS
from datetime import datetime

# paths
ROOT = "EsmamDS"
SAVE_PATH = str(pathlib.Path(__file__).parent.absolute()).split(ROOT)[0] + ROOT + '/EsmamDS_exe{}/'.format(datetime.now().strftime('%Y%m%d'))

# default params
PARAMS_POPULATION = {'alpha': 0.05,
                     'its_to_stagnation': 40,
                     'no_of_ants': 100,
                     'no_rules_converg': 5,
                     'min_size_subgroup': 0.1,
                     'logistic_offset': 5,
                     'weigh_score': 0.9}
PARAMS_COMPLEMENT = {'alpha': 0.05,
                     'its_to_stagnation': 40,
                     'no_of_ants': 100,
                     'no_rules_converg': 5,
                     'min_size_subgroup': 0.05,
                     'logistic_offset': 10,
                     'weigh_score': 0.9}


def __set_params(args):

    if args.baseline == "population":
        config = PARAMS_POPULATION
    else:
        config = PARAMS_COMPLEMENT

    if args.a:
        config['alpha'] = args.a
    if args.maxStag:
        config['its_to_stagnation'] = args.maxStag
    if args.nAnts:
        config['no_of_ants'] = args.nAnts
    if args.nConverg:
        config['no_rules_converg'] = args.nConverg
    if args.minCov:
        config['min_size_subgroup'] = args.minCov
    if args.l:
        config['logistic_offset'] = args.l
    if args.w:
        config['weigh_score'] = args.w

    return config


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='EsmamDS: Exceptional Survival Model Ant Miner - Diverse Search')
    parser.add_argument("--f", type=str, required=True,
                        help="Data set file path")
    parser.add_argument("--dtypes", type=str,
                        help="Json dtypes file path ")
    parser.add_argument("--baseline",
                        choices=["population", "complement"],
                        type=str, default="population",
                        help="Baseline for subgroup comparison")
    parser.add_argument("--a", type=float,
                        help="(Alpha) Level of significance")
    parser.add_argument("--maxStag", type=int,
                        help="Maximum stagnation of the algorithm")
    parser.add_argument("--nAnts", type=int,
                        help="Size of the ant colony")
    parser.add_argument("--nConverg", type=int,
                        help="Number of similar patters for convergence")
    parser.add_argument("--minCov", type=float,
                        help="Minimum subgroup coverage")
    parser.add_argument("--l", type=int,
                        help="Logistic function offset (description attenuation)")
    parser.add_argument("--w", type=float,
                        help="Weight parameter")
    parser.add_argument("--seed", type=int, default=0,
                        help="Numpy random seed")
    parser.add_argument("--log", action='store_true',
                        help="Saves (output) log file")
    args = parser.parse_args()

    params = __set_params(args)
    
    # creates directory for saving results and logs
    if not os.path.exists(os.path.dirname(SAVE_PATH)):
        try:
            os.makedirs(os.path.dirname(SAVE_PATH))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise


    if args.baseline == "population":
        save_name = SAVE_PATH + 'EsmamDS-pop'
    else:
        save_name = SAVE_PATH + 'EsmamDS-cpm'

    alg = EsmamDS(sg_baseline=args.baseline, seed=args.seed, **params)
    alg.read_data(args.f, args.dtypes)
    alg.fit()
    alg.save_results(save_name)
    if args.log:
        alg.save_logs(save_name)
