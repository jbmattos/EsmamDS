import argparse
import errno
import json
import os
import pathlib

from algorithm import EsmamDS
from datetime import datetime

"""
ADJUST PATHS
"""
ROOT = "EsmamDS"
ROOT_PATH = str(pathlib.Path(__file__).parent.absolute()).split(ROOT)[0] + ROOT + '/'
CODE_PATH = ROOT_PATH + 'code/'
DATA_PATH = ROOT_PATH + 'data sets/'
SAVE_PATH = ROOT_PATH + 'EsmamDS_statExe{}/'.format(datetime.now().strftime('%Y%m%d'))
DB_NAMES = ['actg320','breast-cancer','cancer','carcinoma','gbsg2','lung','melanoma','mgus2','mgus','pbc','ptc','uis','veteran','whas500']

PARAMS_COMPLEMENT = {'no_of_ants': 100,
                     'min_size_subgroup': 0.05,
                     'no_rules_converg': 5,
                     'its_to_stagnation': 40,
                     'logistic_offset': 10}

# creates directory for saving results and logs
if not os.path.exists(os.path.dirname(SAVE_PATH)):
    try:
        os.makedirs(os.path.dirname(SAVE_PATH))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise
# read the seeds
with open(CODE_PATH+'_utils/_seeds.json', 'r') as f:
    SEEDS = json.load(f)


def stats_results(sg_baseline, _save_log=False):

    if sg_baseline == 'population': comp = 'pop'
    else: comp = 'cpm'

    for db_name in DB_NAMES:
        print('\n\n>> EsmamDS-{} >> database: {}'.format(comp,db_name))

        # ADJUST FILE PATHS
        db_path = None
        dtypes_path = None

        # select the seeds
        seeds = SEEDS[db_name]

        for exp in range(30):  # statistical experiments
            print('..exp {}'.format(exp))
            save_name = SAVE_PATH + 'EsmamDS-{}_{}_exp{}'.format(comp, db_name, exp)

            if sg_baseline == 'complement':
                run(file_path=db_path, dtypes_path=dtypes_path, sg_baseline=sg_baseline, seed=seeds[exp], save_path=save_name,
                    _save_log=_save_log, **PARAMS_COMPLEMENT)
            else:
                run(file_path=db_path, dtypes_path=dtypes_path, sg_baseline=sg_baseline, seed=seeds[exp], save_path=save_name,
                    _save_log=_save_log)
    return


def run(file_path, dtypes_path, sg_baseline, seed, save_path, _save_log, **kwargs):

    esmam = EsmamDS(sg_baseline=sg_baseline, seed=seed, **kwargs)
    esmam.read_data(file_path, dtypes_path)
    esmam.fit()
    esmam.save_results(save_path)
    if _save_log:
        esmam.save_logs(save_path)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='EsmamDS experiments script')
    parser.add_argument("--log", action='store_true',
                        help="Saves (output) log file")
    args = parser.parse_args()

    print('\n ESMAM-DS EMPIRICAL EVALUATION')
    print('.. this call will execute the EsmamDS algorithm for baseline={complement,population} on 14dbs/30exp')

    stats_results(sg_baseline='complement', _save_log=args.log)
    stats_results(sg_baseline='population', _save_log=args.log)
