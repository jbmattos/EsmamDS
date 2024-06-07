"""
RANDOMISED PARAMETERS' SEARCH

This script executes the randomised search for defining the
EsmamDS configuration for both baselines, population and complement.

"""
import argparse
import errno
import json
import math
import os
import pathlib
import random

from algorithm import EsmamDS
from datetime import datetime

ROOT = "EsmamDS"
ROOT_PATH = str(pathlib.Path(__file__).parent.absolute()).split(ROOT)[0] + ROOT + '/'
CODE_PATH = ROOT_PATH + 'code/'
DATA_PATH = ROOT_PATH + 'data sets/final data sets/'
SAVE_PATH = ROOT_PATH + 'EsmamDS_randExe{}/'.format(datetime.now().strftime('%Y%m%d'))
DB_NAMES = ['actg320', 'breast-cancer', 'ptc']


def __generate_pool_sample():
    '''
    This function is the documentation of the pipeline employed for generating
    the set of EsmamDS parameters' configuration considered in
    the random search < randomised_search_parameters() >.
    '''
    no_of_ants = [100, 200, 500, 1000, 3000]
    min_size_subgroup = [0.01, 0.02, 0.05, 0.1]
    no_rules_converg = [5, 10, 30]
    its_to_stagnation = [20, 30, 40, 50]
    logistic_offset = [1, 3, 5, 10]

    n_pool = 0.1
    params = ['no_of_ants', 'min_size_subgroup', 'no_rules_converg', 'its_to_stagnation', 'logistic_offset']

    pool = [(ants, size, converg, stag, log) for ants in no_of_ants
            for size in min_size_subgroup
            for converg in no_rules_converg
            for stag in its_to_stagnation
            for log in logistic_offset]
    # with open(SAVE_PATH + '_pool_params.json', 'w') as f:
    #    json.dump(pool, f)

    sample_size = math.ceil(len(pool) * n_pool)
    samples = random.sample(pool, sample_size)

    all_configs = {}
    for idx, sample in enumerate(samples):
        all_configs[idx] = dict(zip(params, sample))

    with open(CODE_PATH+'_utils/_random_samples.json', 'w') as f:
        json.dump(all_configs, f)
    return


def randomised_search_parameters():

    with open(CODE_PATH + '_utils/_random_samples.json', 'r') as f:
        params = json.load(f)

    # iterates over params config samples
    for i in range(len(params)):
        print('\n\n>> SAMPLE #{}'.format(i))

        # creates directory for saving results and logs of each sample
        save_path = SAVE_PATH + 'esmamds_sample{}/'.format(i)
        if not os.path.exists(os.path.dirname(save_path)):
            try:
                os.makedirs(os.path.dirname(save_path))
            except OSError as exc:  # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise

        # read sample configuration
        config = params[str(i)]

        # iterates over datasets
        for db_name in DB_NAMES:

            db_path = DATA_PATH + "{}_disc.xz".format(db_name)
            dtypes_path = DATA_PATH + "{}_disc_dtypes.json".format(db_name)
            seeds = SEEDS[db_name]

            # run POPULATION-baseline
            sg_baseline = 'population'
            comp = 'pop'
            print('..EsmamDS-{} >> database: {}'.format(comp, db_name))
            save_name = save_path + 'EsmamDS-{}_{}'.format(comp, db_name)
            run(file_path=db_path, dtypes_path=dtypes_path, sg_baseline=sg_baseline, seed=seeds[0], save_path=save_name,
                _save_log=True, **config)

            # run COMPLEMENT-baseline
            sg_baseline = 'complement'
            comp = 'cpm'
            print('..EsmamDS-{} >> database: {}'.format(comp, db_name))
            save_name = save_path + 'EsmamDS-{}_{}'.format(comp, db_name)
            run(file_path=db_path, dtypes_path=dtypes_path, sg_baseline=sg_baseline, seed=seeds[0], save_path=save_name,
                _save_log=True, **config)
    return


def run(file_path, dtypes_path, sg_baseline, seed, save_path, _save_log, **kwargs):

    esmam = EsmamDS(sg_baseline=sg_baseline, seed=seed, **kwargs)
    esmam.read_data(file_path, dtypes_path)
    esmam.fit()
    esmam.save_results(save_path)
    if _save_log:
        esmam.save_logs(save_path)
    return


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='EsmamDS randomised search of parameters')
    args = parser.parse_args()

    # creates directory for saving results and logs
    if not os.path.exists(os.path.dirname(SAVE_PATH)):
        try:
            os.makedirs(os.path.dirname(SAVE_PATH))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    # read the seeds
    with open(CODE_PATH + '_utils/_seeds.json', 'r') as f:
        SEEDS = json.load(f)

    print('\n ESMAM-DS RANDOMISED PARAMETERS SEARCH')
    print('.. this call will execute the EsmamDS algorithm for baseline={complement,population} on 3 datasets '
          '(actg320, breast-cancer, ptc) and 96 samples of parameters configurations.')

    randomised_search_parameters()
