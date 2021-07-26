'''
SCRIPT FOR GENERATING .ARFF FILES

This script process the entire < EsmamDS/data sets/final data sets > folder
and saves the processed data sets into < EsmamDS/data sets/_arff data sets >

Alternativelly, to process a single data set use --db 
providing the data set name (must be a file in <final data sets> folder)

>> Reference for script implementation:
https://stackoverflow.com/questions/48993918/exporting-dataframe-to-arff-file-python
'''

import arff
import argparse
import glob
import json
import os
import pandas as pd

# Global variables
VAR_TIME_NAME = 'survival_time'
VAR_EVENT_NAME = 'survival_status'

DATASETS = ['actg320','breast-cancer','cancer','carcinoma','gbsg2','lung','melanoma','mgus2','mgus','pbc','ptc','uis','veteran','whas500']

# paths
ROOT = "EsmamDS"
PATH = os.path.dirname(__file__).split(ROOT)[0]+ROOT # read system path up to ROOT
DATA_PATH = PATH+"\\data sets\\"
READ_PATH = DATA_PATH+"final data sets\\"
SAVE_PATH = DATA_PATH+"_arff data sets\\"

ARFF_FILE = '{}_disc.arff'


def __read_data(db_name):

    data_file = READ_PATH+"{}_disc.xz".format(db_name)
    json_file = READ_PATH+"{}_disc_dtypes.json".format(db_name)
    with open(json_file, 'r') as f:
        dtypes = json.load(f)
    return pd.read_csv(data_file, delimiter=',', header=0, index_col=False, compression='xz', dtype=dtypes)


def generate_arff_file(db_name, _sd):
    
    # read database
    db = __read_data(db_name)
    db[VAR_EVENT_NAME] = db[VAR_EVENT_NAME].astype('int') # changes censoring attribute from <bool> to <int (0,1)> (required by LRrules)

    # ADJUSTING DATA TO ARFF FORMAT
    colTime = db[VAR_TIME_NAME].copy()
    colEvent = db[VAR_EVENT_NAME].copy()
    db_arff = db.drop(columns=[VAR_TIME_NAME,VAR_EVENT_NAME])

    # GENERATING ARFF @ATTRIBUTE
    ctgCols = list(db_arff.columns)
    if _sd:
        attributes = [(VAR_TIME_NAME, 'NUMERIC')]
    else:
        attributes = [(VAR_TIME_NAME, 'NUMERIC'),(VAR_EVENT_NAME, 'NUMERIC')]
    attributes += [(col, db[col].unique().astype(str).tolist()) for col in ctgCols]

    # GENERATING ARFF @DATA
    # data is a list in wich each db_arff row is a list
    if _sd:
        data = [ [colTime.loc[i]] + db_arff.loc[i].tolist() for i in range(db_arff.shape[0])]
    else:
        data = [ [colTime.loc[i]] + [colEvent.loc[i]] + db_arff.loc[i].tolist() for i in range(db_arff.shape[0])]

    # SAVING ARFF FILE
    save_name = SAVE_PATH + ARFF_FILE.format(db_name)
    
    arff_dict = {
        'attributes': attributes,
        'data': data,
        'relation': READ_PATH+"{}_disc.xz".format(db_name)}
    with open(save_name, "w", encoding="utf8") as f:
         arff.dump(arff_dict, f)
    print('..saved: {}'.format(save_name))
    return

if __name__ == '__main__':
    
     # parse args setting
    parser = argparse.ArgumentParser(description='Script to convert .xz datasets into .arff files.')
    
    parser.add_argument("--db", type=str,
                        help="Data set (file) name to single process")
    parser.add_argument("--sd", action="store_true",
                        help="Generate .arff files for Subgroup Discovery (remove Survival Event feature)")
    args = parser.parse_args()
    
    if args.db:
        generate_arff_file(args.db, args.sd)
    
    else:
        
        files_path = READ_PATH+'*_disc.xz'
        for file in glob.iglob(files_path):
            
            print('>> process file: {}'.format(file))
            db_name = file.split("\\")[-1].split('.')[0]
            generate_arff_file(db_name, args.sd)
