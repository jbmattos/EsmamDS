## Data processing scripts

### data_analysis.ipynb
Jupyter nootebook that comprises a manual and interactive data analysis proceadure to extract information for the data processing. The notebook generates a json file as output in the format to be used in the "data_preprocessing.py" script. The json files referring to all data sets are provided in the `_prep files` folder.

### data_preprocessing.py
Pipeline employed to process and discretise the data sets.  
The script reads data from `data sets/_original data sets/` and saves the processed data in `data sets/final data sets/`.  
The *_prep files* are used in the processing pipeline.  
A log file containing information from the discretisation is saved in the `_disc logs` folder.  
When provided --db, the script process the single data set. Otherwise, all data sets in the final folder are processed.  
The store-false --nodisc parameter deactivates the discretisation process.
```
usage: data_preprocessing.py [-h] [--db DB] [--nodisc]

Script to perform data preprocessing and discretisation.

optional arguments:
  -h, --help  show this help message and exit
  --db DB     Data set (file) name to single process
  --nodisc    Do not perform data discretisation
```
