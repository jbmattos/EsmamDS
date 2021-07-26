### Arff data sets 

The .arff files regarding the final data sets. For the files provided in `_arff single target (SD)`, the _censoring_ feature was removed (since it is not considered for the SD target concept and does not comprise a descriptive attribute). 

The `_arff_generator.py` file provides the code for generating the .arff files from the final data sets.  
When provided --db, the script process the single <db> data set. Otherwise, all data sets in the final folder are processed.  
The store-true --sd parameter executes the processing for single target.
```
usage: _arff_generator.py [-h] [--db DB] [--sd]

Script to convert .xz datasets into .arff files.

optional arguments:
  -h, --help  show this help message and exit
  --db DB     Data set (file) name to single process
  --sd        Generate .arff files for Subgroup Discovery
```
