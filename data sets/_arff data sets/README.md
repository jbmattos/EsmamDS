### Arff data sets
The .arff files regarding the final data sets.  
The `model target (EMM)` folder provides the data set files for algorithms using the KM model as target concept.  
For the files provided in `single target (SD)` folder, the censoring feature was removed (since it is not considered for the SD target concept and does not comprise a descriptive attribute).

The `_arff_generator.py` file provides the code for generating the .arff files from the [final data sets](https://github.com/jbmattos/EsmamDS/tree/icde2022/data%20sets/final%20data%20sets).  
When provided --db, the script process the single data set. Otherwise, all data sets in the final folder are processed.  
The store-true --sd parameter executes the processing for single target.
```
usage: _arff_generator.py [-h] [--db DB] [--sd]

Script to convert .xz datasets into .arff files.

optional arguments:
  -h, --help  show this help message and exit
  --db DB     Data set (file) name to single process
  --sd        Generate .arff files for Subgroup Discovery
```
