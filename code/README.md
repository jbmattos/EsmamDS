# EsmamDS algorithm
This folder provides the code for running the EsmamDS algorithm. It runs over a data set and provides the descriptive characterisation and KM models of subgroups of data observations that present unusual survival model (when comparing to a baseline model).

## esmamds.py usage
_This file should run inside the `EsmamDS/code/` folder._

```
usage: esmamds.py [-h] --f F [--dtypes DTYPES]
                  [--baseline {population,complement}] [--a A]
                  [--maxStag MAXSTAG] [--nAnts NANTS] [--nConverg NCONVERG]
                  [--minCov MINCOV] [--l L] [--w W] [--seed SEED] [--log]

EsmamDS: Exceptional Survival Model Ant Miner - Diverse Search

required arguments:
  --f F                 Data set .csv.xz file path [type:string]

optional arguments:
  --dtypes DTYPES       Json dtypes file path [type:string]
  --baseline            Baseline for subgroup comparison [default:population]
  {population,complement}
  --a A                 (Alpha) Level of significance [type:float]
  --maxStag MAXSTAG     Maximum stagnation of the algorithm [type:int]
  --nAnts NANTS         Size of the ant colony [type:int]
  --nConverg NCONVERG   Number of similar patters for convergence [type:int]
  --minCov MINCOV       Minimum subgroup (percentage) coverage [type:float]
  --l L                 Logistic function offset (description attenuation) [type:int]
  --w W                 Weight parameter [type:float]
  --seed SEED           Numpy random seed
  --log                 Saves (output) log file [type:store-true]
```
**Considerations**:  
- The data file should be a .csv.xz format
- It is advisible to provide a dtypes file to guarantee proper data representation. All descriptive attributes should be categorical; the survival time feature should be numerical and named 'survival_time'; and the censoring feature should be boolean (with 1 indicating the event) and named 'survival_status'.
- The default configuration of the hyper-parameters (_a, maxStag, nAnts, nConverg, minCov, l, w_) changes according to the --baseline.
- When provided --log, a nested dictionary containing a variety of information regarding the algorithm's execution is saved in json format.

### output files
A folder `EsmamDS/EsmamDS_exe{DATESTAMP}` is automatically created with the following output files:
- **_SurvivalModels.csv**: provides a table with the KM models of the discovered subgroups (and population).
- **_RuleModel.csv**: provides a table with general information on the discovered subgroups.
- **_RuleSet.txt**: provides the subgroups' characterisations grouped by model and description similarity.
- (optional)**_log.json**: nested dictionary containing information regarding the execution
