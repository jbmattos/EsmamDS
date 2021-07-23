# ICDE2022 - Empirical Evaluation

Here are presented the experiments conducted to evaluate the EsmamDS approach for mining local patterns associated with exceptional survival behaviours.
The goal is to provide comprehensible and diverse characterisation of subgroups presenting unusual KM models.
Its results were confronted with state of the art methods in the literature that provide characterisation over unusual survival behaviour.  

## Evaluation metrics

The empirical evaluation was designed to assess the EsmamDS findings with respect to the exceptionality, interpretability, generality and redundancy as global metrics for a set of discovered subgroups. Redundancy is assessed for descriptions, coverage, and survival models.
Additionally, metrics of similarity between pairs of subgroups were employed to compare different sets of subgroups with respect to similarities in description, coverage and survival models.  
All metrics are defined in Table I of the article.  
A brief introduction of the evaluation metrics is provided bellow.  

**Similarity** (between pair of subgroups)  
Description and coverage similarity are assessed based on a modified Jaccard index.  
(survival) Model similarity is assessed based on the log-rank test  
  
**Exceptionality**  
Proportion of exceptional subgroups (in a set)  
  
**Interpretability**  
Number of discovered subgroups  
Average size of subgroup description  
  
**Generality**  
Average (percentage) subgroup coverage  
Data set coverage  
  
**Redundancy**  
The normalised sum of the similarity measures (description, coverage, model) for all (unordered) pairs of subgroups in a set  
CR measure (Van Leeuwen and Knobbe, 2012)


## Experimental proceadure

The EsmamDS algorithm was implemented to mine exceptional subgroups considering two different (comparison) baselines:  
(i)  to compare the subgroup to its complement; or  
(ii) to compare to the population (whole data set)  
Such configuration is defined by the _baseline_ parameter. In the empirical evaluation, both configurations are considered separately.

For each baseline, the EsmamDS was performed 30 times over each data set.  


## EsmamDS hyper-parameters selection

The following hyper-parameters were defined with a randomised search considering the values: 
- _nAnts_ = {100, 200, 500, 1000, 3000};
- _minCov_ = {0.01, 0.02, 0.05, 0.1};
- _nConverg_ = {5, 10, 30};
- _maxStag_ = {20, 30, 40, 50};
- _L_ = {1,3,5,10};
- _W_ = 0.9; and
- _alpha_ = 0.05  

The EsmamDS was executed on actg320, breast-cancer and ptc data sets for (a sample of) 10\% of the total number of parameters' combinations.
The algorithm configurations were chosen by ordering all samples according to their average performance considering the following metrics' order: description redundancy; coverage redundancy; CR; model redundancy; subgroup (average) coverage; data set coverage; subgroup (average) description size; and number of discovered subgroups.  

The following configurations were selected for paramenter(value):
- Population baseline: _alpha_(0.05), _nAnts_(100), _minCov_(0.1), _nConverg_(5), _maxStag_(40), _W_(0.9), _L_(5)
- Complement baseline: _alpha_(0.05), _nAnts_ (100), _minCov_(0.05), _nConverg_(5), _maxStag_(40), _W_(0.9), _L_(10)

The configuration of other approaches confronted in the experimentation were defined for each data set (considering to their baseline) according to the results achieved by the EsmamDS in the experimental proceadure.
The following three parameters were adjusted as follows:
- (_minCov_) _Minimum coverage_: defined by the same parameter value chosen for the EsmamDS;
- (_bs_) _Beam-size_ (or maximum number of discovered subgroups): given by the average number of subgroups discovered by the EsmamDS in the 30 experiments;
- (_maxDepth_) _Rule-depth_ (or refinement/search depth): given by the average of the the maximum description size achieved during the EsmamDS execution in the 30 experiments.

The values used for the above parameters, for each data set and both baselines, are presented in the `_paramsConfig.csv` file.


## Results

The results and statistical analysis of the EsmamDS (and confronted approaches) performance for the _exceptionality_, _interpretability_, _generality_ and _redundancy_ metrics are presented in `metrics results (tables and statistics)` folder.

The `sets similarities (heatmap matrix)` folder presents the results for the analysis of the _similarity_ metrics.

The `survival models` folder present the analysis of the survival models related to the EsmamDS (and confronted approaches) discovered subgroups. 
