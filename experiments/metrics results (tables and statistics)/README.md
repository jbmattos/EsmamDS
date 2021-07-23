## Evaluation Metrics Results

The performances of the EsmamDS and confronted approaches for the _exceptionality_, _interpretability_, _generality_ and _redundancy_ metrics are presented in
```
metrics_baseline-complement.csv
metrics_baseline-population.csv
```
For the EsmamDS and Esmam algorithms, the tables present the complete sample of 420 results for each metric (30 experiments on 14 data sets).
For the remaining deterministic algorithms, the results were paired for each data set by repeating them 30 times.

### Statistical Analysis

For the metrics of _interpretability_, _generality_ and _redundancy_, statistical analysis of the results was performed by a (paired) Friedman test followed by a Nemenyi post-hoc test (for a level of significance of 5%).  
_(the generality metrics were maximised and the metrics of interpretability and redundancy were minimised)_

The `statistical_report.ipynb` notebook presents the p-value of the Friedman test and the results of Nemenyi post-hoc test in the form of Critical Distance plot.
