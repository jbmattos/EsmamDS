## Analisys of redundancy within subgroup set

Redundancy in a final set of subgroups is analysed assessing the three types of redundancy (description, coverage and model) within each resultant set of subgroups
provided by the compared algorithms.  
For each pair of different subgroups in the considered set, the similarity measures are computed. 
The redundancy metrics presented in the [metrics results](https://github.com/jbmattos/EsmamDS/tree/icde2022/experiments/metrics%20results%20(tables%20and%20statistics))
are the normalised sum of the (similarity) values ploted in this analysis.

The results for the three types of redundancy are presented in the files 
```
{baseline}_coverRedundancy.pdf
{baseline}_descrRedundancy.pdf
{baseline}_modelRedundancy.pdf
``` 
where:
- The pages provide the results for each data set;
- In each page, the columns of plots are the results of different algorithms, and the rows are the (30) runs.

In each triangular heatmap matrix, the dimensions are the subgroups provided by each algorithm in their final set os discovered subgroups. 
The main diagonal is supressed from the plot because it represents the comparison of a subgroup with itself.
