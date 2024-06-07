'''
Script to generate the report about the statistical tests performed 
over the evaluation metrics.
'''

import argparse
import math
import numpy as np
import pandas as pd
import pathlib
import warnings
from autorank import autorank, plot_stats
from matplotlib import pyplot as plt
from matplotlib.pylab import rcParams

warnings.filterwarnings("ignore")

THIS_PATH = str(pathlib.Path(__file__).parent.absolute())
TBL_POP_FILE = THIS_PATH+'\\metrics_baseline-population.csv'
TBL_CPM_FILE = THIS_PATH+'\\metrics_baseline-complement.csv'
METRICS = ['#sg', 'length', 'sgCov', 'setCov', 
           'description redundancy', 'cover redundancy',
           'CR', 'model redundancy']
STAT_ORDER = {'#sg': 'ascending',
              'length': 'ascending',
              'sgCov':'descending',
              'setCov': 'descending',
              'description redundancy': 'ascending',
              'cover redundancy': 'ascending',
              'CR': 'ascending',
              'model redundancy': 'ascending'}
ALPHA = 0.05

'''
Functions from autorank package
    Adjustment on < plot_stats > to provide CD plot reversed (from smaller ranks to higher - ascending order)
'''
def get_sorted_rank_groups(result, reverse):
    if reverse:
        names = result.rankdf.iloc[::-1].index.to_list()
        if result.cd is not None:
            sorted_ranks = result.rankdf.iloc[::-1].meanrank
            critical_difference = result.cd
        else:
            sorted_ranks = result.rankdf.iloc[::-1]['mean']
            critical_difference = (result.rankdf.ci_upper[0] - result.rankdf.ci_lower[0]) / 2
    else:
        names = result.rankdf.index.to_list()
        if result.cd is not None:
            sorted_ranks = result.rankdf.meanrank
            critical_difference = result.cd
        else:
            sorted_ranks = result.rankdf['mean']
            critical_difference = (result.rankdf.ci_upper[0] - result.rankdf.ci_lower[0]) / 2

    groups = []
    cur_max_j = -1
    for i in range(len(sorted_ranks)):
        max_j = None
        for j in range(i + 1, len(sorted_ranks)):
            if abs(sorted_ranks[i] - sorted_ranks[j]) <= critical_difference:
                max_j = j
                # print(i, j)
        if max_j is not None and max_j > cur_max_j:
            cur_max_j = max_j
            groups.append((i, max_j))
    return sorted_ranks, names, groups

def cd_diagram(result, reverse, ax, width):
    """
    Creates a Critical Distance diagram.
    """

    def plot_line(line, color='k', **kwargs):
        ax.plot([pos[0] / width for pos in line], [pos[1] / height for pos in line], color=color, **kwargs)

    def plot_text(x, y, s, *args, **kwargs):
        ax.text(x / width, y / height, s, *args, **kwargs)

    sorted_ranks, names, groups = get_sorted_rank_groups(result, reverse)
    cd = result.cd

    lowv = min(1, int(math.floor(min(sorted_ranks))))
    highv = max(len(sorted_ranks), int(math.ceil(max(sorted_ranks))))
    cline = 0.4
    textspace = 1
    scalewidth = width - 2 * textspace

    def rankpos(rank):
        if not reverse:
            relative_rank = rank - lowv
        else:
            relative_rank = highv - rank
        return textspace + scalewidth / (highv - lowv) * relative_rank

    linesblank = 0.2 + 0.2 + (len(groups) - 1) * 0.1

    # add scale
    distanceh = 0.25
    cline += distanceh

    # calculate height needed height of an image
    minnotsignificant = max(2 * 0.2, linesblank)
    height = cline + ((len(sorted_ranks) + 1) / 2) * 0.2 + minnotsignificant

    if ax is None:
        fig = plt.figure(figsize=(width, height))
        fig.set_facecolor('white')
        ax = fig.add_axes([0, 0, 1, 1])  # reverse y axis
    ax.set_axis_off()

    # Upper left corner is (0,0).
    ax.plot([0, 1], [0, 1], c="w")
    ax.set_xlim(0, 1)
    ax.set_ylim(1, 0)

    plot_line([(textspace, cline), (width - textspace, cline)], linewidth=0.7)

    bigtick = 0.1
    smalltick = 0.05

    tick = None
    for a in list(np.arange(lowv, highv, 0.5)) + [highv]:
        tick = smalltick
        if a == int(a):
            tick = bigtick
        plot_line([(rankpos(a), cline - tick / 2),
                   (rankpos(a), cline)],
                  linewidth=0.7)

    for a in range(lowv, highv + 1):
        plot_text(rankpos(a), cline - tick / 2 - 0.05, str(a),
                  ha="center", va="bottom")

    for i in range(math.ceil(len(sorted_ranks) / 2)):
        chei = cline + minnotsignificant + i * 0.2
        plot_line([(rankpos(sorted_ranks[i]), cline),
                   (rankpos(sorted_ranks[i]), chei),
                   (textspace - 0.1, chei)],
                  linewidth=0.7)
        plot_text(textspace - 0.2, chei, names[i], ha="right", va="center")

    for i in range(math.ceil(len(sorted_ranks) / 2), len(sorted_ranks)):
        chei = cline + minnotsignificant + (len(sorted_ranks) - i - 1) * 0.2
        plot_line([(rankpos(sorted_ranks[i]), cline),
                   (rankpos(sorted_ranks[i]), chei),
                   (textspace + scalewidth + 0.1, chei)],
                  linewidth=0.7)
        plot_text(textspace + scalewidth + 0.2, chei, names[i],
                  ha="left", va="center")

    # upper scale
    if not reverse:
        begin, end = rankpos(lowv), rankpos(lowv + cd)
    else:
        begin, end = rankpos(highv), rankpos(highv - cd)

    plot_line([(begin, distanceh), (end, distanceh)], linewidth=0.7)
    plot_line([(begin, distanceh + bigtick / 2),
               (begin, distanceh - bigtick / 2)],
              linewidth=0.7)
    plot_line([(end, distanceh + bigtick / 2),
               (end, distanceh - bigtick / 2)],
              linewidth=0.7)
    plot_text((begin + end) / 2, distanceh - 0.05, "CD",
              ha="center", va="bottom")

    # no-significance lines
    side = 0.05
    no_sig_height = 0.1
    start = cline + 0.2
    for l, r in groups:
        plot_line([(rankpos(sorted_ranks[l]) - side, start),
                   (rankpos(sorted_ranks[r]) + side, start)],
                  linewidth=2.5)
        start += no_sig_height

    return ax

def plot_stats(result, reverse=False, allow_insignificant=False, ax=None, width=None):

    #if not isinstance(result, RankResult):
    #    raise TypeError("result must be of type RankResult and should be the outcome of calling the autorank function.")

    if result.omnibus == 'bayes':
        raise ValueError("ploting results of bayesian analysis not yet supported.")

    if result.pvalue >= result.alpha and not allow_insignificant:
        raise ValueError(
            "result is not significant and results of the plot may be misleading. If you want to create the plot "
            "regardless, use the allow_insignificant parameter to suppress this exception.")

    if ax is not None and width is not None:
        warnings.warn('width may be ignored because ax is defined.')
    if width is None:
        width = 6

    if result.omnibus == 'ttest':
        ax = ci_plot(result, True, ax, width)
    elif result.omnibus == 'wilcoxon':
        warnings.warn('No plot to visualize statistics for Wilcoxon test available. Doing nothing.')
    elif result.posthoc == 'tukeyhsd':
        ax = ci_plot(result, True, ax, width)
    elif result.posthoc == 'nemenyi':
        ax = cd_diagram(result, reverse, ax, width)
    return ax

'''
REPORT
'''
def report(file, _save, _ext, _base):
    
    df = pd.read_csv(file, header=[0,1], index_col=[0,1])
    rcParams["font.size"] = 16
    rcParams["font.family"] = "Times New Roman"
    
    for metric in METRICS:
        stat = autorank(df[metric], alpha=ALPHA, order=STAT_ORDER[metric], verbose=False)
        
        print('\033[1m' + "\n\n>> METRIC: {}".format(metric))
        print("Statistical test: {}".format(stat.omnibus))
        print("p-value: {}".format(stat.pvalue))
        print("Post-hoc test: {}".format(stat.posthoc))
        plot_stats(stat)
        
        if isinstance(_save, str):
            plt.savefig(_save +'\\_esmamds_CDplots_{}_{}.{}'.format(metric, _base, _ext), bbox_inches='tight') 
        plt.show()

if __name__ == '__main__':

    # ARG PARSE SETTINGS
    parser = argparse.ArgumentParser(description="Script to generate report on the statistical test performed over evaluation metrics")
    parser.add_argument("--save_path", type=str, default=None,
                        help="Path to save the CD plots.")
    parser.add_argument("--ext", type=str, default='pdf',
                        help="Extension to save the CD plots.")
        
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--pop', action='store_true', help="Report for baseline-population algorithms")
    group.add_argument('--cpm', action='store_true', help="Report for baseline-complement algorithms")
    
    args = parser.parse_args()
    
    if args.pop:
        report(TBL_POP_FILE, args.save_path, args.ext, 'pop')
    if args.cpm:
        report(TBL_CPM_FILE, args.save_path, args.ext, 'cpm')
