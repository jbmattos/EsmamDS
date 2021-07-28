import numpy as np
import math


class Term:

    def __init__(self, attribute, value, dataset, min_case_per_rule):
        self.attribute = attribute
        self.value = value
        self.covered_cases = self._get_cases(dataset.data, dataset.get_col_index(attribute), value)
        self._heuristic = None
        self._min_cases = min_case_per_rule
        self.set_heuristic(dataset)

    @staticmethod
    def _get_cases(np_data, attribute_idx, value):
        return list(np.asarray(np_data[:, attribute_idx] == value).nonzero()[0])

    def set_heuristic(self, dataset):
        # 2 classes: 0 for survival time < average_survival, 1 otherwise

        data = dataset.get_data().loc[dataset.get_uncovered_cases()].copy()
        survival_average = data[dataset.surv_name].mean()

        data_attr_val = data[data[self.attribute] == self.value].copy()
        covered_rows = list(data_attr_val.index)
        if len(covered_rows) < self._min_cases:
            self._heuristic = float(0)
            return

        # A POSTERIORI PROBABILITY: P(W|A=V)
        n_classes = 2
        classes = (data_attr_val[dataset.surv_name] >= survival_average).astype(str)
        prob_posteriori = classes.value_counts(normalize=True).to_dict()

        # ENTROPY = -SUM_FOR_ALL_CLASSES[ P(W|A=V) * log2(P(W|A=V)) ]
        entropy = 0
        for w in prob_posteriori:
            if prob_posteriori[w] != 0:
                entropy -= prob_posteriori[w] * math.log2(prob_posteriori[w])

        self._heuristic = math.log2(n_classes) - entropy
        return

    def get_heuristic(self):
        return self._heuristic
