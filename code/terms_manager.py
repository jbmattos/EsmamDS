import numpy as np
import pandas as pd
import math
import copy
from functools import reduce

from term import Term


class TermsManager:
    NP_INT = (np.int_, np.intc, np.intp, np.int8,
              np.int16, np.int32, np.int64, np.uint8,
              np.uint16, np.uint32, np.uint64)

    def __init__(self, dataset, min_case_per_rule, seed):
        self._terms = {}
        self._attr_values = {}
        self._availability = {}
        self._pheromone_table = {}
        self._heuristic_table = {}
        self._count_table = {}
        self._logistic_count_table = {}
        self._no_of_terms = 0
        self._Dataset = dataset
        np.random.seed(seed)

        # build object
        self._constructor(dataset, min_case_per_rule)

    def _constructor(self, dataset, min_case_per_rule):

        attr_values = dataset.attr_values.copy()    # Attribute-Values from the entire dataset
        heuristic_accum = 0

        # TERMS:
        # constructs _terms, _availability and _attr_values
        for attr, values in attr_values.items():
            values_terms = {}
            for value in values:
                term_obj = Term(attr, value, dataset, min_case_per_rule)
                values_terms[value] = term_obj
                self._no_of_terms += 1
                heuristic_accum += term_obj.get_heuristic()

            self._terms[attr] = values_terms.copy()
            self._attr_values[attr] = values[:]
            self._availability[attr] = True

        # TABLES:
        # _pheromone_table: {Attr : {Value : Pheromone}} | _heuristic_table: {Attr : {Value : Heuristic}}
        initial_pheromone = 1 / self._no_of_terms
        for attr, values in self._attr_values.items():
            self._pheromone_table[attr] = {}.fromkeys(values, initial_pheromone)
            self._count_table[attr] = {}.fromkeys(values, 0)
            self._logistic_count_table[attr] = {}.fromkeys(values, 0)
            self._heuristic_table[attr] = {}.fromkeys(values, None)
            for value in values:
                self._heuristic_table[attr][value] = (self._terms[attr][value].get_heuristic() / heuristic_accum)
        return

    def _get_prob_accum(self, antecedent):

        data = self._Dataset._original_data
        data_subset = data.loc[np.all(data[list(antecedent)] == pd.Series(antecedent), axis=1)]

        accum = 0
        for attr in self._attr_values.keys():
            if self._availability[attr]:
                if antecedent:
                    values = list(data_subset[attr].unique())
                else:
                    values = self._attr_values[attr]
                for value in values:
                    accum += self._heuristic_table[attr][value] * self._pheromone_table[attr][value]
        return accum

    def _get_pheromone_accum(self):

        accum = 0
        for attr, values in self._attr_values.items():
            for value in values:
                accum += self._pheromone_table[attr][value]
        return accum

    def _reset_availability(self):

        attrs = list(self._attr_values.keys())
        self._availability = {}.fromkeys(attrs, True)
        return

    def _get_probabilities(self, antecedent):

        prob_accum = self._get_prob_accum(antecedent)
        if prob_accum == 0:
            return None

        probabilities = []
        data = self._Dataset._original_data
        data_subset = data.loc[np.all(data[list(antecedent)] == pd.Series(antecedent), axis=1)]

        for attr in self._attr_values.keys():
            if self._availability[attr]:
                if antecedent:
                    values = list(data_subset[attr].unique())
                else:
                    values = self._attr_values[attr]

                for value in values:
                    prob = (self._heuristic_table[attr][value] * self._pheromone_table[attr][value]) / prob_accum
                    probabilities.append((prob, self._terms[attr][value]))
        return probabilities

    def size(self):
        return self._no_of_terms

    def available(self):

        for attr in self._availability:
            if self._availability[attr]:
                return True
        return False

    def sort_term(self, antecedent):

        probabilities = self._get_probabilities(antecedent)
        if not probabilities:
            return None

        nan_check = [math.isnan(p[0]) for p in probabilities]  # very low probabilities result in overflow - NaN values
        if any(nan_check): # !! resolve this problem better
            choice_idx = np.random.choice(len(probabilities), size=1)[0] # random with equal probs
        else:
            probs = [prob[0] for prob in probabilities]
            choice_idx = np.random.choice(len(probabilities), size=1, p=probs)[0]
        return probabilities[choice_idx][1]

    def update_availability(self, attr):
        self._availability[attr] = False
        return

    def get_cases(self, antecedent):

        all_cases = []
        for attr, value in antecedent.items():
            all_cases.append(self._terms[attr][value].covered_cases)

        cases = all_cases.pop()
        for case_set in all_cases:
            cases = list(set(cases) & set(case_set))
        return cases

    def pheromone_updating(self, antecedent, quality):

        # increasing pheromone of used terms
        for attr, value in antecedent.items():
            self._pheromone_table[attr][value] += self._pheromone_table[attr][value] * quality

        # Decreasing not used terms: normalization
        pheromone_normalization = self._get_pheromone_accum()
        for attr, values in self._attr_values.items():
            for value in values:
                self._pheromone_table[attr][value] = self._pheromone_table[attr][value] / pheromone_normalization

        self._reset_availability()
        return

    def pheromone_init(self):
        # Restores pheromone initial configuration
        initial_pheromone = 1 / self._no_of_terms
        self._pheromone_table = {}
        for attr, values in self._attr_values.items():
            self._pheromone_table[attr] = {}.fromkeys(values, initial_pheromone)
            self._count_table[attr] = {}.fromkeys(values, 0)
        return

    def att_discovered_terms(self, att_dic):
        # update terms use
        for attr, value in att_dic.items():
            self._logistic_count_table[attr][value] += 1
        return

    def heuristics_updating(self, dataset, weigh_score, offset):

        # calc heuristic values
        for v_dict in self._terms.values():
            for term in v_dict.values():
                term.set_heuristic(dataset)
                # with both description and cover attenuation (bellow)
                self._heuristic_table[term.attribute][term.value] = term.get_heuristic()\
                                                                    * self.dscrpt_attenuation(self._logistic_count_table[term.attribute][term.value], x0=offset)\
                                                                    * self.cover_attenuation(term, dataset.get_case_count(), weigh_score)
        accum = sum(list(map(lambda dic: sum(dic.values()), self._heuristic_table.values())))
        if accum == 0:
            return False

        # normalize heuristic values
        def att_heur(dic, tpl):
            dic[tpl[0]] = tpl[1] / accum
            return dic
        self._heuristic_table.update({key:
                                          reduce(lambda dic, tpl: att_heur(dic, tpl), value.items(), {})
                                      for key, value in self._heuristic_table.items()})
        return True

    @staticmethod
    def dscrpt_attenuation(x, L=1, k=1, x0=0):
        # inverse logistic function
        return 1 - (L / (1 + math.exp(-k * (x - x0))))

    @staticmethod
    def cover_attenuation(term, all_count, weigh_score):
        cases = term.covered_cases
        score = sum(list(map(lambda case: weigh_score ** all_count[case], cases))) / len(cases)
        return score

    def add_count(self, attr, v):
        self._count_table[attr][v] += 1
        return

    def get_pheromone_table(self):
        return copy.deepcopy(self._pheromone_table)

    def get_heuristic_table(self):
        return copy.deepcopy(self._heuristic_table)

    def get_counts_table(self):
        return copy.deepcopy(self._count_table)

    def get_logistic_table(self):
        return copy.deepcopy(self._logistic_count_table)

    def get_num_terms(self):
        return self._no_of_terms

    def get_attr_values(self):
        return self._attr_values.copy()

    def get_term(self, attr, val):
        return self._terms[attr][val]
