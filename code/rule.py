import pandas as pd
import statsmodels.api as sm
from lifelines import KaplanMeierFitter

from utils import NoIndent


class Rule:

    def __init__(self, dataset, comp):
        self.antecedent = {}
        self.sub_group_cases = dataset.get_instances()
        self.no_covered_cases = len(self.sub_group_cases)
        self.fitness = 0.0
        self.p_value = None
        self._string_repr = ()
        self._KMmodel = {'subgroup': KaplanMeierFitter(), 'complement': KaplanMeierFitter()}
        self._complement_cases = []
        self._Dataset = dataset
        self._comparison = comp

    @property
    def size(self):
        return len(self.antecedent)

    @property
    def mean_survival(self):
        return self._Dataset.survival_times[1].iloc[self.sub_group_cases].mean()

    @property
    def complement_mean_survival(self):
        return self._Dataset.survival_times[1].iloc[self._complement_cases].mean()

    @property
    def id(self):
        return self._string_repr[0]

    @property
    def description(self):
        return self._string_repr[1]

    def is_exceptional(self, rule, alpha):
        # calculate p-value btw (new_rule, rule)
        times = self._Dataset.survival_times[1][self.sub_group_cases].to_list() + self._Dataset.survival_times[1][rule.sub_group_cases].to_list()
        events = self._Dataset.events[1][self.sub_group_cases].to_list() + self._Dataset.events[1][rule.sub_group_cases].to_list()
        group_id = ['self']*self._Dataset.survival_times[1][self.sub_group_cases].shape[0] + ['rule']*self._Dataset.survival_times[1][rule.sub_group_cases].shape[0]
        try:
            _, p_val = sm.duration.survdiff(time=times, status=events, group=group_id)
        except: # equal rules/ rules with no event
            p_val = 1

        if p_val >= alpha:
            return False
        else:
            return True

    def is_in(self, rule):
        attr_intersec = self.get_attributes().intersection(rule.get_attributes())
        if attr_intersec == rule.get_attributes():
            if self.get_terms(attr_intersec).issubset(rule.get_terms()):
                    return True
        return False

    def has_root(self, rule):
        if self.get_terms().intersection(rule.get_terms()):
            return True
        else:
            return False

    def has_merge(self, rule):
        if self.get_attributes() == rule.get_attributes():
            return True
        else:
            return False

    def _construct_from_items(self, items, terms_mng):
        antecedent = {}
        cases = set()
        for (attr,val) in items:
            if attr not in antecedent:
                antecedent[attr] = {val}
            else:
                antecedent[attr] = set(antecedent[attr]).union({val})
            cases = cases.union(set(terms_mng.get_term(attr, val).covered_cases))

        self.antecedent = antecedent.copy()
        self.set_cases(list(cases))
        self.set_fitness()
        return

    def root(self, rule1, rule2, terms_mng):
        items = rule1.get_terms().intersection(rule2.get_terms())
        return self._construct_from_items(items, terms_mng)

    def merge(self, rule1, rule2, terms_mng):
        items = rule1.get_terms().union(rule2.get_terms())
        return self._construct_from_items(items, terms_mng)

    def get_attributes(self):
        return set([term[0] for term in self.get_terms()])

    def get_terms(self, attrs=None):
        terms = []
        if not attrs:
            attrs = self.antecedent.keys()

        for attr in attrs:
            if isinstance(self.antecedent[attr], list) or isinstance(self.antecedent[attr], set):
                for v in self.antecedent[attr]:
                    terms.append((attr, v))
            else:
                terms.append((attr, self.antecedent[attr]))
        return set(terms)

    def set_cases(self, cases):
        self.sub_group_cases = cases
        self.no_covered_cases = len(cases)
        self._complement_cases = list(set(self.sub_group_cases) ^ set(self._Dataset.get_instances()))
        return

    def set_fitness(self):
        # against population
        if self._comparison == 'population':
            times = self._Dataset.survival_times[1][self.sub_group_cases].to_list() + self._Dataset.survival_times[1].to_list()
            events = self._Dataset.events[1][self.sub_group_cases].to_list() + self._Dataset.events[1].to_list()
            group_id = ['sg'] * self._Dataset.survival_times[1][self.sub_group_cases].shape[0] + ['pop'] * self._Dataset.survival_times[1].shape[0]
            try:
                _, self.p_value = sm.duration.survdiff(time=times, status=events, group=group_id)
            except:
                print("!! Raise < sm.duration.survdiff > except rule-fitness:")
                print('...baseline: population')
                print('rule: {}'.format(self.antecedent))
                self.p_value = 1
        # against complement
        if self._comparison == 'complement':
            sg = pd.Series('sub_group', index=self.sub_group_cases)
            cpm = pd.Series('complement', index=self._complement_cases)
            group = pd.concat([sg, cpm], axis=0, ignore_index=False).sort_index()
            try:
                _, self.p_value = sm.duration.survdiff(self._Dataset.survival_times[1], self._Dataset.events[1], group)
            except:
                print("!! Raise < sm.duration.survdiff > except rule-fitness:")
                print('...baseline: complement')
                print('rule: {}'.format(self.antecedent))
                self.p_value = 1
        self.fitness = 1 - self.p_value
        return

    def construct(self, terms_mgr, min_case_per_rule):
        # ANTECEDENT CONSTRUCTION
        while terms_mgr.available():

            term = terms_mgr.sort_term(self.antecedent)
            if not term:
                break

            covered_cases = list(set(term.covered_cases) & set(self.sub_group_cases))

            # if len(covered_cases) > 1: == version with subgroup size quality component
            if len(covered_cases) >= min_case_per_rule:
                self.antecedent[term.attribute] = term.value
                self.set_cases(covered_cases)
                terms_mgr.update_availability(term.attribute)
                terms_mgr.add_count(term.attribute, term.value)
            else:
                break
        self.set_fitness()
        return

    def equals(self, comp_rule):
        return self.get_terms() == comp_rule.get_terms()

    def set_string_repr(self, index):
        rule_id = 'R' + str(index)

        for (key, value) in self.antecedent.items():
            if not isinstance(value, set):
                self.antecedent[key] = {value}

        string = ' & '.join(['({}={})'.format(key, value) for (key, value) in self.antecedent.items()])
        self._string_repr = (rule_id, string)
        return

    def get_full_description(self, rule_comp=None):

        if rule_comp:
            c_intersec = set(self.sub_group_cases).intersection(set(rule_comp.sub_group_cases))
            c_union = set(self.sub_group_cases).union(set(rule_comp.sub_group_cases))
            jaccard_cover = len(c_intersec)/len(c_union)

            d_intersec = self.get_terms().intersection(rule_comp.get_terms())
            d_union = self.get_terms().union(rule_comp.get_terms())
            jaccard_descript = len(d_intersec) / len(d_union)

            info = '[size={}/{}; jaccard-c={:.2f}; jaccard-d={:.2f}]'.format(self.no_covered_cases,
                                                                             self._Dataset.size,
                                                                             jaccard_cover,
                                                                             jaccard_descript)
        else:
            info = '[size={}/{}]'.format(self.no_covered_cases,self._Dataset.size)
        return self.id, self.description, info

    def set_KMmodel(self, alpha):
        # Subgroup and Complement of induced rule
        sub_group_times = self._Dataset.survival_times[1].iloc[self.sub_group_cases]
        sub_group_events = self._Dataset.events[1].iloc[self.sub_group_cases]
        complement_times = self._Dataset.survival_times[1].iloc[self._complement_cases]
        complement_events = self._Dataset.events[1].iloc[self._complement_cases]

        self._KMmodel['subgroup'].fit(sub_group_times, sub_group_events, label='estimates', alpha=alpha)
        self._KMmodel['complement'].fit(complement_times, complement_events, label='KM estimates for complement', alpha=alpha)
        return

    def get_result(self):
        for attr, val in self.antecedent.items():
            if isinstance(val, set):
                self.antecedent[attr] = list(val)
        dic = {
            'string': self._string_repr[1],
            'antecedent': NoIndent(self.antecedent),
            'num_cases': self.no_covered_cases,
            'fitness': self.fitness,
            'pvalue': self.p_value,
            'cases': NoIndent(self.sub_group_cases),
            'survival_fnc': NoIndent(self._KMmodel['subgroup'].survival_function_.reset_index().to_dict(orient='records'))
        }
        return self._string_repr[0], dic
