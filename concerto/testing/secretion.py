import cobra
import pandas as pd
import pytest
from memote.utils import annotate, wrapper, truncate
import numpy as np

def get_objective_name(model: cobra.Model):
    js = model.objective.to_json()
    for arg in js['expression']['args']:
        if arg['args'][0]['value'] == 1.0:
            return arg['args'][1]['name']


def check_exchanges(metabolites_to_secret, media, model):
    opt_each = dict()
    bm_name = get_objective_name(model)
    lower_bound = model.optimize().objective_value * 0.5
    exchanges = [i.id for i in model.exchanges]
    for i in metabolites_to_secret:
        m = model.copy()
        m.reactions.get_by_id(bm_name).lower_bound = lower_bound
        ex = f'EX_{i}_e'
        if ex in exchanges:
            m.objective = ex
            m.medium = media
            opt_each[i] = m.optimize().objective_value
        else:
            opt_each[i] = 0
    optimize_experimental = pd.Series(opt_each).sort_values()
    return optimize_experimental


def get_excreted_metabolites(model, expected_excreted_list):
    excreted_results = check_exchanges(expected_excreted_list, model.medium, model)
    tp = excreted_results.loc[np.abs(excreted_results) > 1e-5]
    fn = excreted_results.loc[np.abs(excreted_results) < 1e-5]
    return list(tp.index.values), list(fn.index.values)