import numpy as np
import seaborn as sns
from memote.experimental.config import ExperimentConfiguration
import pandas as pd


def create_ge_confusion_matrix(model, exp_file_path):
    """
    Runs gene essentiality for a model and creates a confusion matrix plot

    Parameters
    ----------
    model : cobra.Model
    exp_file_path : str
        Path to experimental.yml file used by memote.

    Returns
    -------

    """
    experiment_data = ExperimentConfiguration(exp_file_path).load(model)
    cols = ['gene', 'essential']
    ge_simulated = experiment_data.essentiality['knockouts'].evaluate(model)[cols].copy()
    ge_experimental = experiment_data.essentiality['knockouts'].data[cols].copy()

    ge_simulated.set_index('gene', inplace=True)
    ge_experimental.set_index('gene', inplace=True)
    ge_experimental['predicted'] = ge_experimental['essential'].astype(bool)
    ge_simulated['actual'] = ge_simulated['essential'].astype(bool)
    del ge_simulated['essential']
    del ge_experimental['essential']
    merged = pd.concat([ge_experimental, ge_simulated], axis=1)

    tp = merged[merged.actual & merged.predicted]
    fp = merged[~merged.actual & merged.predicted]
    tn = merged[~merged.actual & ~merged.predicted]
    fn = merged[merged.actual & ~merged.predicted]

    conf_matrix = np.array(
        [tp.shape[0],  fn.shape[0], fp.shape[0], tn.shape[0]]
    ).reshape((2, 2))

    g = sns.heatmap(conf_matrix, annot=True, fmt='0.0f', cmap='Reds', cbar=False)
    g.set_xlabel("Actual Growth");
    g.set_ylabel("Predicted Growth");
    g.set_xticks([0.5, 1.5], ['True', 'False']);
    g.set_yticks([0.5, 1.5], ['True', 'False']);

    return {'TP': tp, 'FP': fp, 'FN': fn, 'TN': tn}