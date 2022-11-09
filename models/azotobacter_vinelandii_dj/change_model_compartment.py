"""
Fixes iDTR1278.xml model compartment naming.

For some reason, the model compartment is labeled as '_u' for all species.
Within the model, the BIGG ids add _* to each ID, ie, atp_p for periplasm.
So with the _u error, there would be atp_p_u, which throws off comparisons.
Fixing this plus exporting model, then will do a comparison between the
models to see what the differences are.
"""
from memote.suite.cli.reports import diff
import cobra
import os

_file_path = os.path.dirname(__file__)
starting_model_f_name = 'iDT1278.xml'
s_model_path = os.path.join(_file_path, starting_model_f_name)

s_model = cobra.io.read_sbml_model(s_model_path)
s_model.id = "AV"


def modify_compartments(model):
    new_model = model.copy()
    for metabolite in new_model.metabolites:
        if metabolite.id[-2:] == "_u":
            # removes last two characters, namely "_u"
            metabolite.id = metabolite.id[:-2]
            # rename compartment
            metabolite.compartment = metabolite.id[-1]
    return new_model


av = modify_compartments(s_model)
output_model_name = 'azo_vine.xml'
output_model_path = os.path.join(_file_path, output_model_name)
cobra.io.write_sbml_model(av, output_model_path)
model_paths = [s_model_path, output_model_path]


if __name__ == '__main__':

    diff(
        [
            *model_paths,
            '--filename', os.path.join(_file_path, 'fix_compartment_diff.html')
         ]
    )
