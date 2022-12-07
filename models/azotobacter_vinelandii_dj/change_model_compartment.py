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
import logging
from models.helpers.id_mappers import update_ids, add_annotations

log = logging.getLogger()

_file_path = os.path.dirname(__file__)
starting_model_f_name = 'iDT1278.xml'
s_model_path = os.path.join(_file_path, starting_model_f_name)

model = cobra.io.read_sbml_model(s_model_path)
model.id = "AV"

output_model_name = 'azo_vine.xml'
output_model_path = os.path.join(_file_path, output_model_name)


def write_model():
    cobra.io.write_sbml_model(model, output_model_path)


def update_1():
    # updates bug in compartment of the model
    log.info("Updating compartments")
    for metabolite in model.metabolites:
        if metabolite.id[-2:] == "_u":
            # removes last two characters, namely "_u"
            metabolite.id = metabolite.id[:-2]
            # rename compartment
            metabolite.compartment = metabolite.id[-1]


def update_2():
    log.info('Adding annotations to metabolites')

    add_annotations(model)


def update_model():
    # Fix compartments
    update_1()
    update_2()
    write_model()


if __name__ == '__main__':
    update_model()
    model_paths = [s_model_path,  'azo_vine2.xml', output_model_path]
    diff(
        [
            *model_paths,
            '--filename', os.path.join(_file_path, 'model_differences.html')
         ]
    )
