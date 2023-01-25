from memote.suite.cli.reports import diff
import cobra
import os
import logging
from cobra.manipulation.modify import rename_genes
log = logging.getLogger()

_file_path = os.path.dirname(__file__)
starting_model_f_name = 'Rt_IFO0880.xml'
s_model_path = os.path.join(_file_path, starting_model_f_name)

model = cobra.io.read_sbml_model(s_model_path)
model.id = "RT"

output_model_name = 'Rhodo_Toru.xml'
output_model_path = os.path.join(_file_path, output_model_name)


def write_model():
    cobra.io.write_sbml_model(model, output_model_path)


def update_1():
    # updates bug in compartment of the model
    log.info("Adding RT to prefix")

    new_names = {}
    for gene in model.genes:
        new_names[gene.id] = f'RT_{gene.id}'
    rename_genes(model, new_names)


def update_model():
    # Fix compartments
    update_1()
    write_model()


if __name__ == '__main__':
    update_model()
    model_paths = [s_model_path, output_model_path]
    diff(
        [
            *model_paths,
            '--filename', os.path.join(_file_path, 'model_differences.html')
        ]
    )
