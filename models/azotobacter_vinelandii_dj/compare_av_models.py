"""
Compare our model, updated iDTR1278.xml, with iAA1300

"""
from memote.suite.cli.reports import diff
import os
import logging

log = logging.getLogger()

_file_path = os.path.dirname(__file__)


if __name__ == '__main__':

    model_paths = ['azo_vine.xml', 'iAA1300.xml']
    diff(
        [
            *model_paths,
            '--experimental', 'data/experiments.yml',
            '--filename', os.path.join(_file_path, 'idt_vs_iAA_model_differences.html'),

         ]
    )
