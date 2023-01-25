from memote.suite.cli.reports import snapshot


if __name__ == '__main__':

    args = [
        '../Rhodo_Toru.xml',
        '--filename', 'essential_only_report.html',

        '--experimental', '../data/experiments.yml',
        '--exclusive', 'test_essentiality',
        '--pytest-args', '--tb=long'


    ]

    snapshot(
       args
    )
