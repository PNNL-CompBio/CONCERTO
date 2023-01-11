from memote.suite.cli.reports import snapshot


if __name__ == '__main__':

    args = [
        '../azo_vine.xml',
        '--filename', 'growth_report.html',

        '--experimental', '../data/experiments.yml',
        '--exclusive', 'test_growth',
        '--pytest-args', '--tb=long'


    ]

    snapshot(
       args
    )
