import pandas as pd
import pathlib

_path = pathlib.Path(__file__).parent
_f_path = _path.joinpath('plate_to_bigg.csv').__str__()
biolog_map = pd.read_csv(_f_path, index_col=False)
