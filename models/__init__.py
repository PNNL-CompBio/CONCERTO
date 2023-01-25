import cobra
import pathlib

path = pathlib.Path(__file__).parent

av = cobra.io.read_sbml_model(
    path.joinpath('azotobacter_vinelandii_dj', 'azo_vine.xml').__str__()
)

syn = cobra.io.read_sbml_model(
    path.joinpath('synechococcus_elongatus_pcc_7942', 'iJB785.xml').__str__()
)
rt = cobra.io.read_sbml_model(
    path.joinpath('rhodosporidium_toruloides_ifo_08804', 'Rhodo_Toru.xml').__str__()
)
