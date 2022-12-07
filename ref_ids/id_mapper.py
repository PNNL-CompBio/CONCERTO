import pandas as pd
import os

_path = os.path.dirname(__file__)

id_map = pd.read_csv(
    os.path.join(_path, 'bigg_models_metabolites.txt'),
    delimiter='\t'
)


_remap_ids = {
    'KEGG_Compound': 'kegg.compound',
    'KEGG_Drug': 'KEGG_Drug',
    'KEGG_Glycan': 'KEGG_Glycan',
    'CHEBI': 'chebi',
    'Human_Metabolome_Database': 'hmdb',
    'BioCyc': 'biocyc',
    'MetaNetX_(MNX)_Chemical': 'metanetx.chemical',
    'InChI_Key': 'inchikey',
    'SEED_Compound': 'seed.compound',
    'Reactome_Compound': 'reactome',
    'LipidMaps': 'LipidMaps',

}
final_mapper = dict()
count = 0
id_map['database_links'] = id_map['database_links'].str.split(';')
for i, j in id_map[['bigg_id', 'database_links']].values:
    local_map = dict()
    if isinstance(j, float):
        continue
    for ref_id in j:
        key, value = ref_id.split(': ')
        value = value.lstrip().rstrip().rsplit('/')[-1]
        key = key.lstrip().rstrip().replace(' ', '_')
        key = _remap_ids[key]
        local_map[key] = value
    final_mapper[i] = local_map
