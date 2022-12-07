# from magine.mappings.chemical_mapper import ChemicalMapper
from ref_ids.id_mapper import final_mapper


# cm = ChemicalMapper()
# new_dict = cm._to_dict('chemical_formula', 'inchikey')
# name_to_inchi = cm._to_dict('name', 'inchikey')
# name_to_accession = cm._to_dict('name', 'accession')
# name_to_kegg = cm._to_dict('name', 'kegg_id')
# name_to_pubchem = cm._to_dict('name', 'pubchem_compound_id')
# name_to_chebi_id = cm._to_dict('name', 'chebi_id')
# name_to_biocyc_id = cm._to_dict('name', 'biocyc_id')

valid_ids = [
    'kegg_id', 'name', 'accession', 'chebi_id', 'inchikey',
    'chemspider_id', 'biocyc_id', 'synonyms', 'iupac_name',
    'pubchem_compound_id', 'protein_associations',
    'ontology', 'drugbank_id', 'chemical_formula',
    'smiles', 'metlin_id'
]


def _add_from_dict(metabolite,  new_id_name, ref_dict):
    chem_name = metabolite.name

    if chem_name in ref_dict:
        id_to_add = ref_dict[chem_name]
        if len(id_to_add) == 1:
            metabolite.annotation[new_id_name] = list(id_to_add)[0]
        else:
            metabolite.annotation[new_id_name] = '|'.join(sorted(id_to_add))


def update_ids(model):

    for met in model.metabolites:
        _add_from_dict(met, 'inchikey', name_to_inchi)
        _add_from_dict(met, 'hmdb', name_to_accession)
        _add_from_dict(met, 'kegg.compound', name_to_kegg)
        _add_from_dict(met, 'pubchem.compound', name_to_pubchem)
        _add_from_dict(met, 'chebi', name_to_chebi_id)
        _add_from_dict(met, 'biocyc', name_to_biocyc_id)


def add_annotations(model):
    # add annotations to model
    for meta in model.metabolites:
        m_id = meta.id
        if m_id in final_mapper:
            for ref_id in final_mapper[m_id]:
                meta.annotation[ref_id] = final_mapper[m_id][ref_id]
