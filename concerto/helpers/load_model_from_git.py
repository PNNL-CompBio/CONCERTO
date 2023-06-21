from cobra.io import read_sbml_model
import requests


def load_model_from_git(model_name: str):
    """Helper function for loading a SBML compliant model into CobraPy directly from Git

    Args:
        model_name (str): The name of the model you wish to load

    Returns:
        model (cobra.core.model.Model): A Cobra Model object

    """
    
    model_url_dict = dict([('Rhodosporidium','https://raw.githubusercontent.com/PNNL-CompBio/RToruGEM/main/rtoru/Rhodo_Toru.xml'),
                           ('Azotobacter','https://raw.githubusercontent.com/PNNL-CompBio/iAzotobacterVinelandiiGEM/main/a_vine/azo_vine.xml'),
                           ('Synechococcus','https://raw.githubusercontent.com/PNNL-CompBio/S-elongatus7942/main/syn_elong/syn_elong.xml')])

    xml = requests.get(model_url_dict[model_name])

    with open('tmp.xml', 'wb') as f:
        f.write(xml.content)
        model = read_sbml_model('tmp.xml')
        
    return model 
