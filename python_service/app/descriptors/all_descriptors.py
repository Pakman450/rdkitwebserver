from rdkit import Chem
from rdkit.Chem import Descriptors



def calc_all_descriptors(smiles: str):
    # Your molecule (example: aspirin)


    mol = Chem.MolFromSmiles(smiles)

    # Calculate all descriptors
    all_desc = {}
    for name, func in Descriptors._descList:
        try:
            all_desc[name] = func(mol)
        except:
            all_desc[name] = None  # in case some descriptor fails
    
    return all_desc
