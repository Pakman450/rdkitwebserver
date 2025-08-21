from rdkit import Chem
from rdkit.Chem import Descriptors


from typing import List


def calc_all_descriptors(smiles: List[str]):
    # Your molecule (example: aspirin)

    desc_list = []


    for smi in smiles:

        mol = Chem.MolFromSmiles(smi)

        # Calculate all descriptors
        all_desc = {}
        for name, func in Descriptors._descList:
            try:
                all_desc[name] = func(mol)
            except:
                all_desc[name] = None  # in case some descriptor fails

        desc_list.append(all_desc)

    return desc_list
