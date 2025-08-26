from rdkit import Chem
from rdkit.Chem import Descriptors


from typing import List

# TODO should I make another one for sdf or keep this?
def calc_all_descriptors_smi(l_mols: List[str]):

    desc_list = []

    for m in l_mols:

        mol = Chem.MolFromSmiles(m)


        # Calculate all descriptors
        all_desc = {}
        for name, func in Descriptors._descList:
            try:
                all_desc[name] = func(mol)
            except:
                all_desc[name] = None  # in case some descriptor fails

        desc_list.append(all_desc)

    return desc_list

def calc_all_descriptors_sdf(l_sdf: List[str]):

    desc_list = []

    for i,m in enumerate(l_sdf):

        mol = Chem.MolFromMolBlock(m)
        if mol is None:
            print(f"Block {i} failed")

        # Calculate all descriptors
        all_desc = {}
        for name, func in Descriptors._descList:
            try:
                all_desc[name] = func(mol)
            except:
                all_desc[name] = None  # in case some descriptor fails

        desc_list.append(all_desc)

    return desc_list
