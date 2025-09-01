from rdkit import Chem
from rdkit.Chem import Descriptors


from typing import List

def calc_all_descriptors_smi(l_mols: List[str]):

    desc_list = []

    for m in l_mols:

        mol = Chem.MolFromSmiles(m)

        inchi_str = Chem.inchi.MolToInchi(mol)

        # Calculate all descriptors
        all_desc = {}

        # place in the smiles that was used by user
        all_desc['smiles'] = m
        all_desc['inchi_str'] = inchi_str

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

        # NOTE: this codeblock print out the iteratira
        # number that is in this l_sdf. so 
        # if you set chunk num to 1.
        # you will get index num 0,0,0,0,0
        if mol is None:
            print(f"Block {i} failed")

        smi = Chem.MolToSmiles(mol)

        inchi_str = Chem.inchi.MolToInchi(mol)

        # Calculate all descriptors
        all_desc = {}
        all_desc["smiles"] = smi
        all_desc['inchi_str'] = inchi_str

        for name, func in Descriptors._descList:
            try:
                all_desc[name] = func(mol)
            except:
                all_desc[name] = None  # in case some descriptor fails

        desc_list.append(all_desc)

    return desc_list
