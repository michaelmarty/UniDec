import glob, os
import comtypes.client as w32
from rdkit import Chem
import pandas as pd


def cdx2smi(cdx_folder, dest_folder):
    files = glob.glob(cdx_folder + '/*.cdx') + glob.glob(cdx_folder + '/*.cdxml')
    ChemDraw = w32.CreateObject("ChemDraw.Application")

    smi_dict = {'file': [], 'SMILES': []}
    for tmp_file in files:
        file_name = os.path.basename(tmp_file)
        try:
            tmp_compound = ChemDraw.Documents.Open(tmp_file)
            smiles = tmp_compound.Objects.Data("chemical/x-smiles")
            mol = Chem.MolFromSmiles(smiles)
            smiles = Chem.MolToSmiles(mol)

        except:
            smiles = ''
        smi_dict['file'].append(file_name)
        smi_dict['SMILES'].append(smiles)
    ChemDraw.quit()
    df = pd.DataFrame.from_dict(smi_dict)
    df.to_csv(dest_folder + '/smile.csv')


if __name__ == '__main__':
    # cdx_folder = input('please input the path which ".cdx" saved in (eg. C:/tmp/): ')
    # dest_folder = input('please input one existed folder to save smiles file (eg. C:/tmp/smiles/): ')
    cdx_folder = r"C:\Data\ChemDrawStds"
    dest_folder = r"C:\Data\ChemDrawStds\smi"
    cdx2smi(cdx_folder, dest_folder)

    ### usage "python cdx2smi"