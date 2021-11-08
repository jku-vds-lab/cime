# %%
#####################################
### import the used functions
#####################################

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools
import pandas as pd

#######################################
###define the used functions
#######################################

def smiles_to_mols(smis):
    """Convert a list of smiles to a list of RDKit Mol objects"""
    if not isinstance(smis, list):
        raise TypeError("Please provide smiles as a list")
    mols = []
    successful_inds = []
    for ind, smi in enumerate(smis):
        # RDKit will sometimes 'fail' but instead of throwing an
        # error it will (sometimes!) print a warning and return None.
        m = Chem.MolFromSmiles(smi, sanitize=True)
        if m is None:
            print("Mol generation failed for", smi, "(index", ind, ")\n")
        else:
            mols.append(m)
            successful_inds.append(ind)
    return mols, successful_inds


#######################################################
##define some smiles
#######################################################

smiles = ['Nc1nc(NC2CC2)c2ncn(C3C=CC(CO)C3)c2n1',
'CC(=O)NCCCS(=O)(=O)O',
'CCCC(=O)Nc1ccc(OCC(O)CNC(C)C)c(C(C)=O)c1',
'CC(=O)Nc1ccc(O)cc1',
'CC(=O)Nc1nnc(S(N)(=O)=O)s1',
'CC(=O)NO'
]

################################################################################
##make sure that the smiles are canonical to have full compatibility with CIME
################################################################################

#canonicalize smiles
mols, successful_inds = smiles_to_mols(smiles)
# get canonical smiles
smiles = [Chem.MolToSmiles(m) for m in mols]

#generate some names for the smiles
smiles_names = [i for i in range(len(smiles))]

##################################################
# put the smiles and info into a data frame
###################################################
sdfExportData = pd.DataFrame({'names': smiles_names,
        'smiles': smiles})


# %%
###################################################
##Calculate the Gasteiger Charges
###################################################

for l in mols:
    AllChem.ComputeGasteigerCharges(l)
charges_contribs = [
    [mols[i].GetAtomWithIdx(j).GetDoubleProp('_GasteigerCharge') for j in range(mols[i].GetNumAtoms())] for i in
    range(len(mols))]


####################################################################################
##Add the partial charges to the sdfExportDataFrame as one string for each compound
####################################################################################

sdfExportData["atom.Contribution.rep_GasteigerCharges"] = [' '.join([str(round(i, 4)) for i in charges_contribs[j]]) for
                                                           j in range(len(charges_contribs))]

# %%
#######################################################################################
#add the chemical info to the data frame
#######################################################################################
# smilesCol: in which column the smiles are located
PandasTools.AddMoleculeColumnToFrame(sdfExportData, smilesCol='smiles', molCol='ROMol', includeFingerprints=False)
# to avoid duplicate smiles column
del sdfExportData['smiles']


# %%
########################################################################
###Export as an sdf file you can read into CIME
########################################################################

PandasTools.WriteSDF(sdfExportData,"SDF_OUTPUT.sdf",
                     properties=list(sdfExportData.columns), idName="RowID")





