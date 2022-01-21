# %%
#####################################
### import the used functions
#####################################

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors
from rdkit.Chem import PandasTools
import pandas as pd
import catboost
import shap
import numpy as np



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


# define a function to calculate count based morgan fingerprints
def calc_morgan_counts_for_smiles(mol, radius=2, nBits=2048):
    try:
        fpCounts = AllChem.GetHashedMorganFingerprint(
            mol, radius, nBits=nBits, includeRedundantEnvironments=False)
        array = np.zeros((0,), dtype=np.int8)
        DataStructs.ConvertToNumpyArray(fpCounts, array)
    except:
        print("Failure in calc_morgan_counts_for_smiles()")
        array = np.repeat(np.nan, nBits)

    return array


def getSubstructIndex(mol,atomID,radius):
    if radius>0:
        env = Chem.FindAtomEnvironmentOfRadiusN(mol,radius,atomID)
        atomsToUse=[]
        for b in env:
            atomsToUse.append(mol.GetBondWithIdx(b).GetBeginAtomIdx())
            atomsToUse.append(mol.GetBondWithIdx(b).GetEndAtomIdx())
        atomsToUse = list(set(atomsToUse))
    else:
        atomsToUse = [atomID]
        env=None
    return(atomsToUse)


def calc_indices_for_fingerprints(smiles,radius = 2,nBits=2048):
   mol = Chem.MolFromSmiles(smiles,sanitize = True)
   bi = {}
   fp = AllChem.GetHashedMorganFingerprint(mol,radius = radius, nBits=nBits, includeRedundantEnvironments=False,bitInfo=bi)
   fp = fp.GetNonzeroElements()
   #weights = [0  for i in mol.GetAtoms()]
   #weights[slice(mol.GetSubstructMatch(Chem.MolFromSmiles('CNC(=O)C(C)(C)C',sanitize = True)))] += 1


   fp_indices = {}
   for fragment_id, count in fp.items():
      sub_indices = []
      for i in range(count):
         #print(i)
         root, radius = bi[fragment_id][i]
         subindex = getSubstructIndex(mol,root,radius)
         sub_indices = sub_indices + subindex

      sub_indices = list(set(sub_indices))

      fp_indices[fragment_id] = sub_indices
   fp_indices = {"fp_morgan_counts_" +  str(k): v for k, v in fp_indices.items() }

   return(fp_indices)


def atom_attributions(mol, indices, shap_vals, featureNames, direction="both"):
    """
    Sum up the attributions to the model predictions (Shap values) for each atom in mol.

    Parameters:
        mol: rdkit molecule object
        indices (dict): a dictionary with keys = feature names and values = list of atom positions in mol (output of calc_indices_for_fingerprints())
        shap_vals: a numpy array of dimensions (len(featureNames))
        featureNames (np.array): numpy array containing the feature names for shap_vals
        direction (str): whether to consider positive (up), negative (down) or both shap_vals contributions

    Returns:
        np.array: each element = 1 atom of mol and the values are the summed contributions from shap_vals
    """
    ## TO DO: add support for MACCS keys

    # get for each fingerprint bit the atom numbers mapping to that bit and initialise weights to 0
    curWeights = np.array([0 for i in range(mol.GetNumAtoms())], 'float64')

    # step through each index (dictionary)
    for j in indices:
        # retrieve index to subset into the featureNames and the corresponding shap values array
        curFPIndex = [i for i, s in enumerate(featureNames) if s == j]
        if len(curFPIndex) == 0:
            continue
        if len(curFPIndex) > 1:
            raise IndexError("Found several matches for feature " + j)
        # do something only if (i) the fingerprint bit was in the selected features
        # and (ii) only if the retrieved shap value is positive/negative (direction argument)
        if (direction == "up" and shap_vals[curFPIndex[0]] > 0) or (direction == "down" and shap_vals[curFPIndex[0]] < 0) or (direction == "both"):
            # use the atom numbers ('positions') to update the curWeights array, where
            # each atom of the current molecule is represented
            # one atom may be part of several fingerprints so the summing has to be done across
            # all fingerprints
            for pos in indices[j]:
                #value = (curWeights[pos] + shap_vals[curFPIndex[0]])
                value = (curWeights[pos] + (shap_vals[curFPIndex[0]]/len(indices[j])))
                curWeights[pos] = value

    return curWeights

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



################################################################################
##calculate the logp for each of them
################################################################################

logp =  [Descriptors.MolLogP(i) for i in mols]


################################################################################
##calculate the morgan fingerprints
################################################################################


# calculate the morgan counts fingerprints
featuresMorganCounts = [calc_morgan_counts_for_smiles(
    mol, radius=2, nBits=2048) for mol in mols]
columns = ["fp_morgan_counts_" + str(i) for i in range(2048)]
featuresMorganCounts = pd.DataFrame(featuresMorganCounts, columns=[
                                  "fp_morgan_counts_" + str(i) for i in range(2048)])

##################################################
##train a model on the values
##################################################

model = catboost.CatBoostRegressor()
# Fit model
model.fit(featuresMorganCounts, logp)

##################################################
##calculate the shap values
##################################################

explainer = shap.TreeExplainer(model)
shap_vals = explainer.shap_values(featuresMorganCounts)


####################################################
#now assign weights based on the model
###################################################
indices =  [calc_indices_for_fingerprints(smi,radius = 2,nBits=2048) for smi in smiles]
#indices = indices[1]
featureNames = featuresMorganCounts.columns
#shap_vals = shap_vals[1]
#mol = mols[1]

modelAttributions = [list(atom_attributions(mols[i], indices[i], shap_vals[i],featureNames , direction="down")) for i in range(len(smiles))]


# %%
###################################################
##Calculate the Gasteiger Charges
###################################################

for l in mols:
    AllChem.ComputeGasteigerCharges(l)
charges_contribs = [
    [mols[i].GetAtomWithIdx(j).GetDoubleProp('_GasteigerCharge') for j in range(mols[i].GetNumAtoms())] for i in
    range(len(mols))]


##################################################
# put the smiles and info into a data frame
###################################################
sdfExportData = pd.DataFrame({'names': smiles_names,
        'smiles': smiles})

####################################################################################
##Add the partial charges to the sdfExportDataFrame as one string for each compound
####################################################################################

sdfExportData["atom.Contribution.rep_modelAttributions"] = [' '.join([str(round(i, 4)) for i in modelAttributions[j]]) for
                                                           j in range(len(modelAttributions))]

####################################################################################
##Add the model attributions for logp to the data frame
####################################################################################
sdfExportData["atom.Contribution.rep_GasteigerCharges"] = [' '.join([str(round(i, 4)) for i in charges_contribs[j]]) for
                                                           j in range(len(charges_contribs))]

####################################################


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





