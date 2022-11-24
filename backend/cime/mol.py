import base64
import logging
import re
from io import BytesIO

from rdkit import Chem
from rdkit.Chem import Draw, rdFMCS
from rdkit.Chem.Draw import SimilarityMaps


def get_mcs(mol_list):

    if len(mol_list) <= 1:
        return Chem.MolFromSmiles("*")

    if type(mol_list[0]) == str:
        # TODO: handle invalid smiles
        mol_list = [Chem.MolFromSmiles(sm) for sm in mol_list]

    # completeRingsOnly=True # there are different settings possible here
    res = rdFMCS.FindMCS(mol_list, timeout=60, matchValences=False, ringMatchesRingOnly=True, completeRingsOnly=True)
    if res.canceled:
        patt = Chem.MolFromSmiles("*")
    else:
        patt = res.queryMol

    return patt


def smiles_to_base64(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m:
        return mol_to_base64(m)
    else:
        return "invalid smiles"


def mol_to_base64(m):
    pil_img = Draw.MolToImage(m)

    buffered = BytesIO()
    pil_img.save(buffered, format="JPEG")
    img_str = base64.b64encode(buffered.getvalue())
    buffered.close()
    return img_str.decode("utf-8")


def mol_to_base64_highlight_substructure(mol, patt, width=250, d=None, showMCS=True):
    width = int(width)

    if d is None:
        d = Chem.Draw.rdMolDraw2D.MolDraw2DCairo(width, width)  # MolDraw2DSVG
    if showMCS:
        hit_ats = list(mol.GetSubstructMatch(patt))
        hit_bonds = []
        for bond in patt.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(mol.GetBondBetweenAtoms(aid1, aid2).GetIdx())

        # specify black color for each atom and bond index
        col = (0, 0, 0, 0.1)
        atom_cols = {}
        for i, at in enumerate(hit_ats):
            atom_cols[at] = col
        bond_cols = {}
        for i, bd in enumerate(hit_bonds):
            bond_cols[bd] = col

        Chem.Draw.rdMolDraw2D.PrepareAndDrawMolecule(
            d, mol, highlightAtoms=hit_ats, highlightBonds=hit_bonds, highlightAtomColors=atom_cols, highlightBondColors=bond_cols
        )

    d.FinishDrawing()
    stream = BytesIO(d.GetDrawingText())
    # image = Image.open(stream).convert("RGBA")
    img_str = base64.b64encode(stream.getvalue())
    stream.close()
    return img_str.decode("utf-8")  # d.GetDrawingText()