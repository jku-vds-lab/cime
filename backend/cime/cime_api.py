import pickle
from flask import Blueprint, request, current_app, abort, jsonify
import copy
from rdkit import Chem
from rdkit.Chem import TemplateAlign, AllChem, Draw
from rdkit.Chem.PropertyMol import PropertyMol
from io import BytesIO, StringIO
import base64
import pandas as pd
import time
import logging
from .mol import mol_to_base64_highlight_importances, get_mcs, mol_to_base64_highlight_substructure, smiles_to_base64

_log = logging.getLogger(__name__)

cime_api = Blueprint('cime', __name__)


fingerprint_modifier = "fingerprint"
descriptor_names_no_lineup = [fingerprint_modifier, "rep"]
descriptor_names_show_lineup = ["pred", "predicted", "measured"]
smiles_prefix = "smiles"
smiles_col = 'SMILES'
mol_col = "Molecule"


def get_cime_dbo():
    return current_app.config['CIME_DBO']


def tryParseFloat(value):
    try:
        return float(value)
    except ValueError:
        return value


def sdf_to_df(file, id_col_name='ID', smiles_col_name="SMILES", mol_col_name='Molecule'):
    # Inspired by https://github.com/rdkit/rdkit/blob/master/rdkit/Chem/PandasTools.py#L452
    # adapt LoadSDF method from rdkit. this version creates a pandas dataframe excluding atom specific properties and does not give every property the object type
    _log.info('Starting processing of SDF file')
    records = []
    rep_list = set()
    # TODO: Check if parallelization is possible?
    with Chem.ForwardSDMolSupplier(file) as suppl:
        i = 0
        for mol in suppl:
            if mol is None:
                continue

            # Define the row to be filled
            row = {}

            # Add prop names
            for k in mol.GetPropNames():
                if k.startswith("atom"):
                    rep_list.add(k)
                else:
                    row[k] = tryParseFloat(mol.GetProp(k))

            # Add name
            if mol.HasProp('_Name'):
                row[id_col_name] = mol.GetProp('_Name')

            # Add Mol from smiles
            if smiles_col_name is not None:
                try:
                    row[smiles_col_name] = Chem.MolToSmiles(mol)
                except:
                    row[smiles_col_name] = None

            # Add property mol
            if mol_col_name is not None:
                # TODO: we only use that once in the code to create the fingerprints/create images, that could be done from the smiles as well?
                row[mol_col_name] = pickle.dumps(PropertyMol(mol), protocol=4)
                pass

            records.append(row)

            i += 1
            if i % 1000 == 0:
                _log.info(f'Processed {i} entries')

    _log.info(f'Finished processing of {i} entries')
    frame = pd.DataFrame(records)
    return frame, list(rep_list)


@cime_api.route('/get_uploaded_files_list', methods=['GET'])
def get_uploaded_files_list():
    # TODO: Check permissions will probably be done in dbo implementation
    return jsonify([{'id': d.id, 'name': d.name} for d in get_cime_dbo().get_datasets()])


@ cime_api.route('/delete_file/<id>', methods=['GET'])
def delete_uploaded_file(id):
    # TODO: Is this required?
    if id == "test.sdf" or id == "test":
        pass
        # return {"deleted": "false", "error": "can't delete default file"}

    # TODO: Rewrite to use id
    deleted = get_cime_dbo().delete_dataset_by(id=id)

    # TODO: Refactor to true and false instead of strings
    return {"deleted": "true" if deleted else "false"}


# TODO
# @cime_api.route('/get_atom_rep_list', methods=["GET"])
@ cime_api.route('/get_atom_rep_list/<id>', methods=["GET"])
def get_atom_rep_list(id):
    dataset = get_cime_dbo().get_dataset_by(id=id)

    if not dataset or not dataset.rep_list:
        abort(404)

    if dataset.rep_list:
        return {"rep_list": dataset.rep_list}

    # TODO: Is this required now that we do that immediately?
    # rep_list = [rep for rep in pickle.loads(dataset.dataframe[mol_col][0]).GetPropNames() if rep.startswith('atom')]
    # return {"rep_list": rep_list}


@ cime_api.route('/upload_sdf', methods=['POST'])
def upload_sdf():
    _log.info(f'Received new file to upload')
    file = request.files.get("myFile")
    # TODO: Check for valid file, i.e. type/size/...
    if not file or file.filename == '':
        abort('No valid file provided')

    _log.info(f'Uploading file {file.filename}')

    frame, rep_list = sdf_to_df(
        file, smiles_col_name=smiles_col, mol_col_name=mol_col)

    # TODO: Define format of saved file, i.e. include user?
    saved = get_cime_dbo().save_file({
        'name': file.filename,
        'dataframe': frame,
        'rep_list': rep_list
    })

    return {
        'filename': saved.name,
        'id': saved.id
    }


# TODO: WE REMOVED THIS ROUTE
# @cime_api.route('/get_csv/', methods=['GET'])
@ cime_api.route('/get_csv/<id>/', methods=['GET'])
@ cime_api.route('/get_csv/<id>/<modifiers>', methods=['GET'])
def sdf_to_csv(id, modifiers=None):
    start_time = time.time()
    if modifiers:
        # split and trim modifier string
        descriptor_names_no_lineup.extend(
            [x.strip() for x in modifiers.split(";")])

    dataset = get_cime_dbo().get_dataset_by(id=id)

    if not dataset:
        abort(404)

    frame = dataset.dataframe
    all_smiles = frame[smiles_col]

    # mols = frame[mol_col]
    # frame = frame.drop(columns=[mol_col])

    # sort such that the name column comes first and the smiles column comes second
    sm = frame[smiles_col]
    name = frame["ID"]
    frame.drop(labels=[smiles_col, "ID"], axis=1, inplace=True)
    frame.insert(0, "ID", name)
    frame.insert(1, smiles_col, sm)

    has_fingerprint = False
    new_cols = []
    for col in frame.columns:
        modifier = ""
        col_name = col


        if col.startswith(tuple(descriptor_names_no_lineup)):
            # this modifier tells lineup that the column should not be viewed at all (remove this modifier, if you want to be able to add the column with the sideview of lineup)
            modifier = '%s"noLineUp":true,' % modifier
            # this modifier tells lineup that the columns belong to a certain group
            modifier = '%s"featureLabel":"%s",' % (modifier, col.split("_")[0])
            split_col = col.split("_")
            col_name = col.replace(
                split_col[0]+"_", "") + " (" + split_col[0] + ")"
        elif col.startswith(tuple(descriptor_names_show_lineup)):
            # modifier = '%s"showLineUp":true,'%modifier # this modifier tells lineup that the column should be initially viewed
            # this modifier tells lineup that the columns belong to a certain group
            modifier = '%s"featureLabel":"%s",' % (modifier, col.split("_")[0])
            split_col = col.split("_")
            col_name = col.replace(
                split_col[0]+"_", "") + " (" + split_col[0] + ")"
        # else:
            # modifier = '%s"showLineUp":true,'%modifier # this modifier tells lineup that the column should be initially viewed

        elif col == smiles_col or col.startswith(smiles_prefix):
            # this modifier tells lineup that a structure image of this smiles string should be loaded
            modifier = '%s"hideLineUp":true,"imgSmiles":true,' % modifier

            split_col = col.split("_")
            if len(split_col) >= 2:
                col_name = col.replace(
                    split_col[0]+"_", "") + " (" + split_col[0] + ")"

        if col.startswith(fingerprint_modifier):
            has_fingerprint = True
        else:
            modifier = '%s"project":false,'%modifier
        
        if col == "ID":
            modifier = '%s"dtype":"string",' % modifier

        # remove the last comma
        new_cols.append("%s{%s}" % (col_name, modifier[0:-1]))

    frame.columns = new_cols

    if not has_fingerprint:  # when there are no morgan fingerprints included in the dataset, calculate them now
        # TODO: We should change that to MolFromSmiles instead of the pickled Mol
        # fps = pd.DataFrame([list(AllChem.GetMorganFingerprintAsBitVect(
        #     pickle.loads(mol), 5, nBits=256)) for mol in mols])
        fps = pd.DataFrame([list(AllChem.GetMorganFingerprintAsBitVect(
            Chem.MolFromSmiles(smi), 5, nBits=256)) for smi in all_smiles])
        fps.columns = [
            'fingerprint_%s{"noLineUp":true,"featureLabel": "fingerprint"}' % fp for fp in fps]
        frame = frame.join(fps)

    csv_buffer = StringIO()
    frame.to_csv(csv_buffer, index=False)

    delta_time = time.time()-start_time
    print("took", time.strftime('%H:%M:%S', time.gmtime(
        delta_time)), f"to load file {id}")
    print(f"took {delta_time/60} min {delta_time % 60} s to load file {id}")
    # print("get_csv time elapsed [s]:", time.time()-start_time)

    return csv_buffer.getvalue()


@ cime_api.route('/get_difference_highlight', methods=['OPTIONS', 'POST'])
def smiles_to_difference_highlight():
    if request.method == 'POST':
        smilesA = request.form.get("smilesA").split(",")
        smilesB = request.form.get("smilesB").split(",")

        if type(smilesA) == list:
            if len(smilesA) <= 1:
                smilesA = smilesA[0]
                molA = Chem.MolFromSmiles(smilesA)
            else:
                molA = get_mcs(smilesA)
                smilesA = Chem.MolToSmiles(molA)
        else:
            molA = Chem.MolFromSmiles(smilesA)

        if type(smilesB) == list:
            if len(smilesB) <= 1:
                smilesB = smilesB[0]
                molB = Chem.MolFromSmiles(smilesB)
            else:
                molB = get_mcs(smilesB)
                smilesB = Chem.MolToSmiles(molB)
        else:
            molB = Chem.MolFromSmiles(smilesB)

        # need a copy because otherwise the structure is messed up
        molA_cpy = copy.deepcopy(molA)  # molA#Chem.MolFromSmarts(smilesA)
        molB_cpy = copy.deepcopy(molB)  # molB#Chem.MolFromSmarts(smilesB)
        # mol_cpy = molB_cpy
        patt = get_mcs([molA_cpy, molB_cpy])

        mol = molB

        highlight_atoms = set(range(len(mol.GetAtoms()))) - \
            set(mol.GetSubstructMatch(patt))
        highlight_bonds = [bond.GetIdx() for bond in mol.GetBonds() if bond.GetBeginAtomIdx(
        ) in highlight_atoms or bond.GetEndAtomIdx() in highlight_atoms]
        print(highlight_atoms, mol.GetSubstructMatch(patt))

        highlight_atom_colors = {i: (0, 0.49, 0.68) for i in highlight_atoms}
        highlight_bond_colors = {i: (0, 0.49, 0.68) for i in highlight_bonds}
        # print(mol_cpy.GetSubstructMatch(patt), mol.GetSubstructMatch(patt), highlight_atoms)

        d = Chem.Draw.rdMolDraw2D.MolDraw2DCairo(200, 200)
        Chem.Draw.rdMolDraw2D.PrepareAndDrawMolecule(d, mol, highlightAtoms=highlight_atoms, highlightBonds=highlight_bonds,
                                                     highlightAtomColors=highlight_atom_colors, highlightBondColors=highlight_bond_colors)

        d.FinishDrawing()

        stream = BytesIO(d.GetDrawingText())

        img_str = base64.b64encode(stream.getvalue())
        stream.close()
        return {"data": img_str.decode("utf-8")}
        # return {"data": mol_to_base64(mol)}
    else:
        return {}


@ cime_api.route('/get_mol_img', methods=['OPTIONS', 'POST'])
def smiles_to_img_post():
    if request.method == 'POST':
        smiles = request.form.get("smiles")
        img = smiles_to_base64(smiles)
        return {"data": img}
    else:
        return {}


@ cime_api.route('/get_mol_imgs', methods=['OPTIONS', 'POST'])
def smiles_list_to_imgs():
    if request.method == 'POST':
        smiles_list = request.form.getlist("smiles_list")
        current_rep = request.form.get("current_rep")
        contourLines = request.form.get("contourLines")
        scale = request.form.get("scale")
        sigma = request.form.get("sigma")
        showMCS = request.form.get("showMCS")
        width = request.form.get("width")
        doAlignment = request.form.get("doAlignment") == "true"

        if len(smiles_list) == 0:
            return {"error": "empty SMILES list"}
        # if len(smiles_list) == 1:
        #    return {"img_lst": [smiles_to_base64(smiles_list[0])]}

        mol_lst = []
        error_smiles = []
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mol_lst.append(mol)
            else:
                error_smiles.append(smiles)

        if len(mol_lst) > 1:
            patt = get_mcs(mol_lst)

            # TemplateAlign.rdDepictor.Compute2DCoords(patt)
        else:
            patt = Chem.MolFromSmiles("*")

        df = None
        if current_rep != "Common Substructure":
            # TODO: Use id instead of filename
            id = request.form.get("filename")
            dataset = get_cime_dbo().get_dataset_by(id=id)
            df = dataset.dataframe if dataset else None

        if doAlignment:
            # set the coordinates of the core pattern based on the coordinates of the first molecule     
            match = mol_lst[0].GetSubstructMatch(patt)
            conf = Chem.Conformer(patt.GetNumAtoms())
            TemplateAlign.rdDepictor.Compute2DCoords(mol_lst[0])
            mconf = mol_lst[0].GetConformer(0)
            for i,aidx in enumerate(match):
                conf.SetAtomPosition(i,mconf.GetAtomPosition(aidx))
            patt.AddConformer(conf)

        img_lst = []
        for mol in mol_lst:
            if doAlignment:  # if user disables alignment, skip
                # if no common substructure was found, skip the alignment
                if(patt and Chem.MolToSmiles(patt) != "*"):
                    TemplateAlign.rdDepictor.Compute2DCoords(mol)
                    match = mol.GetSubstructMatch(patt)
                    TemplateAlign.AlignMolToTemplate2D(mol,patt,match=match,clearConfs=True)
            if current_rep == "Common Substructure":
                img_lst.append(mol_to_base64_highlight_substructure(
                    mol, patt, width=width))
            else:
                img_lst.append(mol_to_base64_highlight_importances(
                    df, mol, patt, current_rep, contourLines, scale, sigma, showMCS, smiles_col, mol_col, width))

        return {"img_lst": img_lst, "error_smiles": error_smiles}
    else:
        return {}


@ cime_api.route('/get_common_mol_img', methods=['OPTIONS', 'POST'])
def smiles_list_to_common_substructure_img():
    if request.method == 'POST':
        smiles_list = request.form.getlist("smiles_list")
        if len(smiles_list) == 0:
            return {"error": "empty SMILES list"}
        if len(smiles_list) == 1:
            ret = smiles_to_base64(smiles_list[0])
            return {"data": ret, "smiles": smiles_list[0]}

        mol_lst = []
        error_smiles = []
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mol_lst.append(mol)
            else:
                error_smiles.append(smiles)

        m = get_mcs(mol_lst)
        pil_img = Draw.MolToImage(m)

        buffered = BytesIO()
        pil_img.save(buffered, format="JPEG")
        img_str = base64.b64encode(buffered.getvalue())
        return {"data": img_str.decode("utf-8"), "smiles": Chem.MolToSmiles(m)}
    else:
        return {}


# --------- search & filter ---------
@ cime_api.route('/get_substructure_count', methods=['OPTIONS', 'POST'])
def smiles_list_to_substructure_count():
    if request.method == 'POST':
        smiles_list = request.form.get("smiles_list").split(",")
        filter_smiles = request.form.get("filter_smiles")

        if len(smiles_list) == 0:
            return {"error": "empty SMILES list"}

        patt = Chem.MolFromSmiles(filter_smiles)
        if patt:
            substructure_counts = [(smiles, len(Chem.MolFromSmiles(smiles).GetSubstructMatch(
                patt))) for smiles in smiles_list if Chem.MolFromSmiles(smiles) is not None]
            return {"substructure_counts": substructure_counts}
        return {"error": "invalid SMILES filter"}
    else:
        return {}
