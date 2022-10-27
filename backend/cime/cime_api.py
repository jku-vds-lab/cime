import base64
import copy
import json
import logging
import multiprocessing
import pickle
import time
from io import BytesIO, StringIO
from typing import Callable, Dict, List

import pandas as pd
from flask import Blueprint, abort, current_app, jsonify, request
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, MACCSkeys, TemplateAlign
from rdkit.Chem.PropertyMol import PropertyMol

from .constants import (
    descriptor_names_no_lineup,
    descriptor_names_show_lineup,
    fingerprint_modifier,
    id_col_name,
    maccs_modifier,
    mol_col,
    morgan_modifier,
    smiles_col,
)
from .db import ACimeDBO
from .mol import get_mcs, mol_to_base64_highlight_importances, mol_to_base64_highlight_substructure, smiles_to_base64

_log = logging.getLogger(__name__)

cime_api = Blueprint("cime", __name__)


def get_cime_dbo() -> ACimeDBO:
    return current_app.config["CIME_DBO"]


def try_parse_float(value):
    try:
        return float(value)
    except ValueError:
        return value


def sdf_to_df_generator(file, smiles_col_name=smiles_col, column_metadata=None):
    # Inspired by https://github.com/rdkit/rdkit/blob/master/rdkit/Chem/PandasTools.py#L452
    # adapt LoadSDF method from rdkit. this version creates a pandas dataframe excluding atom specific properties and does not give every property the object type
    _log.info("Starting processing of SDF file")
    rep_list = set()

    morgan_fp = None
    maccs_fp = None
    flattened_columns = None

    if column_metadata is not None:
        morgan_fp = column_metadata.get("morganFingerprint")
        maccs_fp = column_metadata.get("maccsFingerprint")
        flattened_columns = column_metadata.get("columns")

    _log.info(flattened_columns)

    with Chem.ForwardSDMolSupplier(file) as suppl:
        i = 0
        for mol in suppl:
            if mol is None:
                continue

            # Define the row to be filled
            row = {}

            # TODO: Make sure that *all* structures have either, otherwise conflicts arise
            id: str = mol.GetProp("_Name") if mol.HasProp("_Name") else str(i)

            smiles = Chem.MolToSmiles(mol)
            if smiles is None:
                continue

            # Add prop names
            for k in mol.GetPropNames():
                if k.startswith("atom"):
                    rep_list.add(k)
                elif flattened_columns is None or (k in flattened_columns and flattened_columns[k]["include"]):
                    row[k] = try_parse_float(mol.GetProp(k))

            # Add Mol from smiles
            if smiles_col_name is not None:
                try:
                    row[smiles_col_name] = smiles
                except:
                    row[smiles_col_name] = None

            if morgan_fp is not None and morgan_fp["include"]:
                radius = morgan_fp["radius"]
                fps = {
                    f"{morgan_modifier}{radius}_{i}": key
                    for i, key in enumerate(AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=morgan_fp["bits"]))
                }
                row.update(fps)

            if maccs_fp is not None and maccs_fp["include"]:
                fps = {f"{maccs_modifier}_{i}": key for i, key in enumerate(MACCSkeys.GenMACCSKeys(mol))}
                row.update(fps)

            yield id, smiles, row, {"mol": PropertyMol(mol), "column": smiles_col_name}, list(rep_list)

            i += 1
            if i % 1000 == 0:
                _log.info(f"Processed {i} entries")

    _log.info(f"Finished processing of {i} entries")


def sdf_to_df(file, id_col_name=id_col_name, smiles_col_name=smiles_col, mol_col_name=mol_col, first_only=False):
    # Inspired by https://github.com/rdkit/rdkit/blob/master/rdkit/Chem/PandasTools.py#L452
    # adapt LoadSDF method from rdkit. this version creates a pandas dataframe excluding atom specific properties and does not give every property the object type
    _log.info("Starting processing of SDF file")
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
                    row[k] = try_parse_float(mol.GetProp(k))

            # Add name
            if mol.HasProp("_Name"):
                row[id_col_name] = mol.GetProp("_Name")

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

            records.append(row)

            i += 1
            if i % 1000 == 0:
                _log.info(f"Processed {i} entries")

            if first_only:
                break

    has_fingerprint = False
    frame = pd.DataFrame(records)
    for col in frame.columns:
        if col.startswith(fingerprint_modifier):
            has_fingerprint = True
            break

    if not has_fingerprint:
        all_smiles = frame[smiles_col]
        fps = pd.DataFrame([list(AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi), 5, nBits=256)) for smi in all_smiles])
        fps = fps.rename(mapper=lambda a: f"fingerprint_{a}", axis="columns")
        frame = frame.join(fps)

    _log.info(f"Finished processing of {i} entries")
    return frame, list(rep_list)


@cime_api.route("/get_uploaded_files_list", methods=["GET"])
def get_uploaded_files_list():
    # TODO: Check permissions will probably be done in dbo implementation
    return jsonify([{"id": d.id, "name": d.name} for d in get_cime_dbo().get_datasets()])


@cime_api.route("/delete_file/<id>", methods=["GET"])
def delete_uploaded_file(id):
    # TODO: Is this required?
    if id == "test.sdf" or id == "test":
        pass
        # return {"deleted": "false", "error": "can't delete default file"}

    # TODO: Rewrite to use id
    deleted = get_cime_dbo().delete_dataset_by(id=id)

    # TODO: Refactor to true and false instead of strings
    return {"deleted": "true" if deleted else "false"}


@cime_api.route("/get_atom_rep_list/<id>", methods=["GET"])
def get_atom_rep_list(id):
    dataset = get_cime_dbo().get_dataset_by(id=id)

    if not dataset or not dataset.rep_list:
        abort(404)

    if dataset.rep_list:
        return {"rep_list": dataset.rep_list}


@cime_api.route("/upload_sdf", methods=["POST"])
def upload_sdf():
    _log.info(f"Received new file to upload")
    file = request.files.get("myFile")
    # TODO: Check for valid file, i.e. type/size/...
    if not file or file.filename == "":
        abort("No valid file provided")

    _log.info(f"Uploading file {file.filename}")

    frame, rep_list = sdf_to_df(file, smiles_col_name=smiles_col, mol_col_name=mol_col)

    # TODO: Define format of saved file, i.e. include user?
    saved = get_cime_dbo().save_file({"name": file.filename, "dataframe": frame, "rep_list": rep_list})

    return {"filename": saved.name, "id": saved.id}


@cime_api.route("/get_csv/<id>/", methods=["GET"])
@cime_api.route("/get_csv/<id>/<modifiers>", methods=["GET"])
def sdf_to_csv(id, modifiers=None):
    start_time = time.time()
    if modifiers:
        # split and trim modifier string
        # TODO: This modifies the original list!!
        descriptor_names_no_lineup.extend([x.strip() for x in modifiers.split(";")])

    dataset = get_cime_dbo().get_dataset_by(id=id)

    domains = dataset.domains
    flattened_columns = dataset.column_metadata["columns"]
    smiles_column_metadata = next((c for c in flattened_columns.values() if c.get("smiles")), None)
    smiles_col = smiles_column_metadata["id"] if smiles_column_metadata else None
    views = dataset.column_metadata.get("views")

    if not dataset:
        abort(404)

    if hasattr(dataset, "get_chunked_dataframe"):
        _log.info("Dataset supports chunked streaming")

        # @stream_with_context
        def generate():
            for i, frame in enumerate(dataset.get_chunked_dataframe(250)):
                # sort such that the name column comes first and the smiles column comes second
                name = frame[id_col_name]
                frame.drop(labels=[id_col_name], axis=1, inplace=True)
                frame.insert(0, id_col_name, name)

                if smiles_col:
                    sm = frame[smiles_col]
                    frame.drop(labels=[smiles_col], axis=1, inplace=True)
                    frame.insert(1, smiles_col, sm)

                new_cols = []
                for col in frame.columns:
                    metadata = flattened_columns.get(col) or {}
                    modifier = {
                        # TODO: Unify this with the actual metadata from the dataset?
                        "noLineUp": not metadata.get("lineup"),
                        "smiles": metadata.get("smiles"),
                        "featureLabel": metadata.get("group") or "",
                        "dtype": "string" if col == id_col_name else None,
                        "project": col.startswith(maccs_modifier)
                        or col.startswith(morgan_modifier)
                        or col.startswith(fingerprint_modifier),
                    }

                    if views and len(views) > 0:
                        modifier["view"] = []
                        for c, view in enumerate(views):
                            if "x" in view and view["x"] == col:
                                modifier["view"].append(f"x{c + 1}")
                            if "y" in view and view["y"] == col:
                                modifier["view"].append(f"y{c + 1}")

                    if col in domains:
                        modifier["domain"] = domains[col]

                    # remove the last comma
                    new_cols.append(f"{col}{json.dumps(modifier)}")

                frame.columns = new_cols
                csv_buffer = StringIO()
                frame.to_csv(csv_buffer, index=False, header=(i == 0))
                yield csv_buffer.getvalue()

        return "".join(list(generate()))
        # return Response(generate(), mimetype='text/csv', headers={'X-Accel-Buffering': 'no'})
    else:
        frame = dataset.dataframe

        # sort such that the name column comes first and the smiles column comes second
        name = frame[id_col_name]
        frame.drop(labels=[id_col_name], axis=1, inplace=True)
        frame.insert(0, id_col_name, name)

        if smiles_col:
            sm = frame[smiles_col]
            frame.drop(labels=[smiles_col], axis=1, inplace=True)
            frame.insert(1, smiles_col, sm)

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
                col_name = col.replace(split_col[0] + "_", "") + " (" + split_col[0] + ")"
            elif col.startswith(tuple(descriptor_names_show_lineup)):
                # modifier = '%s"showLineUp":true,'%modifier # this modifier tells lineup that the column should be initially viewed
                # this modifier tells lineup that the columns belong to a certain group
                modifier = '%s"featureLabel":"%s",' % (modifier, col.split("_")[0])
                split_col = col.split("_")
                col_name = col.replace(split_col[0] + "_", "") + " (" + split_col[0] + ")"
            # else:
            # modifier = '%s"showLineUp":true,'%modifier # this modifier tells lineup that the column should be initially viewed

            elif smiles_col and col == smiles_col:
                # this modifier tells lineup that a structure image of this smiles string should be loaded
                modifier = '%s"hideLineUp":true,"imgSmiles":true,' % modifier

                split_col = col.split("_")
                if len(split_col) >= 2:
                    col_name = col.replace(split_col[0] + "_", "") + " (" + split_col[0] + ")"

            if not col.startswith(maccs_modifier) and not col.startswith(morgan_modifier) and not col.startswith(fingerprint_modifier):
                modifier = '%s"project":false,' % modifier

            if col == id_col_name:
                modifier = '%s"dtype":"string",' % modifier

            # remove the last comma
            new_cols.append("%s{%s}" % (col_name, modifier[0:-1]))

        frame.columns = new_cols

        csv_buffer = StringIO()
        frame.to_csv(csv_buffer, index=False)

        delta_time = time.time() - start_time
        print("took", time.strftime("%H:%M:%S", time.gmtime(delta_time)), f"to load file {id}")
        print(f"took {delta_time/60} min {delta_time % 60} s to load file {id}")
        # print("get_csv time elapsed [s]:", time.time()-start_time)

        return csv_buffer.getvalue()


@cime_api.route("/get_difference_highlight", methods=["OPTIONS", "POST"])
def smiles_to_difference_highlight():
    if request.method == "POST":
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

        highlight_atoms = set(range(len(mol.GetAtoms()))) - set(mol.GetSubstructMatch(patt))
        highlight_bonds = [
            bond.GetIdx() for bond in mol.GetBonds() if bond.GetBeginAtomIdx() in highlight_atoms or bond.GetEndAtomIdx() in highlight_atoms
        ]
        print(highlight_atoms, mol.GetSubstructMatch(patt))

        highlight_atom_colors = {i: (0, 0.49, 0.68) for i in highlight_atoms}
        highlight_bond_colors = {i: (0, 0.49, 0.68) for i in highlight_bonds}
        # print(mol_cpy.GetSubstructMatch(patt), mol.GetSubstructMatch(patt), highlight_atoms)

        d = Chem.Draw.rdMolDraw2D.MolDraw2DCairo(200, 200)
        Chem.Draw.rdMolDraw2D.PrepareAndDrawMolecule(
            d,
            mol,
            highlightAtoms=highlight_atoms,
            highlightBonds=highlight_bonds,
            highlightAtomColors=highlight_atom_colors,
            highlightBondColors=highlight_bond_colors,
        )

        d.FinishDrawing()

        stream = BytesIO(d.GetDrawingText())

        img_str = base64.b64encode(stream.getvalue())
        stream.close()
        return {"data": img_str.decode("utf-8")}
        # return {"data": mol_to_base64(mol)}
    else:
        return {}


@cime_api.route("/get_mol_img", methods=["OPTIONS", "POST"])
def smiles_to_img_post():
    if request.method == "POST":
        smiles = request.form.get("smiles")
        img = smiles_to_base64(smiles)
        return {"data": img}
    else:
        return {}


# Utility function parallizing array operations
def parallelized(func: Callable, l: List):
    def wrapper(arg_with_index):
        i, arg = arg_with_index
        return (i, func(arg_with_index))

    num_cores = multiprocessing.cpu_count()
    unsorted_result = filter(None, Parallel(n_jobs=num_cores)(delayed(wrapper)(m) for m in list(enumerate(l))))
    return [x[1] for x in sorted(unsorted_result, key=lambda x: x[0])]


@cime_api.route("/get_mol_imgs", methods=["OPTIONS", "POST"])
def smiles_list_to_imgs():
    if request.method == "POST":
        ids = request.form.getlist("ids")
        current_rep = request.form.get("current_rep")
        contourLines = request.form.get("contourLines")
        scale = request.form.get("scale")
        sigma = request.form.get("sigma")
        showMCS = request.form.get("showMCS")
        width = request.form.get("width")
        doAlignment = request.form.get("doAlignment") == "true"

        if len(ids) == 0:
            return {"error": "empty ids list"}
        # if len(smiles_list) == 1:
        #    return {"img_lst": [smiles_to_base64(smiles_list[0])]}

        # TODO: Use id instead of filename
        id = request.form.get("filename")
        dataset = get_cime_dbo().get_dataset_by(id=id)
        if not dataset:
            return {"error": "dataset not found"}

        id_to_mol: Dict[str, str] = dataset.get_mols_by_ids(ids)
        mol_lst = list(id_to_mol.values())

        # TODO: This entire endpoint could be returning a single image if we would fetch the MCS and alignment beforehand?
        # This would allow us to "lazy" load the images and only load the ones that are actually needed
        if len(mol_lst) > 1:
            patt = get_mcs(mol_lst)
        else:
            patt = Chem.MolFromSmiles("*")

        if doAlignment:
            first = mol_lst[0]
            # set the coordinates of the core pattern based on the coordinates of the first molecule
            match = first.GetSubstructMatch(patt)
            conf = Chem.Conformer(patt.GetNumAtoms())
            TemplateAlign.rdDepictor.Compute2DCoords(first)
            mconf = first.GetConformer(0)
            for i, aidx in enumerate(match):
                conf.SetAtomPosition(i, mconf.GetAtomPosition(aidx))
            patt.AddConformer(conf)

        img_lst = []

        def generate_image(arg):
            i, mol = arg
            if doAlignment:  # if user disables alignment, skip
                # if no common substructure was found, skip the alignment
                if patt and Chem.MolToSmiles(patt) != "*":
                    TemplateAlign.rdDepictor.Compute2DCoords(mol)
                    match = mol.GetSubstructMatch(patt)
                    TemplateAlign.AlignMolToTemplate2D(mol, patt, match=match, clearConfs=True)
            if current_rep == "Common Substructure":
                return mol_to_base64_highlight_substructure(mol, patt, width=width)
            else:
                return mol_to_base64_highlight_importances(mol, patt, current_rep, contourLines, scale, sigma, showMCS, width)

        # parallelize the image creation for significant speedups
        img_lst = parallelized(generate_image, id_to_mol.values())

        return {"img_lst": img_lst, "error_smiles": []}
    else:
        return {}


@cime_api.route("/get_common_mol_img", methods=["OPTIONS", "POST"])
def smiles_list_to_common_substructure_img():
    if request.method == "POST":
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
@cime_api.route("/get_substructure_count", methods=["OPTIONS", "POST"])
def smiles_list_to_substructure_count():
    if request.method == "POST":
        smiles_list = request.form.get("smiles_list").split(",")
        filter_smiles = request.form.get("filter_smiles")

        if len(smiles_list) == 0:
            return {"error": "empty SMILES list"}

        patt = Chem.MolFromSmiles(filter_smiles)
        if patt:
            substructure_counts = [
                (smiles, len(Chem.MolFromSmiles(smiles).GetSubstructMatch(patt)))
                for smiles in smiles_list
                if Chem.MolFromSmiles(smiles) is not None
            ]
            return {"substructure_counts": substructure_counts}
        return {"error": "invalid SMILES filter"}
    else:
        return {}
