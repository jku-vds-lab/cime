from typing import Dict, List

from .cime_api import cime_api
from .constants import id_col_name, mol_col
from .db import ACimeDBO


def create_app():
    import logging
    import os

    from flask import Flask
    from flask_cors import CORS
    from flask_sqlalchemy import SQLAlchemy

    app = Flask(__name__)
    CORS(app)

    app.config["SQLALCHEMY_DATABASE_URI"] = "sqlite://"
    app.config["CIME_FILES_DIRECTORY"] = os.path.join(os.path.dirname(__file__), "temp-files")
    app.logger.setLevel(logging.INFO)

    db = SQLAlchemy(app)

    from sqlalchemy import Column, Integer, PickleType, String

    class Dataset(db.Model):
        __tablename__ = "dataset"

        id = Column(Integer, primary_key=True)
        name = Column(String, nullable=False)
        dataframe = Column(PickleType, nullable=False)
        rep_list = Column(PickleType, nullable=False)

        def get_mols_by_ids(self, ids: List[str]) -> Dict[str, object]:
            subset = self.dataframe[self.dataframe[id_col_name].isin(ids)]
            return dict(zip(subset[id_col_name], subset[mol_col]))

    # db.drop_all()
    db.create_all()

    def setup_cime_dbo(db):
        class CimeDBO(ACimeDBO):
            def get_datasets_by(self, **kwargs):
                return Dataset.query.filter_by(**kwargs).all()

            def get_dataset_by(self, **kwargs):
                return Dataset.query.filter_by(**kwargs).first()

            def save_dataset(self, dataset):
                db.session.add(dataset)
                db.session.commit()

            def save_file(self, file):
                dataset = Dataset()
                dataset.name = file.get("name")
                dataset.dataframe = file.get("dataframe")
                dataset.rep_list = file.get("rep_list")
                db.session.add(dataset)
                db.session.commit()
                return dataset

            def delete_dataset_by(self, **kwargs):
                Dataset.query.filter_by(**kwargs).delete()
                db.session.commit()

        return CimeDBO()

    app.config["CIME_DBO"] = setup_cime_dbo(db)

    # dataset1 = Dataset()
    # dataset1.name = 'Test 1'
    # app.config['CIME_DBO'].save_dataset(dataset1)
    # dataset2 = Dataset()
    # dataset2.name = 'Test 2'
    # app.config['CIME_DBO'].save_dataset(dataset2)

    # app.register_blueprint(pse_api)
    app.register_blueprint(cime_api)
    return app
