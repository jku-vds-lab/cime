from setuptools import setup

setup(
    name='cime',
    packages=['cime'],
    include_package_data=True,
    install_requires=[
        'flask',
        'flask_sqlalchemy',
        'projection-space-explorer @ git+https://github.com/jku-vds-lab/projection-space-explorer.git@mschoefl/testmigration#egg=projection-space-explorer&subdirectory=backend',
        'rdkit-pypi',
        'joblib',
        'flask-cors'
    ],
)
