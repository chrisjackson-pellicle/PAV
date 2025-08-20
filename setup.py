#!/usr/bin/env python

import setuptools

# Read version number into a dictionary
version = {}
with open('plastid_annotation_validator/version.py') as fp:
    exec(fp.read(), version)


pav_scripts = ['pav']
pav_description = 'Plastid Annotation Validator'
pav_url = 'https://github.com/chrisjackson-pellicle/plDNA_annotation'
pav_entry_points = {'console_scripts': ['pav = plastid_annotation_validator.pav:main']}

setuptools.setup(name='pav',
                 version=version['__version__'],
                 include_package_data=True,
                 packages=setuptools.find_packages(),
                 author='Chris Jackson',
                 author_email='chris.jackson@rbg.vic.gov.au',
                 description=pav_description,
                 keywords='plastid annotation validator',
                 url=pav_url,
                 entry_points=pav_entry_points)
