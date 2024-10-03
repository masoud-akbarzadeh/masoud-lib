from setuptools import setup, find_packages

setup(
    name='masoud_lib',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        pandas,
        numpy,
        matplotlib,
        pyfluids,
        scipy,
        xarray,
        netCDF4,
        glob,
        os,
        datetime,
        seaborn,
    ],
    author='Masoud Akbarzadeh',
    author_email='masoud.akbarzadeh.edu@gmail.com',
    description='A collection of custom functions for data processing, visualization, aerosol analysis, and utilities.',
)