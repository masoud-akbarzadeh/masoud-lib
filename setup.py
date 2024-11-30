from setuptools import setup, find_packages

setup(
    name='masoud_lib',
    version='0.1.0',
    packages=['masoud_lib'],
    include_package_data=True,
    package_data={'masoud_lib': ['*.py']},
    install_requires=[
        'pandas',
        'numpy',
        'matplotlib',
        'pyfluids',
        'scipy',
        'seaborn',
    ],
    author='Masoud Akbarzadeh',
    author_email='masoud.akbarzadeh.edu@gmail.com',
    description='A collection of custom functions for data processing, visualization, aerosol analysis, and utilities.',
    url='https://github.com/masoud-akbarzadeh/masoud_lib',
    license='MIT',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)
