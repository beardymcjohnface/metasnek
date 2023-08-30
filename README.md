# MetaSnek

[![](https://img.shields.io/static/v1?label=Licence&message=MIT&color=black)](https://opensource.org/license/mit/)
[![install with PyPI](https://img.shields.io/badge/Install%20with-PyPI-brightgreen.svg?style=flat-square)](https://pypi.org/project/metasnek/)
![GitHub last commit (branch)](https://img.shields.io/github/last-commit/beardymcjohnface/MetaSnek/main)
[![Documentation Status](https://readthedocs.org/projects/metasnek/badge/?version=latest)](https://metasnek.readthedocs.io/en/latest/?badge=latest)
[![Unit Tests](https://github.com/beardymcjohnface/metasnek/actions/workflows/unit-tests.yml/badge.svg)](https://github.com/beardymcjohnface/metasnek/actions/workflows/unit-tests.yml)
[![codecov](https://codecov.io/gh/beardymcjohnface/metasnek/branch/main/graph/badge.svg?token=lCyqJhuiCN)](https://codecov.io/gh/beardymcjohnface/metasnek)

Misc functions for metagenomic pipelines

## fastq_finder

Return a python dictionary of sample names and fasta/q files from a directory or a TSV file.
If given a directory, sample names will be everything upto the first match for "_R1|_R2|_S".
Write a TSV file based on your dictionary.
[More information and examples on here](https://gist.github.com/beardymcjohnface/bb161ba04ae1042299f48a4849e917c8)
