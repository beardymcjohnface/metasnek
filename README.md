# MetaSnek

[![](https://img.shields.io/static/v1?label=Licence&message=MIT&color=black)](https://opensource.org/license/mit/)
[![install with PyPI](https://img.shields.io/badge/Install%20with-PyPI-brightgreen.svg?style=flat-square)](https://pypi.org/project/metasnek/)
[![Documentation Status](https://readthedocs.org/projects/metasnek/badge/?version=latest)](https://metasnek.readthedocs.io/en/latest/?badge=latest)
[![Unit Tests](https://github.com/beardymcjohnface/metasnek/actions/workflows/unit-tests.yml/badge.svg)](https://github.com/beardymcjohnface/metasnek/actions/workflows/unit-tests.yml)
[![codecov](https://codecov.io/gh/beardymcjohnface/metasnek/branch/main/graph/badge.svg?token=lCyqJhuiCN)](https://codecov.io/gh/beardymcjohnface/metasnek)

Misc functions for metagenomic pipelines

## Examples

### fastq_finder

Return a python dictionary of sample names and fasta/q files from a directory or a TSV file.
If given a directory, sample names will be everything upto the first match for "_R1|_R2".

__samples.tsv__
```text
sample1     reads/sample1_R1.fastq.gz   reads/sample1_R2.fastq.gz
sample2     reads/sample1_R2.fastq      reads/sample2_R2.fastq
sample3     reads/sample3.fasta.gz      None
```

__reads/__

```text
reads
  ├── sample1_R1.fastq.gz
  ├── sample1_R2.fastq.gz
  ├── sample2_R1.fastq
  ├── sample2_R2.fastq
  └── sample3.fasta.gz

```

```python
from metasnek import fastq_finder

samples = fastq_finder.parse_samples_to_dictionary("reads")
samples = fastq_finder.parse_samples_to_dictionary("samples.tsv")
samples
```

```text
{
'sample1': 
    {
    'R1': 'reads/sample1_R1.fastq.gz', 
    'R2': 'reads/sample1_R2.fastq.gz', 
    }, 
'sample2': 
    {
    'R1': 'reads/sample2_R1.fastq', 
    'R2': 'reads/sample2_R2.fastq'
    }
'sample3': 
    {
    'R1': 'reads/sample3.fasta.gz', 
    'R2': None
    }
}
```
