# UMCCR reference data

Reference data version tracking for internal pipelines like https://github.com/umccr/umccrise and https://github.com/umccr/RNAsum


## Versioning

The reference data is versioned via a manifest file, containing md5 sums for each refernece file, create via:

```
hashdeep -r hg38 > hg38-manifest.txt
```


## Python API

This package also provides a python API to access the reference data.

[reference_data/paths.yml](reference_data/paths.yml) contains default data paths and settings for common UMCCR clusters.

`from reference_data import api as refdata` is a python API that can detect the machine based on `hostname`.

Usage:

```
>>> from reference_data import api as refdata
>>> refdata.name
'spartan'

>>> refdata.set_genomes_dir('/g/data3/gx8/extras/umccrise/genomes')

>>> refdata.get_ref_file(genome='hg38', key='fa')
'/g/data3/gx8/extras/umccrise/genomes/hg38.fa'

>>> refdata.get_ref_file('hg38', 'gnomad')
'/g/data3/gx8/extras/umccrise/genomes/gnomad_genome.r2.1.common_pass_clean.norm.vcf.gz'

>>> refdata.get_ref_file('hg38', ['truth_sets', 'giab', 'bed'])
'/g/data3/gx8/local/development/bcbio/genomes/Hsapiens/hg38/validation/giab-NA12878/truth_regions.bed'

>>> refdata.get_genomes_dict('hg38')['truth_sets']['giab']['bed']
'/g/data3/gx8/local/development/bcbio/genomes/Hsapiens/hg38/validation/giab-NA12878/truth_regions.bed'
```

Available genomes:

- `"hg38"`
- `"GRCh37"`

Available keys: see [reference_data/paths.yml](reference_data/paths.yml)

Installation:

```
git clone https://github.com/umccr/reference_data
pip install reference_data
```