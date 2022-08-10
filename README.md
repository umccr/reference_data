# UMCCR reference data

Reference data and API for UMCCR workflows and tools.

> Note that Umccrise reference data will be further reorganised in the future for consistency. This will be performed in
> concert with updates to Umccrise component software and modularisation.

## Contents

* [Quick start](#quick-start)
* [Installation](#installation)
* [Requirements](#requirements)
* [Usage](#usage)
* [Limitations](#limitations)
* [Data](#data)
* [License](#license)

## Quick start

```bash
# Pull the Umccrise reference data bundle
umccr_refdata pull --bundle umccrise --output_dir ./refdata_umccrise/

# Get location of the genome FASTA file
umccr_refdata locate --identifier genome_fasta
```

A Python API is also available, see [Usage](#usage).

## Installation

```bash
# Conda package
conda install -c umccr umccr_refdata

# Directly from the GitHub repo
pip install git+https://github.com/umccr/reference_data@dev
```

## Requirements

Python, and the following packages:

* PyYAML
* DVC
* DVC[s3] (*required only for S3 remotes*)
* DVC[gdrive] (*required only for gdrive remotes*)

Other software requirements:

* git

## Usage

The `umccr_refdata` package provides two major functionalities:

1. download reference data bundles for specific tools or workflows
1. retrieving path for reference files from a symbolic identifier

Both can be invoked either via the CLI or API.

### Downloading bundles

A reference data bundle represents the required set of reference data files for a particular tool or workflow e.g.
Umccrise or GRIDSS. Abstracting reference data files in this way simplifies retrieval from a DVC remote; all required
reference data files can be pulled with a single command and no unnecessary files are downloaded. See
[`refdata_information.yaml`](umccr_refdata/refdata_information.yaml) for defined bundles.

Behaviour of the bundle download process can be configured with several options:

* `--git_tag`/`git_tag`: Git tag used to checkout the DVC repository [*also accepts branch names*]
* `--git_repo_url`/`git_repo_url`: GitHub URL containing the DVC repository
* `--dvc_remote_name`/`dvc_remote_name`: DVC remote name e.g. `storage-s3`, `storage-gdrive`
* `--cache_dir`/`cache_dir`: Local filepath for the DVC cache directory

#### CLI

```bash
# Basic usage
umccr_refdata pull --bundle GRIDSS --output_dir ./reference_data/

# Set specific git tag and DVC remote
umccr_refdata pull --bundle GRIDSS --output_dir ./reference_data/ --git_tag 2.0.0 --dvc_remote_name storage-gdrive
```

#### API

```python
import pathlib


import umccr_refdata.util
import umccr_refdata.api


output_dir = pathlib.Path('reference_data')
refdata_info_fp = umccr_refdata.util.get_refdata_information_fp()
refdata_info = umccr_refdata.api.read_refdata_information(refdata_info_fp)

# Basic usage
umccr_refdata.api.pull_bundle('GRIDSS', output_dir, refdata_info)

# Set specific git tag and DVC remote
umccr_refdata.api.pull_bundle('GRIDSS', output_dir, refdata_info, git_tag='2.0.0', dvc_remote_name='storage-gdrive')
```

### Locating reference files

Individual reference data files are associated with a symbolic identifier e.g. `genome_fasta` maps to the filepath of
the hg38 FASTA file. This mapping enables external programs to refer to files with a stable identifier and hence
minimises the impact of reorganisation of reference files.

Multiple entries under a single identifier (e.g. `genome_fasta` could mapping to different genome versions) can be
differentiated by defining key-value annotations for entries (e.g. `version: hg38`) and then selecting these with the
appropriate CLI/API options.

#### CLI

```bash
# Basic usage
umccr_refdata locate --identifier genome_fasta

# Selecting different entries [demonstration example only]
umccr_refdata locate --identifier genome_fasta --match_dict '{ "version": "my_genome_version" }'
```

#### API

```python
import umccr_refdata.util
import umccr_refdata.api

refdata_info_fp = umccr_refdata.util.get_refdata_information_fp()
refdata_info = umccr_refdata.api.read_refdata_information(refdata_info_fp)

# Basic usage
umccr_refdata.api.locate_file_paths('genome_fasta', refdata_info)

# Selecting different entries [demonstration example only]
match_dict = {'version': 'my_genome_version'}
umccr_refdata.api.locate_file_paths('genome_fasta', refdata_info, match_dict=match_dict)
```

## Limitations

For gdrive remotes, you must set `GDRIVE_CREDENTIALS_DATA` prior to execution.

## Data

### Genomes

Top level directory: `./reference_data/genomes/`

| Data          | Description   |
| --            | --            |
| hg38          | 1000 Genomes Project hg38 reference genome ([link](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/)). |

### Umccrise

Top level directory: `./reference_data/umccrise/`

| Data              | Description   |
| --                | --            |
| hg38              | 1000 Genomes Project hg38 reference genome ([link](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/)). |
| hg38_noalt        | Reference hg38 genome with chr1-22, chrX, chrY, and chrM. No documentation. |
| hg38_noversion    | Appears to be the reference hg38 genome transcripts. No documentation. |
| bwa               | BWA indices for hg38. |
| gnomad            | Processed gnomAD variants ([link](https://github.com/umccr/umccrise/#gnomad)). |
| cacao             | cacao reference data ([link](https://github.com/sigven/cacao)). |
| hmf               | Selected HMFtools data (see [source](https://nextcloud.hartwigmedicalfoundation.nl/s/LTiKTd8XxBqwaiC?path=%2FHMFTools-Resources) and [processing](https://github.com/umccr/umccrise/#building-reference-data)). |
| hotspots          | Compiled mutational hotspots ([link](https://github.com/umccr/umccrise/#hotspots)). |
| panel_of_normals  | UMCCR panel of normals ([link](https://github.com/umccr/vcf_stuff/tree/master/vcf_stuff/panel_of_normals)). |
| pcgr              | PCGR reference data ([link](https://sigven.github.io/pcgr/articles/installation.html#step-1-download-data-bundle)). |
| problem_regions   | Merged problem regions ([link](https://github.com/umccr/umccrise/#problem-regions)). |
| pyensembl         | Ensembl data for PyEnsembl ([link](https://github.com/openvax/pyensembl#installation)). |
| snpeff            | snpEff database ([link](https://github.com/umccr/umccrise/#snpeff)). |
| viral             | Genomic data commons viral sequences obtained from CloudBioLinux ([link](https://s3.amazonaws.com/biodata/collections/GRCh37/viral/gdc-viral.fa)). No documentation. |

### HMFtools

Top level directory: `./reference_data/hmftools/`

> Downloaded from the Hartwig Medical Foundation Nextcloud instance on 15/02/2022
> ([link](https://nextcloud.hartwigmedicalfoundation.nl/s/LTiKTd8XxBqwaiC?path=%2FHMFTools-Resources)).

| Data                  | Description   |
| --                    | --            |
| amber                 | Selected germline heterozygous SNP sites for AMBER. |
| cobalt                | COBALT GC profile. |
| ensemble_data_cache   | Derivative data generated from the Ensembl database. |
| gene_panel            | Driver gene panel. |
| gridss                | GRIDSS PON and problem region list ([link](https://github.com/PapenfussLab/gridss/blob/bd7da/example/ENCFF356LFX.bed)). |
| known_fusions         | Curated known fusion data. |
| lilac                 | LILAC resource files. |
| linx                  | Curated known fragile sites and LINE source regions. |
| mappability           | Mappability values arcoss the reference genome. |
| repeatmasker          | RepeatMasker database. |
| sage                  | Known SAGE hotspots. |

## License

Software and code in this repository are provided under the [GNU General Public License
v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html) unless otherwise indicated.
