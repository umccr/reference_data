# UMCCR reference data

Reference data and API for UMCCR workflows and tools.

> Note that Umccrise reference data will be further reorganised in the future for consistency. This will be performed in
> concert with updates to Umccrise component software and modularisation.

## Contents

* [Quick start](#quick-start)
* [Installation](#installation)
* [Requirements](#requirements)
* [Data](#data)
* [Usage](#usage)
* [License](#license)

## Quick start

> *To be completed*

```bash
umccr_refdata pull --bundle umccrise --output_dir ./refdata_umccrise/
```

## Installation

```bash
conda install -c umccr umccr_refdata
pip install git+https://github.com/umccr/reference_data
```

## Requirements

Python, and the following packages:

* PyYAML

The following packages are required to download data:

* DVC
* DVC[s3] (*required only for S3 remotes*)
* DVC[gdrive] (*required only for gdrive remotes*)

## Usage

> *To be completed*

```bash
umccr_refdata pull --bundle <predefined_bundle> --output_dir ./refdata_bundle/
```

See [`refdata_information.yaml`](umccr_refdata/refdata_information.yaml) for defined bundles.

```bash
import umccr_refdata

reference_genome_fp = umccr_refdata.get_genome(identifier='hg38')
```

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
| panel_of_normals  | UMCCR panel of normals ([link](https://github.com/umccr/vcf_stuff/tree/master/vcf_stuff/panel_of_normals). |
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
| ensemble_data_cache   | Derivated data generated from the Ensembl database. |
| gene_panel            | Driver gene panel. |
| gridss                | GRIDSS PON and problem region list ([link](https://github.com/PapenfussLab/gridss/blob/bd7da/example/ENCFF356LFX.bed)). |
| known_fusions         | Curated known fusion data. |
| linx                  | Curated known fragile sites and LINE source regions. |
| repeatmasker          | RepeatMasker database. |
| sage                  | Known SAGE hotspots. |

## License

Software and code in this repository are provided under the [GNU General Public License
v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html) unless otherwise indicated.
