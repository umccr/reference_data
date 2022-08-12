# hg38 reference genome (no alt contigs)

The `hg38_noalt.fa` reference genome contains only the main contigs from `hg38.fa`. It is generally used to generate
non-alt only alignments that are required for some downstream tools.

To create this reference genome, the main contigs are taken from the existing hg38 assembly and then indexed with BWA:

```bash
dvc pull -r storage-s3 reference_data/genomes/hg38.fa{,.fai}
samtools faidx > reference_data/genomes/hg38_noalt.fa \
  reference_data/genomes/hg38.fa \
  chr{1..22} chr{X,Y,M}
bwa index reference_data/genomes/hg38_noalt.fa
```
