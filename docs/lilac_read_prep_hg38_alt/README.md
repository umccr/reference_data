# LILAC hg38 HLA BED

> This README.md describes how the `reference_data/hmftools/lilac/hla.38.alt.umccr.bed` file was generated.

LILAC only considers reads aligned to HLA-A, HLA-B, and HLA-C on the main chr6 contig and will ignore reads mapped to
these genes present on ALT contigs. Hence to use LILAC with a standard hg38 assembly reads must be extracted and then
realigned to a genome without ALT contigs. Here I define a set of regions across the hg38 assembly from which reads
potentially derived from HLA-ABC are to be extracted.

Two approaches to obtain the final list of regions are taken:

1. processing and selection of GENCODE annotations
1. identification of regions with homology to HLA-ABC

## GENCODE-defined regions

> Specifically decided to use exons + slop rather than transcripts as some HLA transcripts on ALT contigs at 10kb+ in
> length. Moreover, LILAC specifically uses exons to select reads for processing. The complete HLA genes on the main
> chr6 contig are fully covered by this approach, so these regions should capture all relevant reads mapping elsewhere.
> However, we also explicitly add the HLA-ABC genes from the main chr6 contig with a 1000 bp slop to align with LILAC
> read slicing procedure.

Download required files

```bash
mkdir -p data/

urls='
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.chr_patch_hapl_scaff.annotation.gff3.gz
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromAlias.txt
http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
'
parallel wget -P data/ ::: ${urls}
gzip -d data/gencode.v41.chr_patch_hapl_scaff.annotation.gff3.gz
```

Select all HLA genes, transcripts, and CDS across all hg38 contigs then merge regions

```bash
mkdir -p 1_regions_gencode/

for symbol in hla_{a,b,c}; do
  hgnc_symbol=$(tr '_' '-' <<< ${symbol^^});
  resp=$(
    curl \
        -s \
        -H 'Accept:application/json' \
        "https://rest.genenames.org/search/symbol/${hgnc_symbol}"
  );
  hgnc_id=$(jq -r <<< ${resp} '.response.docs[] | .hgnc_id');
  echo ${symbol} ${hgnc_id};
done | tee 1_regions_gencode/hgnc_ids.txt

./scripts/get_regions.py > 1_regions_gencode/regions_gencode_unmerged.bed
./scripts/convert_gencode_ucsc_contigs.py | sort | uniq > 1_regions_gencode/regions_ucsc_unmerged.bed

bedtools sort -i 1_regions_gencode/regions_ucsc_unmerged.bed | \
  bedtools slop -b 500 -g data/hg38.chrom.sizes -i - | \
  bedtools merge -c 4 -o distinct -i - > 1_regions_gencode/regions_ucsc.bed

# Add in full main HLA-ABC genes with 1000 bp slop to mirror HMF read selection
grep '^chr6.\+\tgene\t.\+gene_name=HLA-[ABC]' data/gencode.v41.chr_patch_hapl_scaff.annotation.gff3 | \
  awk '{ print $1, $4-1, $5 }' | \
  tr ' ' '\t' | \
  bedtools sort -i - | \
  bedtools slop -b 1000 -g data/hg38.chrom.sizes -i - > 1_regions_gencode/regions_main_genes.bed

# Combine
cat \
  1_regions_gencode/regions_main_genes.bed \
  <(cut -f1-3 -d$'\t' 1_regions_gencode/regions_ucsc.bed) | \
  bedtools sort -i - | \
  bedtools merge -i - > 1_regions_gencode/final.bed
```

## Homologous regions

> Done under the assumption that in the context of LILAC analysis we don't really care about reads mapped to
> non-homologous regions of HLA transcripts present elsewhere in hg38 since these reads should never align back to the
> HLA-ABC on the main chr6 contig by definition.

Search subject

```bash
mkdir -p 2_regions_homologous/

contigs=$(awk '$1 ~ /chr6|chrUn/ { print $1 }' ~/repos/reference_data/reference_data/genomes/hg38.fa.fai)

samtools faidx > 2_regions_homologous/genome_contigs.fasta \
  ~/repos/reference_data/reference_data/genomes/hg38.fa ${contigs}
```

Search query, rename FASTA ids

```bash
awk '{ print $1 ":" $2 "-" $3 }' ~/repos/reference_data/reference_data/hmftools/lilac/hla.38.bed |
    samtools faidx > 2_regions_homologous/hmf_hla_genes.fasta \
    ~/repos/reference_data/reference_data/genomes/hg38.fa \
    --region-file -

rename_map=$(cat <<EOF | tr ' ' '\t'
chr6:29940260-29946884 HLA-A
chr6:31352872-31358188 HLA-B
chr6:31267749-31273130 HLA-C
EOF
)
tmp_fp=$(mktemp $(pwd -P)/tmp.XXXXXX)
seqkit replace > ${tmp_fp} \
  -p '(.+)$' \
  -r '{kv}' \
  -k <(echo "${rename_map}") \
  2_regions_homologous/hmf_hla_genes.fasta
mv ${tmp_fp} 2_regions_homologous/hmf_hla_genes.fasta
```

Run search

```bash
blat \
  -out=blast9 2_regions_homologous/genome_contigs.fasta \
  2_regions_homologous/hmf_hla_genes.fasta \
  2_regions_homologous/blat_results.tsv
```

Process results

```bash
grep -v '#' 2_regions_homologous/blat_results.tsv | \
  awk '
    BEGIN {
      OFS="\t";
    }
    $4 >= 50 {
      d[1] = $9;
      d[2] = $10;
      asort(d);
      print $2, d[1]-1, d[2], "gene:" $1 ";type=blat_homologous"
    }
  ' | \
  bedtools sort -i - | \
  bedtools merge -c 4 -o distinct -i - > 2_regions_homologous/final.bed
```

## Merge regions

Require regions to be at least 50 bp. Add in all HLA contigs.

```bash
hla_contigs=$(
  awk < ~/repos/reference_data/reference_data/genomes/hg38.fa.fai \
    '
      BEGIN { OFS="\t" }
      /^HLA/ { print $1, "0", $2 }
    '
)

cat \
 1_regions_gencode/final.bed \
 <(cut -f1-3 -d$'\t' 2_regions_homologous/final.bed) \
 <(echo "${hla_contigs}") | \
 bedtools sort -i - | \
 bedtools merge -i - | \
 awk '{ s=($3+1)- $2; if (s >= 50) { print } }' > hla.38.alt.umccr.bed
```
