# Compile HMFTools reference data

Download HMFTools DNA-Resources hg38 v5.30 tarball, decompress, and extract tarball

```bash
mkdir -p data/
NC_URL_BASE='https://nextcloud.hartwigmedicalfoundation.nl/s/LTiKTd8XxBqwaiC/download'
NC_URL_PATH='%2FHMFTools-Resources%2FDNA-Resources&files=hmf_pipeline_resources.38_v5.30.gz'
wget -O data/hmf_pipeline_resources.38_v5.30.tar.gz "${NC_URL_BASE}?path=${NC_URL_PATH}"

mkdir -p data/hmf_pipeline_resources/
tar -zxvf data/hmf_pipeline_resources.38_v5.30.tar.gz -C data/hmf_pipeline_resources/
```

Set up destination directory

```bash
subdirs="
amber
cobalt
ensembl_data_cache
gene_panel
gridss
known_fusions
lilac
linx
mappability
purple
repeatmasker
sage
"
for subdir in ${subdirs}; do mkdir -p reference_data/hmftools/${subdir}/; done
```

Organise reference files in destination directory

```bash
mv data/hmf_pipeline_resources/copy_number/GermlineHetPon.38.vcf.gz reference_data/hmftools/amber/GermlineHetPon.38.vcf.gz
mv data/hmf_pipeline_resources/copy_number/GC_profile.1000bp.38.cnp reference_data/hmftools/cobalt/GC_profile.1000bp.38.cnp
mv data/hmf_pipeline_resources/common/ensembl_data/ensembl_gene_data.csv reference_data/hmftools/ensembl_data_cache/ensembl_gene_data.csv
mv data/hmf_pipeline_resources/common/ensembl_data/ensembl_protein_features.csv reference_data/hmftools/ensembl_data_cache/ensembl_protein_features.csv
mv data/hmf_pipeline_resources/common/ensembl_data/ensembl_trans_exon_data.csv reference_data/hmftools/ensembl_data_cache/ensembl_trans_exon_data.csv
mv data/hmf_pipeline_resources/common/ensembl_data/ensembl_trans_splice_data.csv reference_data/hmftools/ensembl_data_cache/ensembl_trans_splice_data.csv
mv data/hmf_pipeline_resources/common/DriverGenePanel.38.tsv reference_data/hmftools/gene_panel/DriverGenePanel.38.tsv
mv data/hmf_pipeline_resources/sv/gridss_pon_breakpoint.38.bedpe.gz reference_data/hmftools/gridss/gridss_pon_breakpoint.38.bedpe.gz
mv data/hmf_pipeline_resources/sv/gridss_pon_single_breakend.38.bed.gz reference_data/hmftools/gridss/gridss_pon_single_breakend.38.bed.gz
mv data/hmf_pipeline_resources/sv/known_fusion_data.38.csv reference_data/hmftools/known_fusions/known_fusion_data.38.csv
mv data/hmf_pipeline_resources/sv/known_fusions.38.bedpe reference_data/hmftools/known_fusions/known_fusions.38.bedpe
mv data/hmf_pipeline_resources/immune/hla_ref_aminoacid_sequences.csv reference_data/hmftools/lilac/hla_ref_aminoacid_sequences.csv
mv data/hmf_pipeline_resources/immune/hla_ref_nucleotide_sequences.csv reference_data/hmftools/lilac/hla_ref_nucleotide_sequences.csv
mv data/hmf_pipeline_resources/immune/lilac_allele_frequencies.csv reference_data/hmftools/lilac/lilac_allele_frequencies.csv
mv data/hmf_pipeline_resources/sv/fragile_sites_hmf.38.csv reference_data/hmftools/linx/fragile_sites_hmf.38.csv
mv data/hmf_pipeline_resources/sv/line_elements.38.csv reference_data/hmftools/linx/line_elements.38.csv
mv data/hmf_pipeline_resources/variants/mappability_150.38.bed.gz reference_data/hmftools/mappability/mappability_150.38.bed.gz
mv data/hmf_pipeline_resources/copy_number/cohort_germline_del_freq.38.csv reference_data/hmftools/purple/cohort_germline_del_freq.38.csv
mv data/hmf_pipeline_resources/sv/38.fa.out.gz reference_data/hmftools/repeatmasker/38.fa.out.gz
mv data/hmf_pipeline_resources/variants/ActionableCodingPanel.38.bed.gz reference_data/hmftools/sage/ActionableCodingPanel.38.bed.gz
mv data/hmf_pipeline_resources/variants/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.gz reference_data/hmftools/sage/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.gz
mv data/hmf_pipeline_resources/variants/KnownBlacklist.germline.38.bed reference_data/hmftools/sage/KnownBlacklist.germline.38.bed
mv data/hmf_pipeline_resources/variants/KnownBlacklist.germline.38.vcf.gz reference_data/hmftools/sage/KnownBlacklist.germline.38.vcf.gz
mv data/hmf_pipeline_resources/variants/KnownHotspots.germline.38.vcf.gz reference_data/hmftools/sage/KnownHotspots.germline.38.vcf.gz
mv data/hmf_pipeline_resources/variants/KnownHotspots.somatic.38.vcf.gz reference_data/hmftools/sage/KnownHotspots.somatic.38.vcf.gz
mv data/hmf_pipeline_resources/variants/SageGermlinePon.98x.38.tsv.gz reference_data/hmftools/sage/SageGermlinePon.98x.38.tsv.gz
mv data/hmf_pipeline_resources/variants/clinvar.38.vcf.gz reference_data/hmftools/sage/clinvar.38.vcf.gz
```

Download additional required reference files

```bash
wget -P reference_data/hmftools/gridss/ 'https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz'
wget -O reference_data/hmftools/lilac/hla.38.bed "${NC_URL_BASE}?path=%2FHMFTools-Resources%2FLilac%2F38&files=hla.38.bed"
```
