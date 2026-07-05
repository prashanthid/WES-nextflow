# WES somatic pipeline: BQSR → Mutect2 + Strelka2 → Funcotator

Paired tumor/normal WES calling. Input BAMs are assumed BWA-MEM aligned and
Picard-deduplicated, but **not** base-recalibrated — this pipeline runs BQSR
first.

## Pipeline steps

1. **BaseRecalibrator + ApplyBQSR** (GATK4) — run independently on each
   tumor and normal BAM.
2. **Mutect2** (tumor vs. matched normal), with the full best-practices
   filtering chain: `LearnReadOrientationModel`, `GetPileupSummaries`,
   `CalculateContamination`, `FilterMutectCalls`.
3. **Manta** — generates candidate small indels, fed into Strelka2 to boost
   its indel sensitivity (GATK/Illumina best practice for Strelka2).
4. **Strelka2** somatic workflow — outputs separate SNV and indel VCFs,
   merged and sorted with `bcftools`.
5. **Funcotator** — annotates both the filtered Mutect2 VCF and the merged
   Strelka2 VCF, output as MAF (default) or VCF.

## Requirements

- Nextflow ≥ 23.04
- Docker or Singularity/Apptainer (containers pinned in `nextflow.config`)
- Reference FASTA + `.fai` + `.dict`
- A Funcotator data source bundle (download via
  `gatk FuncotatorDataSourceDownloader --somatic --validate-integrity`)
- WES bait/target intervals (`.interval_list` or BED, **bgzipped + tabixed**
  BED if used by Manta/Strelka2 — see note below)
- Known-sites VCFs for BQSR (e.g. dbSNP, Mills & 1000G gold-standard indels)
- gnomAD germline resource VCF for Mutect2 contamination estimation

## Before running

Manta/Strelka2 expect the callRegions BED to be **bgzip-compressed and
tabix-indexed** (`your_targets.bed.gz` + `.tbi`), separate from the GATK
`.interval_list`/BED used elsewhere. Prepare it once:

```bash
sort -k1,1 -k2,2n targets.bed | bgzip -c > targets.bed.gz
tabix -p bed targets.bed.gz
```

The pipeline as written passes `${params.intervals}.gz` to Manta/Strelka2 —
adjust that path if your bgzipped BED lives somewhere else.

Also confirm your BAM `@RG SM:` tags. `main.nf` currently assumes
`SM:<patient_id>_tumor` / `SM:<patient_id>_normal` for the Mutect2
`-tumor`/`-normal` flags — change the `mutect2_input_ch` mapping in the
workflow block if your naming differs.

## Run

```bash
nextflow run main.nf \
  --samplesheet samplesheet_example.csv \
  --fasta /ref/Homo_sapiens_assembly38.fasta \
  --intervals /ref/wes_targets.interval_list \
  --known_sites '["/ref/dbsnp146.vcf.gz","/ref/Mills_and_1000G_gold_standard.indels.vcf.gz"]' \
  --germline_resource /ref/af-only-gnomad.hg38.vcf.gz \
  --pon /ref/1000g_pon.hg38.vcf.gz \
  --funcotator_sources /ref/funcotator_dataSources.v1.7.20200521g \
  --funcotator_ref_version hg38 \
  --outdir results \
  -resume
```

`--known_sites` is a Groovy list literal on the CLI — easier to manage via a
`params.yaml`/`params.json` file with `-params-file`:

```yaml
samplesheet: samplesheet_example.csv
fasta: /ref/Homo_sapiens_assembly38.fasta
intervals: /ref/wes_targets.interval_list
known_sites:
  - /ref/dbsnp146.vcf.gz
  - /ref/Mills_and_1000G_gold_standard.indels.vcf.gz
germline_resource: /ref/af-only-gnomad.hg38.vcf.gz
pon: /ref/1000g_pon.hg38.vcf.gz
funcotator_sources: /ref/funcotator_dataSources.v1.7.20200521g
funcotator_ref_version: hg38
outdir: results
```

```bash
nextflow run main.nf -params-file params.yaml -resume
```

## Output layout

```
results/
├── bqsr/<patient_id>/            recal tables + recalibrated BAMs
├── mutect2/<patient_id>/         unfiltered + filtered VCFs, stats, pileups
├── manta/<patient_id>/           candidate indels
├── strelka2/<patient_id>/        somatic.snvs/indels + merged VCF
├── funcotator/<patient_id>/      MAF (or VCF) for both callers
└── pipeline_info/                trace/report/timeline
```

## Notes on WES-specific settings

- `-L`/`--callRegions` intervals are applied throughout (BaseRecalibrator,
  Mutect2, GetPileupSummaries, Funcotator) — always restrict to bait/target
  regions ± padding for exome data to keep runtime and false positives down.
- `--exome` flag is set for both Manta and Strelka2, which adjusts their
  depth filters appropriately for targeted (non-WGS) data.
- Consider adding padding (e.g. 100 bp) around your WES targets when
  generating `--intervals` if you haven't already — this catches spanning
  indels near exon boundaries.

