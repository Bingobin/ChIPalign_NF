# ChIPalign_NF

ChIPalign_NF is an end-to-end ChIP-seq Nextflow pipeline for read QC, alignment, duplicate removal, bigWig generation, peak calling, peak annotation, and motif analysis.

The current main entry point is [chipseq_merged.nf](chipseq_merged.nf). It combines the two original sub-pipelines:

- [chip_align.nf](chip_align.nf): FASTQ to BAM and bigWig
- [chip_callpeak.nf](chip_callpeak.nf): BAM to peak calling and motif analysis

## Workflow Overview

The merged pipeline runs the following steps:

1. `FASTQC`
2. `fastp`
3. `bowtie2` alignment
4. `samtools sort/index`
5. `GATK MarkDuplicates`
6. `bamCoverage` to generate RPGC-normalized bigWig
7. Optional `BamPairBalancer`
8. `MACS2 callpeak`
9. `HOMER annotatePeaks`
10. Optional `HOMER findMotifsGenome`
11. `MultiQC`

The pipeline is designed for two common ChIP-seq scenarios:

- factor ChIP-seq with a matched input control
- no-control peak calling workflows such as `NoCtr`

## Input Samplesheet

The merged pipeline supports two input modes:

- `--input_mode fastq`: start from FASTQ files and run QC, alignment, duplicate removal, bigWig generation, and peak calling
- `--input_mode bam`: start from existing BAM files and run peak calling, peak annotation, motif analysis, and MultiQC

FASTQ mode samplesheet columns:

`ID,R1,R2,Layout,PeakMode,ControlID,ControlBam`

BAM mode samplesheet columns:

`ID,BAM,PeakMode,ControlID,ControlBam`

Column definitions:

- `ID`: unique sample name
- `R1`: read 1 FASTQ path
- `R2`: read 2 FASTQ path; leave blank for single-end data
- `Layout`: `PE` or `SE`
- `BAM`: treatment/input BAM path for BAM input mode
- `PeakMode`: `TF`, `Histone`, or `NoCtr`
- `ControlID`: matched input/control sample ID; leave blank for `NoCtr` or control-only rows
- `ControlBam`: matched input/control BAM path; used when `ControlID` is blank

FASTQ mode example:

```csv
ID,R1,R2,Layout,PeakMode,ControlID,ControlBam
Input_1,/path/to/Input_1_R1.fastq.gz,/path/to/Input_1_R2.fastq.gz,PE,,,
CTCF_1,/path/to/CTCF_1_R1.fastq.gz,/path/to/CTCF_1_R2.fastq.gz,PE,TF,Input_1,
H3K27ac_1,/path/to/H3K27ac_1_R1.fastq.gz,/path/to/H3K27ac_1_R2.fastq.gz,PE,Histone,Input_1,
CTCF_2,/path/to/CTCF_2_R1.fastq.gz,/path/to/CTCF_2_R2.fastq.gz,PE,TF,,/path/to/Input_2.rmdup.bam
CUTRUN_1,/path/to/CUTRUN_1_R1.fastq.gz,/path/to/CUTRUN_1_R2.fastq.gz,PE,NoCtr,,
```

BAM mode example:

```csv
ID,BAM,PeakMode,ControlID,ControlBam
Input_1,/path/to/Input_1.rmdup.bam,,,
CTCF_1,/path/to/CTCF_1.rmdup.bam,TF,Input_1,
H3K27ac_1,/path/to/H3K27ac_1.rmdup.bam,Histone,Input_1,
CTCF_2,/path/to/CTCF_2.rmdup.bam,TF,,/path/to/Input_2.rmdup.bam
CUTRUN_1,/path/to/CUTRUN_1.rmdup.bam,NoCtr,,
```

Reference files in this repository:

- [samplesheet.chipseq.csv](assets/samplesheet.chipseq.csv)
- [samplesheet.chipseq.bam.csv](assets/samplesheet.chipseq.bam.csv)

## Important Metadata Rules

- `Layout` and `PeakMode` are different concepts and should not be merged into one column.
- `Layout` only describes read structure: `SE` or `PE`.
- `PeakMode` only describes the peak-calling strategy.
- Samples used only as controls can keep `PeakMode` and `ControlID` empty.
- If a treatment sample uses a control, `ControlID` must exactly match the control sample `ID`.
- If the control BAM already exists, leave `ControlID` empty and set `ControlBam`.
- For `NoCtr`, leave both `ControlID` and `ControlBam` empty.

## Requirements

This pipeline assumes the following software is available in your runtime environment:

- `Nextflow`
- `fastqc`
- `fastp`
- `bowtie2`
- `samtools`
- `gatk`
- `bamCoverage` from deepTools
- `macs2`
- `annotatePeaks.pl` from HOMER
- `findMotifsGenome.pl` from HOMER
- `multiqc`
- `BamPairBalancer/bam_pair_balancer.py` if `--balance_bam true`

The current repository does not yet define a container, Conda environment, or module file. That means the pipeline currently depends on your local HPC or workstation software environment being set correctly.

## Running the Pipeline

Basic run:

```bash
nextflow run chipseq_merged.nf \
  --input assets/samplesheet.chipseq.csv \
  --outdir results
```

Start directly from BAM files:

```bash
nextflow run chipseq_merged.nf \
  --input_mode bam \
  --input assets/samplesheet.chipseq.bam.csv \
  --outdir results
```

Run with motif analysis disabled:

```bash
nextflow run chipseq_merged.nf \
  --input assets/samplesheet.chipseq.csv \
  --run_motif false
```

Run with BAM balancing enabled:

```bash
nextflow run chipseq_merged.nf \
  --input assets/samplesheet.chipseq.csv \
  --balance_bam true \
  --balance_pairs 10000000
```

You can also override the reference and genome settings:

```bash
nextflow run chipseq_merged.nf \
  --input assets/samplesheet.chipseq.csv \
  --ref /path/to/bowtie2/index/prefix \
  --genome hg38 \
  --effective_genome_size 2913022398
```

## Main Parameters

Default parameter values are defined in [nextflow.config](nextflow.config). They can still be overridden on the command line:

- `--input`: input samplesheet path
- `--input_mode`: input mode, `fastq` or `bam`, default `fastq`
- `--outdir`: output directory, default `results`
- `--project`: project name used in report naming
- `--ref`: Bowtie2 reference index prefix
- `--genome`: genome label for HOMER, default `hg38`
- `--effective_genome_size`: genome size for RPGC normalization
- `--balance_bam`: whether to run BAM balancing before peak calling
- `--balance_pairs`: number of read pairs for BAM balancing
- `--run_motif`: whether to run HOMER motif discovery, default `false`
- `--macs2_tf_qvalue`: MACS2 `-q` value for TF mode, default `0.01`
- `--macs2_histone_pvalue`: MACS2 `-p` value for Histone mode, default `1e-9`
- `--macs2_noctrl_qvalue`: MACS2 `-q` value for NoCtr mode, default `0.01`
- `--macs2_noctrl_keep_dup`: MACS2 `--keep-dup` value for NoCtr mode, default `1`
- `--macs2_noctrl_extsize`: MACS2 `--extsize` value for NoCtr mode, default `250`
- `--motif_size`: HOMER `findMotifsGenome.pl -size` value, default `200`
- `--motif_len`: HOMER `findMotifsGenome.pl -len` value, default `6,8,12,16,20`
- `--motif_mknown`: HOMER `findMotifsGenome.pl -mknown` motif database path

Default process resources are currently defined in [nextflow.config](nextflow.config):

- `cpus = 4`
- `memory = 4 GB x cpus`

## Output Structure

By default, results are written to:

`results`

Main output subdirectories:

- `fastqc/`: raw FASTQ QC reports
- `clean/`: fastp reports and cleaned FASTQ files
- `align/`: Bowtie2 alignment logs
- `bam/`: duplicate-removed BAM files and metrics
- `bigwig/`: normalized bigWig tracks
- `bam_balance/`: balanced BAMs if enabled
- `callpeaks/`: MACS2 peak calling results
- `annoPeaks/`: HOMER peak annotation results
- `motif/`: HOMER motif results if enabled
- `reports/`: MultiQC reports

## Peak Calling Behavior

The pipeline currently uses two peak-calling branches:

- `PeakMode = TF`: `macs2 callpeak -q`, controlled by `--macs2_tf_qvalue`, default `0.01`
- `PeakMode = Histone`: `macs2 callpeak -p`, controlled by `--macs2_histone_pvalue`, default `1e-9`
- `PeakMode = NoCtr`: control-free MACS2 mode with `--SPMR`, `--nomodel`, and configurable `--macs2_noctrl_qvalue`, `--macs2_noctrl_keep_dup`, and `--macs2_noctrl_extsize`

This means your samplesheet drives downstream behavior directly. If the metadata are wrong, peak calling will also be wrong.

## Notes and Current Limitations

- `--ref` is passed to `bowtie2 -x`, so it should be a Bowtie2 index prefix, not a raw FASTA file path.
- The current MACS2 outputs are defined as `narrowPeak` and `summits.bed` for all modes. That may not be ideal for broad histone marks if you later switch to a true broad peak strategy.
- The motif database path can be changed with `--motif_mknown`.
- The BAM balancer script path is currently hard-coded as `BamPairBalancer/bam_pair_balancer.py`.
- The pipeline has not yet been packaged with profiles such as `standard`, `slurm`, `docker`, or `conda`.
- The current environment in this session does not have `nextflow` installed, so this repository update was done by static inspection rather than an actual pipeline run.

## Repository Files

- [chipseq_merged.nf](chipseq_merged.nf): main merged pipeline
- [chip_align.nf](chip_align.nf): original alignment-only workflow
- [chip_callpeak.nf](chip_callpeak.nf): original peak-calling workflow
- [nextflow.config](nextflow.config): default process settings and Nextflow manifest
- [samplesheet.chipseq.csv](assets/samplesheet.chipseq.csv): example samplesheet
- [samplesheet.chipseq.bam.csv](assets/samplesheet.chipseq.bam.csv): example BAM-mode samplesheet
