# ChIPalign_NF

ChIPalign_NF is an end-to-end ChIP-seq Nextflow pipeline for read QC, alignment, duplicate removal, bigWig generation, peak calling, peak annotation, and motif analysis.

The current main entry point is [chipseq_merged.nf](/Users/liuyabin/Desktop/ChIPalign_NF/chipseq_merged.nf). It combines the two original sub-pipelines:

- [chip_align.nf](/Users/liuyabin/Desktop/ChIPalign_NF/chip_align.nf): FASTQ to BAM and bigWig
- [chip_callpeak.nf](/Users/liuyabin/Desktop/ChIPalign_NF/chip_callpeak.nf): BAM to peak calling and motif analysis

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

The merged pipeline expects a CSV or TSV samplesheet with the following columns:

`ID,R1,R2,Layout,PeakMode,ControlID`

Column definitions:

- `ID`: unique sample name
- `R1`: read 1 FASTQ path
- `R2`: read 2 FASTQ path; leave blank for single-end data
- `Layout`: `PE` or `SE`
- `PeakMode`: `TF`, `Histone`, or `NoCtr`
- `ControlID`: matched input/control sample ID; leave blank for `NoCtr` or control-only rows

Example:

```csv
ID,R1,R2,Layout,PeakMode,ControlID
Input_1,/path/to/Input_1_R1.fastq.gz,/path/to/Input_1_R2.fastq.gz,PE,,
CTCF_1,/path/to/CTCF_1_R1.fastq.gz,/path/to/CTCF_1_R2.fastq.gz,PE,TF,Input_1
H3K27ac_1,/path/to/H3K27ac_1_R1.fastq.gz,/path/to/H3K27ac_1_R2.fastq.gz,PE,Histone,Input_1
CUTRUN_1,/path/to/CUTRUN_1_R1.fastq.gz,/path/to/CUTRUN_1_R2.fastq.gz,PE,NoCtr,
```

Reference file in this repository:

- [samplesheet.chipseq.csv](/Users/liuyabin/Desktop/ChIPalign_NF/assets/samplesheet.chipseq.csv)

## Important Metadata Rules

- `Layout` and `PeakMode` are different concepts and should not be merged into one column.
- `Layout` only describes read structure: `SE` or `PE`.
- `PeakMode` only describes the peak-calling strategy.
- Samples used only as controls can keep `PeakMode` and `ControlID` empty.
- If a treatment sample uses a control, `ControlID` must exactly match the control sample `ID`.
- For `NoCtr`, leave `ControlID` empty.

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
  --outdir results/ChIPalign_NF
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

Key parameters currently exposed in [chipseq_merged.nf](/Users/liuyabin/Desktop/ChIPalign_NF/chipseq_merged.nf):

- `--input`: input samplesheet path
- `--outdir`: output directory, default `results/ChIPalign_NF`
- `--project`: project name used in report naming
- `--ref`: Bowtie2 reference index prefix
- `--genome`: genome label for HOMER, default `hg38`
- `--effective_genome_size`: genome size for RPGC normalization
- `--balance_bam`: whether to run BAM balancing before peak calling
- `--balance_pairs`: number of read pairs for BAM balancing
- `--run_motif`: whether to run HOMER motif discovery

Default process resources are currently defined in [nextflow.config](/Users/liuyabin/Desktop/ChIPalign_NF/nextflow.config):

- `cpus = 4`
- `memory = 4 GB x cpus`

## Output Structure

By default, results are written to:

`results/ChIPalign_NF`

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

- `PeakMode = TF`: `macs2 callpeak -q 0.01` with matched control
- `PeakMode = Histone`: `macs2 callpeak -p 1e-9` with matched control
- `PeakMode = NoCtr`: control-free MACS2 mode with `--SPMR --keep-dup 1 --extsize=250 --nomodel -g hs`

This means your samplesheet drives downstream behavior directly. If the metadata are wrong, peak calling will also be wrong.

## Notes and Current Limitations

- `--ref` is passed to `bowtie2 -x`, so it should be a Bowtie2 index prefix, not a raw FASTA file path.
- The current MACS2 outputs are defined as `narrowPeak` and `summits.bed` for all modes. That may not be ideal for broad histone marks if you later switch to a true broad peak strategy.
- The motif database path in `HOMER_findMotifs` is currently hard-coded in the pipeline.
- The BAM balancer script path is currently hard-coded as `BamPairBalancer/bam_pair_balancer.py`.
- The pipeline has not yet been packaged with profiles such as `standard`, `slurm`, `docker`, or `conda`.
- The current environment in this session does not have `nextflow` installed, so this repository update was done by static inspection rather than an actual pipeline run.

## Repository Files

- [chipseq_merged.nf](/Users/liuyabin/Desktop/ChIPalign_NF/chipseq_merged.nf): main merged pipeline
- [chip_align.nf](/Users/liuyabin/Desktop/ChIPalign_NF/chip_align.nf): original alignment-only workflow
- [chip_callpeak.nf](/Users/liuyabin/Desktop/ChIPalign_NF/chip_callpeak.nf): original peak-calling workflow
- [nextflow.config](/Users/liuyabin/Desktop/ChIPalign_NF/nextflow.config): default process settings and Nextflow manifest
- [samplesheet.chipseq.csv](/Users/liuyabin/Desktop/ChIPalign_NF/assets/samplesheet.chipseq.csv): example samplesheet

## Suggested Next Improvements

- add `profiles` for local and cluster execution
- add Conda or container support
- move hard-coded database paths into parameters
- split the merged pipeline into `main.nf`, `modules/`, and `subworkflows/`
- add a real test dataset and a minimal smoke test
