#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input = "$projectDir/assets/samplesheet.raw_fq.csv"
params.outdir = "results"
params.project = "ChIPalign_NF"
params.ref = "/lustre/home/acct-medkkw/medlyb/database/annotation/gatk_ann/hg38/bowtie2index2/Homo_sapiens_assembly38.fasta"

process FASTQC {
    tag "fastqc in $ID"
    publishDir "$params.outdir/fastqc", pattern: "*.{html,zip}", mode: 'copy'

    input:
    tuple val(ID), val(R1), val(R2), val(Type)

    output:
    tuple val(ID), path("*.{html,zip}")

    script:
    if ( Type == "SE" )
        """
        ln -s  $R1 ${ID}.fastq.gz
        fastqc -t $task.cpus ${ID}.fastq.gz
        """
    else
        """
        ln -s  $R1 ${ID}_R1.fastq.gz
        ln -s  $R2 ${ID}_R2.fastq.gz
        fastqc -t $task.cpus ${ID}_R1.fastq.gz ${ID}_R2.fastq.gz
        """
}

process FastpFilter {
    tag "fastp QC in $ID"
    publishDir "$params.outdir/clean", pattern: "${ID}.*.{html,json}", mode: 'copy'

    input:
    tuple val(ID), val(R1), val(R2), val(Type)

    output:
    tuple val(ID), path("${ID}.*.fq.gz"), val(Type), emit: fastq
    tuple val(ID), path("*.{html,json}"), emit: stats

    script:
    if ( Type == "SE" )
        """
        fastp -i $R1 -o ${ID}.clean.fq.gz -h ${ID}.fastp.html -j ${ID}.fastp.json -w $task.cpus
        """
    else 
        """
        fastp -i $R1 -I $R2 -o ${ID}.clean_R1.fq.gz -O ${ID}.clean_R2.fq.gz -h ${ID}.fastp.html -j ${ID}.fastp.json -w $task.cpus
        """
}

process FASTQ_ALIGN {
    tag "bowtie2_align in $ID"
    publishDir "$params.outdir/align", pattern: "${ID}.bowtie2.out", mode: 'copy'

    input:
    tuple val(ID), path(READ), val(Type)

    output:
    tuple val(ID), val(Type), path("${ID}.sort.bam"), path("${ID}.sort.bam.bai"), emit: bam
    tuple val(ID), path("${ID}.bowtie2.out"), emit: stats

    script:
    if ( Type == "SE" )
        """
        bowtie2 -x ${params.ref} -U ${READ[0]} --rg-id ${ID} --rg "PL:ILLUMINA"  --rg "SM:${ID}" -p $task.cpus 2>${ID}.bowtie2.out | samtools view -@ $task.cpus -Sb -o ${ID}.bam -
        samtools sort -@ $task.cpus -o ${ID}.sort.bam ${ID}.bam && rm ${ID}.bam && samtools index ${ID}.sort.bam
        """
    else
        """
        bowtie2 -x ${params.ref} -1 ${READ[0]} -2 ${READ[1]} --rg-id ${ID} --rg "PL:ILLUMINA"  --rg "SM:${ID}" -p $task.cpus 2>${ID}.bowtie2.out | samtools view -@ $task.cpus -Sb -o ${ID}.bam -
        samtools sort -@ $task.cpus -o ${ID}.sort.bam ${ID}.bam && rm ${ID}.bam && samtools index ${ID}.sort.bam
        """
}


process PICARD_RMDUP{
    tag "picard rmdup in $ID"
    publishDir "$params.outdir/bam", pattern: "${ID}*", mode: 'copy'

    input:
    tuple val(ID), val(Type), path(BAM), path(BAI)

    output:
    tuple val(ID), val(Type), path("${ID}.rmdup.bam"), path("${ID}.rmdup.bam.bai"), emit: bam
    tuple val(ID), path("${ID}.rmdup.metrics"), emit: stats

    script:
    """
    gatk --java-options "-Xmx${task.cpus * 4}g -XX:+UseParallelGC" MarkDuplicates -I ${BAM} -O ${ID}.rmdup.bam -M ${ID}.rmdup.metrics --REMOVE_SEQUENCING_DUPLICATES true -ASO coordinate --VALIDATION_STRINGENCY LENIENT --MAX_RECORDS_IN_RAM ${task.cpus * 4 * 500000}
    samtools index ${ID}.rmdup.bam
    """
}

process BAM_TO_BIGWIG{
    tag "bamCoverage in $ID"
    publishDir "$params.outdir/bigwig", pattern: "${ID}*", mode: 'copy'

    input:
    tuple val(ID), val(Type), path(BAM), path(BAI)

    output:
    tuple val(ID), val(Type), path("${ID}.RPGC.bw")

    script:
    if ( Type == "SE" )
        """
        bamCoverage --bam ${BAM} -o ${ID}.RPGC.bw  --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX --extendReads 146 -p $task.cpus
        """
    else
        """
        bamCoverage --bam ${BAM} -o ${ID}.RPGC.bw  --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX --extendReads -p $task.cpus
        """
}

process MultiQC {
    tag "MultiQC"
    publishDir "$params.outdir/reports", pattern: "*", mode: 'copy'

    input:
    val done

    output:
    path "*"

    script:
    """
    multiqc -n ${params.project}.reports ${workflow.launchDir}/${params.outdir} -f
    """

}

workflow {
    log.info """\
        ChIPalign_NF Align
        ==================
        Sample Info : ${params.input}
        Out Dir     : ${params.outdir}
        """
        .stripIndent(true)

    def input_file = file(params.input)
    def sep_char = input_file.name.toLowerCase().endsWith('.csv') ? ',' : '\t'

    channel.fromPath(params.input)
    .splitCsv(header: true, sep: sep_char)
    .map { row -> ["${row.ID}","${row.R1}","${row.R2}", "${row.Type}"] }
    .set { ch_sample }
    //ch_sample.view()
    FASTQC(ch_sample)
    ch_qc = FastpFilter(ch_sample)
    ch_align = FASTQ_ALIGN(ch_qc.fastq)
    ch_rmdup = PICARD_RMDUP(ch_align.bam)
    ch_bigwig = BAM_TO_BIGWIG(ch_rmdup.bam)

    MultiQC(ch_bigwig)
}
