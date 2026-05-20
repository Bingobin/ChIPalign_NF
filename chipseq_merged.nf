#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
Required samplesheet columns:
ID,R1,R2,Layout,PeakMode,ControlID

- Layout: PE or SE
- PeakMode: TF / Histone / NoCtr
- ControlID: matched input sample ID; leave blank for NoCtr or control-only rows
*/

params.input = "$projectDir/assets/samplesheet.chipseq.csv"
params.outdir = "results/ChIPalign_NF"
params.project = "ChIPalign_NF"
params.ref = "/lustre/home/acct-medkkw/medlyb/database/annotation/gatk_ann/hg38/bowtie2index2/Homo_sapiens_assembly38.fasta"
params.genome = "hg38"
params.effective_genome_size = 2913022398
params.balance_bam = false
params.balance_pairs = 10000000
params.run_motif = true

def normValue(v) {
    def s = v == null ? '' : v.toString().trim()
    return (s in ['', 'NA', 'na', 'null', 'NULL', 'None']) ? '' : s
}

process FASTQC {
    tag "fastqc in $ID"
    publishDir "$params.outdir/fastqc", pattern: "*.{html,zip}", mode: 'copy'

    input:
    tuple val(ID), val(R1), val(R2), val(LAYOUT)

    output:
    tuple val(ID), path("*.{html,zip}")

    script:
    if ( LAYOUT == "SE" )
        """
        ln -s $R1 ${ID}.fastq.gz
        fastqc -t $task.cpus ${ID}.fastq.gz
        """
    else
        """
        ln -s $R1 ${ID}_R1.fastq.gz
        ln -s $R2 ${ID}_R2.fastq.gz
        fastqc -t $task.cpus ${ID}_R1.fastq.gz ${ID}_R2.fastq.gz
        """
}

process FastpFilter {
    tag "fastp QC in $ID"
    publishDir "$params.outdir/clean", pattern: "${ID}.*.{html,json}", mode: 'copy'

    input:
    tuple val(ID), val(R1), val(R2), val(LAYOUT)

    output:
    tuple val(ID), path("${ID}.*.fq.gz"), val(LAYOUT), emit: fastq
    tuple val(ID), path("*.{html,json}"), emit: stats

    script:
    if ( LAYOUT == "SE" )
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
    tuple val(ID), path(READ), val(LAYOUT)

    output:
    tuple val(ID), val(LAYOUT), path("${ID}.sort.bam"), path("${ID}.sort.bam.bai"), emit: bam
    tuple val(ID), path("${ID}.bowtie2.out"), emit: stats

    script:
    if ( LAYOUT == "SE" )
        """
        bowtie2 -x ${params.ref} -U ${READ[0]} --rg-id ${ID} --rg "PL:ILLUMINA" --rg "SM:${ID}" -p $task.cpus 2>${ID}.bowtie2.out | samtools view -@ $task.cpus -Sb -o ${ID}.bam -
        samtools sort -@ $task.cpus -o ${ID}.sort.bam ${ID}.bam
        rm ${ID}.bam
        samtools index ${ID}.sort.bam
        """
    else
        """
        bowtie2 -x ${params.ref} -1 ${READ[0]} -2 ${READ[1]} --rg-id ${ID} --rg "PL:ILLUMINA" --rg "SM:${ID}" -p $task.cpus 2>${ID}.bowtie2.out | samtools view -@ $task.cpus -Sb -o ${ID}.bam -
        samtools sort -@ $task.cpus -o ${ID}.sort.bam ${ID}.bam
        rm ${ID}.bam
        samtools index ${ID}.sort.bam
        """
}

process PICARD_RMDUP {
    tag "picard rmdup in $ID"
    publishDir "$params.outdir/bam", pattern: "${ID}*", mode: 'copy'

    input:
    tuple val(ID), val(LAYOUT), path(BAM), path(BAI)

    output:
    tuple val(ID), val(LAYOUT), path("${ID}.rmdup.bam"), path("${ID}.rmdup.bam.bai"), emit: bam
    tuple val(ID), path("${ID}.rmdup.metrics"), emit: stats

    script:
    """
    gatk --java-options "-Xmx${task.cpus * 4}g -XX:+UseParallelGC" MarkDuplicates -I ${BAM} -O ${ID}.rmdup.bam -M ${ID}.rmdup.metrics --REMOVE_SEQUENCING_DUPLICATES true -ASO coordinate --VALIDATION_STRINGENCY LENIENT --MAX_RECORDS_IN_RAM ${task.cpus * 4 * 500000}
    samtools index ${ID}.rmdup.bam
    """
}

process BAM_TO_BIGWIG {
    tag "bamCoverage in $ID"
    publishDir "$params.outdir/bigwig", pattern: "${ID}*", mode: 'copy'

    input:
    tuple val(ID), val(LAYOUT), path(BAM), path(BAI)

    output:
    tuple val(ID), val(LAYOUT), path("${ID}.RPGC.bw")

    script:
    if ( LAYOUT == "SE" )
        """
        bamCoverage --bam ${BAM} -o ${ID}.RPGC.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize ${params.effective_genome_size} --ignoreForNormalization chrX --extendReads 146 -p $task.cpus
        """
    else
        """
        bamCoverage --bam ${BAM} -o ${ID}.RPGC.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize ${params.effective_genome_size} --ignoreForNormalization chrX --extendReads -p $task.cpus
        """
}

process BamPairBalancerWithCtrl {
    tag "BamPairBalancer in $ID"
    publishDir "$params.outdir/bam_balance", pattern: "${ID}*", mode: 'copy'

    input:
    tuple val(ID), path(BAM_t), path(BAM_c), val(PEAK_MODE)

    output:
    tuple val(ID), path("${ID}.treat.bal.sort.bam"), path("${ID}.ctrl.bal.sort.bam"), val(PEAK_MODE)

    script:
    """
    BamPairBalancer/bam_pair_balancer.py -i ${BAM_t} -o ${ID}.treat.bal.bam -n ${params.balance_pairs} -p $task.cpus
    samtools sort -@ $task.cpus -o ${ID}.treat.bal.sort.bam ${ID}.treat.bal.bam
    samtools index ${ID}.treat.bal.sort.bam

    BamPairBalancer/bam_pair_balancer.py -i ${BAM_c} -o ${ID}.ctrl.bal.bam -n ${params.balance_pairs} -p $task.cpus
    samtools sort -@ $task.cpus -o ${ID}.ctrl.bal.sort.bam ${ID}.ctrl.bal.bam
    samtools index ${ID}.ctrl.bal.sort.bam
    """
}

process MACS2_callpeaks_withCtrl {
    tag "macs2 callpeaks in $ID"
    publishDir "$params.outdir/callpeaks", pattern: "${ID}*", mode: 'copy'

    input:
    tuple val(ID), path(BAM_t), path(BAM_c), val(PEAK_MODE)

    output:
    tuple val(ID), path("${ID}_peaks.xls"), emit: xls
    tuple val(ID), path("${ID}_peaks.narrowPeak"), emit: peaks
    tuple val(ID), path("${ID}_summits.bed"), emit: summits

    script:
    if ( PEAK_MODE == "Histone" )
        """
        macs2 callpeak -t $BAM_t -c $BAM_c -f BAMPE -n ${ID} -p 1e-9
        """
    else
        """
        macs2 callpeak -t $BAM_t -c $BAM_c -f BAMPE -n ${ID} -q 0.01
        """
}

process MACS2_callpeaks_noCtrl {
    tag "macs2 callpeaks in $ID"
    publishDir "$params.outdir/callpeaks", pattern: "${ID}*", mode: 'copy'

    input:
    tuple val(ID), path(BAM_t), val(PEAK_MODE)

    output:
    tuple val(ID), path("${ID}_peaks.xls"), emit: xls
    tuple val(ID), path("${ID}_peaks.narrowPeak"), emit: peaks
    tuple val(ID), path("${ID}_summits.bed"), emit: summits

    script:
    """
    macs2 callpeak -t $BAM_t -f BAMPE -n ${ID} --SPMR -q 0.01 --keep-dup 1 --extsize=250 --nomodel -g hs
    """
}

process HOMER_annotatePeaks {
    tag "homer annotatePeaks in $ID"
    publishDir "$params.outdir/annoPeaks", pattern: "${ID}_peaks.narrowPeak.anno.*", mode: 'copy'

    input:
    tuple val(ID), path(PEAK)

    output:
    tuple val(ID), path("${ID}_peaks.narrowPeak.anno.*")

    script:
    """
    annotatePeaks.pl ${PEAK} ${params.genome} -annStats ${ID}_peaks.narrowPeak.anno.stats > ${ID}_peaks.narrowPeak.anno.txt
    """
}

process HOMER_findMotifs {
    tag "homer findMotifs in $ID"
    publishDir "$params.outdir/motif", pattern: "${ID}.MotifOutput_j", mode: 'copy'

    input:
    tuple val(ID), path(SUMMIT)

    output:
    tuple val(ID), path("${ID}.MotifOutput_j")

    script:
    """
    findMotifsGenome.pl ${SUMMIT} ${params.genome} ${ID}.MotifOutput_j -size 200 -mask -p $task.cpus -len 6,8,12,16,20 -mknown /lustre/home/acct-medkkw/medlyb/wl_proj/WL234_Lib/database/JASPAR/JASPAR2020_CORE_vertebrates_non-redundant_pfms_homer.txt
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
    multiqc -n ${params.project}.fastqc.reports ${workflow.launchDir}/${params.outdir}/fastqc/* -f
    multiqc -n ${params.project}.clean.reports ${workflow.launchDir}/${params.outdir}/clean/* -f
    multiqc -n ${params.project}.align.reports ${workflow.launchDir}/${params.outdir}/align/* -f
    multiqc -n ${params.project}.rmdup.reports ${workflow.launchDir}/${params.outdir}/bam/* -f
    """
}

workflow ALIGNMENT {
    take:
    ch_sample

    main:
    FASTQC(ch_sample)
    ch_qc = FastpFilter(ch_sample)
    ch_align = FASTQ_ALIGN(ch_qc.fastq)
    ch_rmdup = PICARD_RMDUP(ch_align.bam)
    ch_bigwig = BAM_TO_BIGWIG(ch_rmdup.bam)

    emit:
    bam = ch_rmdup.bam
    bigwig = ch_bigwig
}

workflow PEAK_CALLING {
    take:
    ch_peak_with_ctrl
    ch_peak_no_ctrl

    main:
    def ch_with_ctrl_for_macs2 = params.balance_bam ? BamPairBalancerWithCtrl(ch_peak_with_ctrl) : ch_peak_with_ctrl

    ch_callpeak_with_ctrl = MACS2_callpeaks_withCtrl(ch_with_ctrl_for_macs2)
    ch_callpeak_no_ctrl = MACS2_callpeaks_noCtrl(ch_peak_no_ctrl)

    ch_peak_files = ch_callpeak_with_ctrl.peaks.mix(ch_callpeak_no_ctrl.peaks)
    ch_summit_files = ch_callpeak_with_ctrl.summits.mix(ch_callpeak_no_ctrl.summits)

    HOMER_annotatePeaks(ch_peak_files)
    if ( params.run_motif ) {
        HOMER_findMotifs(ch_summit_files)
    }

    emit:
    peaks = ch_peak_files
}

workflow {
    log.info """\
        ChIPalign_NF
        ============
        Sample Info   : ${params.input}
        Out Dir       : ${params.outdir}
        Reference     : ${params.ref}
        Genome        : ${params.genome}
        balance_bam   : ${params.balance_bam}
        balance_pairs : ${params.balance_pairs}
        run_motif     : ${params.run_motif}
        """.stripIndent(true)

    def input_file = file(params.input)
    def sep_char = input_file.name.toLowerCase().endsWith('.csv') ? ',' : '\t'

    channel
        .fromPath(params.input)
        .splitCsv(header: true, sep: sep_char)
        .map { row ->
            def id = normValue(row.ID)
            def r1 = normValue(row.R1)
            def r2 = normValue(row.R2)
            def layout = normValue(row.Layout ?: row.TYPE ?: row.Type).toUpperCase()
            def peakMode = normValue(row.PeakMode)
            def controlId = normValue(row.ControlID)
            tuple(id, r1, r2, layout, peakMode, controlId)
        }
        .set { ch_sheet }

    ch_align_input = ch_sheet.map { id, r1, r2, layout, peakMode, controlId ->
        tuple(id, r1, r2, layout)
    }

    ch_peak_meta = ch_sheet
        .filter { id, r1, r2, layout, peakMode, controlId -> peakMode }
        .map { id, r1, r2, layout, peakMode, controlId ->
            tuple(id, peakMode, controlId)
        }

    ch_aligned_bam = ALIGNMENT(ch_align_input).bam

    ch_treat_with_meta = ch_peak_meta
        .join(ch_aligned_bam.map { id, layout, bam, bai -> tuple(id, bam) }, by: 0)
        .map { id, peakMode, controlId, bam ->
            tuple(id, bam, peakMode, controlId)
        }

    ch_control_bam = ch_aligned_bam
        .map { id, layout, bam, bai -> tuple(id, bam) }

    ch_peak_no_ctrl = ch_treat_with_meta
        .filter { id, bam, peakMode, controlId -> !controlId || peakMode == 'NoCtr' }
        .map { id, bam, peakMode, controlId -> tuple(id, bam, peakMode ?: 'NoCtr') }

    ch_peak_with_ctrl = ch_treat_with_meta
        .filter { id, bam, peakMode, controlId -> controlId && peakMode != 'NoCtr' }
        .combine(ch_control_bam)
        .filter { id, bam_t, peakMode, controlId, ctrlId, bam_c -> controlId == ctrlId }
        .map { id, bam_t, peakMode, controlId, ctrlId, bam_c ->
            tuple(id, bam_t, bam_c, peakMode)
        }

    PEAK_CALLING(ch_peak_with_ctrl, ch_peak_no_ctrl)

    MultiQC(ch_aligned_bam.map { ignored -> true }.first())

    workflow.onComplete = {
        log.info(workflow.success ? "\nDone! See results --> $params.outdir\n" : "Oops.. something went wrong")
    }
}
