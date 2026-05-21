#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
Required samplesheet columns depend on --input_mode:

- For --input_mode fastq: ID,R1,R2,Layout,PeakMode,ControlID,ControlBam
- For --input_mode bam: ID,BAM,PeakMode,ControlID,ControlBam
- Layout: PE or SE
- BAM: treatment/input BAM path for BAM input mode
- PeakMode: TF / Histone / NoCtr
- ControlID: matched input sample ID; leave blank for NoCtr or control-only rows
- ControlBam: matched input/control BAM path; used when ControlID is blank
*/

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
    publishDir "$params.outdir/clean", pattern: "*.{html,json}", mode: 'copy'

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
    publishDir "$params.outdir/align", pattern: "*.bowtie2.out", mode: 'copy'

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
    publishDir "$params.outdir/bam", pattern: "*", mode: 'copy'

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
    publishDir "$params.outdir/bigwig", pattern: "*.bw", mode: 'copy'

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
    publishDir "$params.outdir/bam_balance", pattern: "*", mode: 'copy'

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
    cpus 2
    publishDir "$params.outdir/callpeaks", pattern: "*", mode: 'copy'

    input:
    tuple val(ID), path(BAM_t), path(BAM_c), val(PEAK_MODE)

    output:
    tuple val(ID), path("${ID}_peaks.xls"), emit: xls
    tuple val(ID), path("${ID}_peaks.narrowPeak"), emit: peaks
    tuple val(ID), path("${ID}_summits.bed"), emit: summits

    script:
    if ( PEAK_MODE == "Histone" )
        """
        macs2 callpeak -t $BAM_t -c $BAM_c -f BAMPE -n ${ID} -p ${params.macs2_histone_pvalue}
        """
    else
        """
        macs2 callpeak -t $BAM_t -c $BAM_c -f BAMPE -n ${ID} -q ${params.macs2_tf_qvalue}
        """
}

process MACS2_callpeaks_noCtrl {
    tag "macs2 callpeaks in $ID"
    cpus 2
    publishDir "$params.outdir/callpeaks", pattern: "*", mode: 'copy'

    input:
    tuple val(ID), path(BAM_t), val(PEAK_MODE)

    output:
    tuple val(ID), path("${ID}_peaks.xls"), emit: xls
    tuple val(ID), path("${ID}_peaks.narrowPeak"), emit: peaks
    tuple val(ID), path("${ID}_summits.bed"), emit: summits

    script:
    """
    macs2 callpeak -t $BAM_t -f BAMPE -n ${ID} --SPMR -q ${params.macs2_noctrl_qvalue} --keep-dup ${params.macs2_noctrl_keep_dup} --extsize=${params.macs2_noctrl_extsize} --nomodel -g hs
    """
}

process HOMER_annotatePeaks {
    tag "homer annotatePeaks in $ID"
    cpus 2
    publishDir "$params.outdir/annoPeaks", pattern: "*.anno.*", mode: 'copy'

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
    publishDir "$params.outdir/motif", pattern: "*.MotifOutput_j", mode: 'copy'

    input:
    tuple val(ID), path(SUMMIT)

    output:
    tuple val(ID), path("${ID}.MotifOutput_j")

    script:
    """
    findMotifsGenome.pl ${SUMMIT} ${params.genome} ${ID}.MotifOutput_j -size ${params.motif_size} -mask -p $task.cpus -len ${params.motif_len} -mknown ${params.motif_mknown}
    """
}

process MultiQC {
    tag "MultiQC"
    cpus 2
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

    ch_peak_files = channel.empty()
        .mix(ch_callpeak_with_ctrl.peaks)
        .mix(ch_callpeak_no_ctrl.peaks)
    ch_summit_files = channel.empty()
        .mix(ch_callpeak_with_ctrl.summits)
        .mix(ch_callpeak_no_ctrl.summits)

    HOMER_annotatePeaks(ch_peak_files)
    if ( params.run_motif ) {
        HOMER_findMotifs(ch_summit_files)
    }

    emit:
    peaks = ch_peak_files
}

workflow {
    def input_file = file(params.input)
    def sep_char = input_file.name.toLowerCase().endsWith('.csv') ? ',' : '\t'
    def input_mode = normValue(params.input_mode).toLowerCase()

    if ( !(input_mode in ['fastq', 'bam']) ) {
        error "Invalid --input_mode '${params.input_mode}'. Use 'fastq' or 'bam'."
    }

    log.info """\
        ChIPalign_NF
        ============
        Sample Info   : ${params.input}
        Input Mode    : ${input_mode}
        Out Dir       : ${params.outdir}
        Reference     : ${params.ref}
        Genome        : ${params.genome}
        balance_bam   : ${params.balance_bam}
        balance_pairs : ${params.balance_pairs}
        run_motif     : ${params.run_motif}
        """.stripIndent(true)

    channel
        .fromPath(params.input)
        .splitCsv(header: true, sep: sep_char)
        .map { row ->
            def id = normValue(row.ID)
            def r1 = normValue(row.R1)
            def r2 = normValue(row.R2)
            def layout = normValue(row.Layout ?: row.TYPE ?: row.Type).toUpperCase()
            def bam = normValue(row.BAM ?: row.Bam ?: row.bam)
            def peakMode = normValue(row.PeakMode)
            def controlId = normValue(row.ControlID)
            def controlBam = normValue(row.ControlBam)
            tuple(id, r1, r2, layout, bam, peakMode, controlId, controlBam)
        }
        .set { ch_sheet }

    ch_peak_meta = ch_sheet
        .filter { id, r1, r2, layout, bam, peakMode, controlId, controlBam -> peakMode }
        .map { id, r1, r2, layout, bam, peakMode, controlId, controlBam ->
            tuple(id, peakMode, controlId, controlBam)
        }

    def ch_bam_for_peak

    if ( input_mode == 'fastq' ) {
        ch_align_input = ch_sheet.map { id, r1, r2, layout, bam, peakMode, controlId, controlBam ->
            tuple(id, r1, r2, layout)
        }

        ch_aligned_bam = ALIGNMENT(ch_align_input).bam

        ch_bam_for_peak = ch_aligned_bam
            .map { id, layout, bam, bai -> tuple(id, bam) }
    } else {
        ch_bam_for_peak = ch_sheet
            .map { id, r1, r2, layout, bam, peakMode, controlId, controlBam ->
                tuple(id, file(bam))
            }
    }

    ch_treat_with_meta = ch_peak_meta
        .join(ch_bam_for_peak, by: 0)
        .map { id, peakMode, controlId, controlBam, bam ->
            tuple(id, bam, peakMode, controlId, controlBam)
        }

    ch_control_bam = ch_bam_for_peak

    ch_peak_no_ctrl = ch_treat_with_meta
        .filter { id, bam, peakMode, controlId, controlBam -> peakMode == 'NoCtr' || (!controlId && !controlBam) }
        .map { id, bam, peakMode, controlId, controlBam -> tuple(id, bam, peakMode ?: 'NoCtr') }

    ch_peak_with_ctrl_by_id = ch_treat_with_meta
        .filter { id, bam, peakMode, controlId, controlBam -> controlId && peakMode != 'NoCtr' }
        .combine(ch_control_bam)
        .filter { id, bam_t, peakMode, controlId, controlBam, ctrlId, bam_c -> controlId == ctrlId }
        .map { id, bam_t, peakMode, controlId, controlBam, ctrlId, bam_c ->
            tuple(id, bam_t, bam_c, peakMode)
        }

    ch_peak_with_ctrl_by_bam = ch_treat_with_meta
        .filter { id, bam, peakMode, controlId, controlBam -> !controlId && controlBam && peakMode != 'NoCtr' }
        .map { id, bam, peakMode, controlId, controlBam ->
            tuple(id, bam, file(controlBam), peakMode)
        }

    ch_peak_with_ctrl = ch_peak_with_ctrl_by_id.mix(ch_peak_with_ctrl_by_bam)

    ch_peaks = PEAK_CALLING(ch_peak_with_ctrl, ch_peak_no_ctrl).peaks

    MultiQC(ch_peaks.collect().map { ignored -> true })

    workflow.onComplete = {
        log.info(workflow.success ? "\nDone! See results --> $params.outdir\n" : "Oops.. something went wrong")
    }
}
