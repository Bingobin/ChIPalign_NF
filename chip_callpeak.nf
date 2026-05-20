#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input = "$projectDir/assets/samplesheet.bam.csv"
params.outdir = "results/ChIPalign_NF"
params.balance_bam = false
params.balance_pairs = 10000000


process BamPairBalancer {
    tag "BamPairBalancer in $ID"
    publishDir "$params.outdir/bam", pattern: "${ID}*", mode: 'copy'

    input:
    tuple val(ID), path(BAM_t), path(BAM_c), val(Type)

    output:
    tuple val(ID),
          path("${BAM_t.baseName}.bal.sort.bam"),
          path("${BAM_t.baseName}.bal.sort.bam.bai"),
          path("${BAM_t.baseName}.bal.bw"),
          path("${BAM_c ? BAM_c.baseName + '.bal.sort.bam' : null}", optional: true),
          path("${BAM_c ? BAM_c.baseName + '.bal.sort.bam.bai' : null}", optional: true),
          path("${BAM_c ? BAM_c.baseName + '.bal.bw' : null}", optional: true),
          val(Type)

    script:
    def N = params.balance_pairs

    """
    BamPairBalancer/bam_pair_balancer.py -i ${BAM_t} -o ${BAM_t.baseName}.bal.bam -n 10000000 -p $task.cpus
    samtools sort -@ $task.cpus -i ${BAM_t.baseName}.bal.bam -o  ${BAM_t.baseName}.bal.sort.bam
    samtools index ${BAM_t.baseName}.bal.sort.bam
    bamCoverage --bam ${BAM_t.baseName}.bal.sort.bam -o ${BAM_t.baseName}.bal.bw  --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX --extendReads -p $task.cpus

    ${ BAM_c ?
        """
        BamPairBalancer/bam_pair_balancer.py -i ${BAM_c} -o ${BAM_c.baseName}.bal.bam -n 10000000 -p $task.cpus
        samtools sort -@ $task.cpus -i ${BAM_c.baseName}.bal.bam -o  ${BAM_c.baseName}.bal.sort.bam
        samtools index ${BAM_c.baseName}.bal.sort.bam
        bamCoverage --bam ${BAM_t.baseName}.bal.sort.bam -o ${BAM_t.baseName}.bal.bw  --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX --extendReads -p $task.cpus
        """ : ""}
    """
}

process MACS2_callpeaks {
    tag "macs2 callpeaks in $ID"
    publishDir "$params.outdir/callpeaks", pattern: "${ID}*", mode: 'copy'

    input:
    tuple val(ID), path(BAM_t), path(BAM_c), val(Type)

    output:
    tuple val(ID), path("${ID}_peaks.xls"), emit: xls
//    tuple val(ID), path("${ID}_peaks.${Type == 'Histone' ? 'broadPeak' : 'narrowPeak'}"), emit: peaks
//    tuple val(ID), path("${ID}_${Type == 'Histone' ? 'peaks.gappedPeak' : 'summits.bed'}"), emit: summits
    tuple val(ID), path("${ID}_peaks.narrowPeak"), emit: peaks
    tuple val(ID), path("${ID}_summits.bed"), emit: summits

    script:
    if( Type == "Histone" )
        """
        macs2 callpeak -t $BAM_t -c $BAM_c -f BAMPE -n ${ID} -p 1e-9 #--broad --broad-cutoff 1e-4 
        """
    else if ( Type == "NoCtr" )
        """
        macs2 callpeak -t $BAM_t -f BAMPE  -n ${ID} --SPMR -q 0.01 --keep-dup 1  --extsize=250 --nomodel -g hs 
        """
    else 
        """
        macs2 callpeak -t $BAM_t -c $BAM_c -f BAMPE -n ${ID} -q 0.01
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
    annotatePeaks.pl ${PEAK} hg38 -annStats ${ID}_peaks.narrowPeak.anno.stats > ${ID}_peaks.narrowPeak.anno.txt
    """
}

process HOMER_findMotifs {
    tag "homer findMotifs in $ID"
    publishDir "$params.outdir/motif", pattern: "${ID}.MotifOutput_j", mode: 'copy'

    input:
    tuple val(ID), path(SUMMIT)

    output:
    tuple val(ID), path("${ID}.MotifOutput_j"), emit: motif

    script:
    """
    findMotifsGenome.pl $SUMMIT hg38 ${ID}.MotifOutput_j -size 200 -mask -p $task.cpus -len 6,8,12,16,20 -mknown /lustre/home/acct-medkkw/medlyb/wl_proj/WL234_Lib/database/JASPAR/JASPAR2020_CORE_vertebrates_non-redundant_pfms_homer.txt
    """
}

process HOMER_chipTag {
    tag "homer makeTagDirectory in $ID"
    publishDir "$params.outdir/chipTag", pattern: "${ID}_tag", mode: 'copy'

    input:
    tuple val(ID), path(BAM_t), path(BAM_c), val(Type)

    output:
    tuple val(ID), path("${ID}_tag")

    script:
    """
    makeTagDirectory ${ID}_tag/ $BAM_t  
    """
}

workflow {
    log.info """\
        ChIPalign_NF Peak Calling
        =========================
        Sample Info   : ${params.input}
        Out Dir       : ${params.outdir}
        balance_bam   : ${params.balance_bam}
        balance_pairs : ${params.balance_pairs}
        """
        .stripIndent(true)

    channel.fromPath(params.input)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def ctrl = (row.BAM_c && row.BAM_c != 'NA' && row.BAM_c != 'null' && row.BAM_c != '') ? file(row.BAM_c) : []
            tuple(row.ID, file(row.BAM_t), ctrl, row.Type)
        }
        .set { ch_bam_raw }

    def ch_for_macs2

    if ( params.balance_bam ) {
        log.info "balance_bam=true: running BamPairBalancer then MACS2"
        ch_for_macs2 = BamPairBalancer(ch_bam_raw)
            .map { ID, bt_bam, bt_bai, bt_bw, bc_bam, bc_bai, bc_bw, Type -> tuple(ID, bt_bam, bc_bam, Type) }
    } else {
        log.info "balance_bam=false: skipping BamPairBalancer, running MACS2 directly"
        ch_for_macs2 = ch_bam_raw
    }

    ch_callpeak = MACS2_callpeaks(ch_for_macs2)

    HOMER_annotatePeaks(ch_callpeak.peaks)
    // HOMER_findMotifs(ch_callpeak.summits)
    // HOMER_chipTag(ch_for_macs2)

    workflow.onComplete = {
        log.info(workflow.success ? "\nDone! See results --> $params.outdir\n" : "Oops.. someting went wrong")
    }

}
