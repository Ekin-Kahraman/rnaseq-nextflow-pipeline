#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * RNA-seq Pipeline
 * FastQC → fastp → HISAT2 → featureCounts → DESeq2 → MultiQC
 */

params.samplesheet = "${projectDir}/assets/samplesheet.csv"
params.genome_index = null
params.gtf = null
params.outdir = "results"
params.strandedness = 2  // 0=unstranded, 1=forward, 2=reverse

// Validate required params
if (!params.genome_index) { error "Provide --genome_index (path to HISAT2 index prefix)" }
if (!params.gtf) { error "Provide --gtf (path to gene annotation GTF)" }

/*
 * Read samplesheet: sample_id, fastq_1, fastq_2, condition
 */
Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row -> tuple(row.sample_id, file(row.fastq_1), file(row.fastq_2), row.condition) }
    .set { reads_ch }

/*
 * STEP 1: FastQC on raw reads
 */
process FASTQC_RAW {
    tag "$sample_id"
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir "${params.outdir}/fastqc_raw", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2), val(condition)

    output:
    path("*.html"), emit: html
    path("*.zip"), emit: zip

    script:
    """
    fastqc --threads 2 --quiet ${fastq1} ${fastq2}
    """
}

/*
 * STEP 2: Trim with fastp
 */
process FASTP {
    tag "$sample_id"
    container 'quay.io/biocontainers/fastp:0.23.4--hadf994f_3'
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: '*.json'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2), val(condition)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq.gz"), path("${sample_id}_trimmed_R2.fastq.gz"), val(condition), emit: trimmed
    path("${sample_id}_fastp.json"), emit: json

    script:
    """
    fastp \
        -i ${fastq1} -I ${fastq2} \
        -o ${sample_id}_trimmed_R1.fastq.gz -O ${sample_id}_trimmed_R2.fastq.gz \
        --json ${sample_id}_fastp.json \
        --thread 4 \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 36
    """
}

/*
 * STEP 3: Align with HISAT2
 */
process HISAT2_ALIGN {
    tag "$sample_id"
    container 'quay.io/biocontainers/hisat2:2.2.1--hdbdd923_6'
    publishDir "${params.outdir}/hisat2", mode: 'copy', pattern: '*.log'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2), val(condition)
    path(index_files)

    output:
    tuple val(sample_id), path("${sample_id}.sam"), val(condition), emit: sam
    path("${sample_id}_hisat2.log"), emit: log

    script:
    def index_prefix = index_files[0].toString().replaceAll(/\.\d+\.ht2$/, '')
    """
    hisat2 \
        -x ${index_prefix} \
        -1 ${fastq1} -2 ${fastq2} \
        --threads 4 \
        --summary-file ${sample_id}_hisat2.log \
        --new-summary \
        --dta \
        -S ${sample_id}.sam
    """
}

/*
 * STEP 3b: Sort and index BAM
 */
process SAMTOOLS_SORT {
    tag "$sample_id"
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(sam), val(condition)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), val(condition), emit: bam

    script:
    """
    samtools sort -@ 2 -o ${sample_id}.bam ${sam}
    samtools index ${sample_id}.bam
    """
}

/*
 * STEP 4: Count with featureCounts
 */
process FEATURECOUNTS {
    container 'quay.io/biocontainers/subread:2.0.6--he4a0461_2'
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    path(bams)
    path(gtf)

    output:
    path("gene_counts.txt"), emit: counts
    path("gene_counts.txt.summary"), emit: summary

    script:
    """
    featureCounts \
        -a ${gtf} \
        -o gene_counts.txt \
        -T 4 \
        -p --countReadPairs \
        -s ${params.strandedness} \
        -t exon \
        -g gene_name \
        ${bams}
    """
}

/*
 * STEP 5: Basic DESeq2 differential expression
 */
process DESEQ2 {
    container 'quay.io/biocontainers/bioconductor-deseq2:1.42.0--r43hf17093f_2'
    publishDir "${params.outdir}/deseq2", mode: 'copy'

    input:
    path(counts)
    path(samplesheet)

    output:
    path("deseq2_results.csv"), emit: results
    path("volcano_plot.png"), emit: volcano
    path("pca_plot.png"), emit: pca

    script:
    """
    #!/usr/bin/env Rscript
    library(DESeq2)
    library(ggplot2)

    # Read counts (skip featureCounts header, drop annotation columns)
    raw <- read.delim("${counts}", skip = 1, check.names = FALSE)
    genes <- raw[, 1]
    counts <- as.matrix(raw[, 7:ncol(raw)])
    rownames(counts) <- genes

    # Clean sample names (strip .bam suffix)
    colnames(counts) <- sub("\\\\.bam\$", "", basename(colnames(counts)))

    # Read samplesheet for condition labels
    ss <- read.csv("${samplesheet}")
    ss <- ss[match(colnames(counts), ss\$sample_id), ]
    stopifnot(all(colnames(counts) == ss\$sample_id))

    # DESeq2
    dds <- DESeqDataSetFromMatrix(counts, ss, design = ~ condition)
    dds\$condition <- relevel(dds\$condition, ref = "negative")

    # Filter low-count genes (keep genes with at least 1 count in 2+ samples)
    keep <- rowSums(counts >= 1) >= 2
    counts <- counts[keep, , drop = FALSE]
    cat(sprintf("Genes passing filter: %d / %d\\n", sum(keep), length(keep)))

    dds <- DESeqDataSetFromMatrix(counts, ss, design = ~ condition)
    dds\$condition <- relevel(dds\$condition, ref = "negative")

    tryCatch({
        dds <- DESeq(dds, sfType = "poscounts")
    }, error = function(e) {
        # Fallback for small/test datasets where dispersion fitting fails
        cat("Note: using gene-wise dispersion estimates (expected for small test data)\\n")
        dds <<- estimateSizeFactors(dds, type = "poscounts")
        dds <<- estimateDispersionsGeneEst(dds)
        dispersions(dds) <<- mcols(dds)\$dispGeneEst
        dds <<- nbinomWaldTest(dds)
    })
    res <- results(dds, contrast = c("condition", "positive", "negative"), alpha = 0.05)
    res <- res[order(res\$padj), ]

    # Save results
    res_df <- as.data.frame(res)
    res_df\$gene <- rownames(res_df)
    res_df <- res_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
    write.csv(res_df, "deseq2_results.csv", row.names = FALSE)

    n_sig <- sum(res_df\$padj < 0.05, na.rm = TRUE)
    cat(sprintf("Significant genes (FDR < 0.05): %d\\n", n_sig))

    # Volcano plot
    res_df\$sig <- ifelse(!is.na(res_df\$padj) & res_df\$padj < 0.05, "Significant", "NS")
    png("volcano_plot.png", width = 800, height = 600, res = 150)
    print(ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), colour = sig)) +
        geom_point(alpha = 0.5, size = 1) +
        scale_colour_manual(values = c("NS" = "grey70", "Significant" = "#E53935")) +
        theme_minimal() +
        labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10(p-value)"))
    dev.off()

    # PCA
    vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
    pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
    pct_var <- round(100 * attr(pca_data, "percentVar"))
    png("pca_plot.png", width = 800, height = 600, res = 150)
    print(ggplot(pca_data, aes(x = PC1, y = PC2, colour = condition)) +
        geom_point(size = 4) +
        theme_minimal() +
        labs(title = "PCA", x = paste0("PC1 (", pct_var[1], "%)"), y = paste0("PC2 (", pct_var[2], "%)")))
    dev.off()
    """
}

/*
 * STEP 6: MultiQC
 */
process MULTIQC {
    container 'quay.io/biocontainers/multiqc:1.27.1--pyhdfd78af_0'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path(fastqc_zips)
    path(fastp_jsons)
    path(hisat2_logs)
    path(fc_summary)

    output:
    path("multiqc_report.html"), emit: report

    script:
    """
    multiqc . --filename multiqc_report
    """
}

/*
 * WORKFLOW
 */
workflow {
    // Raw QC
    FASTQC_RAW(reads_ch)

    // Trim
    FASTP(reads_ch)

    // Align
    index_ch = Channel.fromPath("${params.genome_index}*.ht2").collect()
    HISAT2_ALIGN(FASTP.out.trimmed, index_ch)

    // Sort
    SAMTOOLS_SORT(HISAT2_ALIGN.out.sam)

    // Count — collect all BAMs
    all_bams = SAMTOOLS_SORT.out.bam.map { id, bam, cond -> bam }.collect()
    gtf_ch = Channel.fromPath(params.gtf)
    FEATURECOUNTS(all_bams, gtf_ch)

    // DESeq2
    samplesheet_ch = Channel.fromPath(params.samplesheet)
    DESEQ2(FEATURECOUNTS.out.counts, samplesheet_ch)

    // MultiQC
    MULTIQC(
        FASTQC_RAW.out.zip.collect(),
        FASTP.out.json.collect(),
        HISAT2_ALIGN.out.log.collect(),
        FEATURECOUNTS.out.summary
    )
}
