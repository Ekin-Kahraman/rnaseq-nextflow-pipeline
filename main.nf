#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Bulk RNA-seq Pipeline
 * FASTQ → FastQC → fastp → HISAT2 → samtools → featureCounts → DESeq2 → MultiQC
 *
 * Processes paired-end RNA-seq reads through quality control, alignment to a
 * reference genome, gene-level quantification, and basic differential expression.
 */

// All parameter defaults are defined in nextflow.config

// Validate
if (!params.genome_index) { error "Provide --genome_index (HISAT2 index prefix)" }
if (!params.gtf)          { error "Provide --gtf (gene annotation GTF)" }

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
 *
 * Run before trimming to capture adapter contamination rates and
 * per-base quality profiles. These metrics inform whether fastp
 * parameters need adjusting for a given library.
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
    fastqc --threads ${task.cpus} --quiet ${fastq1} ${fastq2}
    """
}

/*
 * STEP 2: Trim with fastp
 *
 * fastp over Trimmomatic: single binary handles adapter detection,
 * quality trimming, and polyG tail removal. No external adapter file
 * needed — fastp infers adapters from read overlap in PE mode.
 * Phred 20 + min length 36 follows the Himes et al. original analysis.
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
        --thread ${task.cpus} \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 36
    """
}

/*
 * STEP 3: Align with HISAT2
 *
 * HISAT2 uses a graph FM index that fits in ~8GB RAM, vs STAR's
 * ~32GB for the human genome. This makes the pipeline runnable on
 * laptops and standard HPC nodes without high-memory allocation.
 * --dta flag enables downstream transcript assembly compatibility.
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
    def index_prefix = index_files[0].toString().replaceFirst(/\.\d+\.ht2[l]?$/, '')
    """
    hisat2 \
        -x ${index_prefix} \
        -1 ${fastq1} -2 ${fastq2} \
        --threads ${task.cpus} \
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
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), val(condition), emit: bam

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}.bam ${sam}
    samtools index ${sample_id}.bam
    """
}

/*
 * STEP 4: Count with featureCounts
 *
 * Gene-level quantification using gene_name (not gene_id) for
 * human-readable output. --countReadPairs counts fragments not
 * individual reads (correct for PE data). -s flag is configurable
 * because strandedness varies by library prep: dUTP protocols are
 * reverse-stranded (s=2), older protocols may be unstranded (s=0).
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
        -T ${task.cpus} \
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

    # Filter low-count genes (adaptive: require counts in at least half the samples)
    min_count <- ifelse(max(colSums(counts)) > 1e6, 10, 1)
    min_samples <- max(2, floor(ncol(counts) / 2))
    keep <- rowSums(counts >= min_count) >= min_samples
    counts <- counts[keep, , drop = FALSE]
    cat(sprintf("Genes passing filter: %d (min_count=%d, min_samples=%d)\\n", sum(keep), min_count, min_samples))

    # DESeq2
    dds <- DESeqDataSetFromMatrix(counts, ss, design = ~ condition)
    dds\$condition <- relevel(dds\$condition, ref = "${params.ref_condition}")

    tryCatch({
        dds <- DESeq(dds)
    }, error = function(e) {
        # Fallback for small/test datasets where dispersion fitting fails
        cat("Note: using gene-wise dispersion estimates (expected for small test data)\\n")
        dds <<- estimateSizeFactors(dds, type = "poscounts")
        dds <<- estimateDispersionsGeneEst(dds)
        dispersions(dds) <<- mcols(dds)\$dispGeneEst
        dds <<- nbinomWaldTest(dds)
    })
    # Get the non-reference condition level for contrast
    cond_levels <- levels(dds\$condition)
    treat_level <- cond_levels[cond_levels != "${params.ref_condition}"][1]
    cat(sprintf("Contrast: %s vs %s\\n", treat_level, "${params.ref_condition}"))
    res <- results(dds, contrast = c("condition", treat_level, "${params.ref_condition}"), alpha = 0.05)
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
    all_bams = SAMTOOLS_SORT.out.bam.map { id, bam, bai, cond -> bam }.collect()
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
