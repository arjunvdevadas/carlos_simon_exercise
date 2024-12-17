# author Arjun

# libraries to load
{
  library(GenomicAlignments)
  library(rtracklayer)
  library(GenomicRanges)
  library(Rsamtools)
  library(BiocParallel)
  library(DESeq2)
  library(EnhancedVolcano)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  library(ggplot2)
  library(ggpubr)
  library(tidyr)
  library(dplyr)
  library(doParallel)
  library(parallel)
}

# datasets and variables
{
  # fastq file details
  fastq_dir = "/usr/users/adevada/Projects/cs/rawdata/fastq/" # path to location of extracted fastq files after fastq-dump
  samples = c("SRR5223500", "SRR5223522", "SRR5223543", "SRR5223505", "SRR5223547", "SRR5223570") # sample names
  r1_r2_sep = c("_1", "_2")
  extension = "fastq"
  
  # fastq file names
  fastq_names = lapply(samples, function(id) { list(paste0(fastq_dir, id, r1_r2_sep[1], ".", extension), paste0(fastq_dir, id, r1_r2_sep[2], ".", extension)) })
  
  # mapping index
  hg38_ercc_index = "/usr/users/adevada/Projects/cs/rawdata/STARindices/" # combined genome index for hg38 and ercc spike in genome
  star = "/usr/users/adevada/bin/STAR_2.7.11b/Linux_x86_64_static/STAR"
  
  # annotation
  anno = makeGRangesFromDataFrame(get(load("/usr/users/adevada/Projects/cs/rawdata/genes_spikeins.RData")), keep.extra.columns = T) # canonical gene annotation
  
  # output directories
  fastqc_out_dir = "/usr/users/adevada/Projects/cs/output/fastqc_raw/"
  cutadapt_out_dir = "/usr/users/adevada/Projects/cs/output/cutadapt_raw/"
  star_out_dir = "/usr/users/adevada/Projects/cs/output/star_qc/"
  
  # bam files (STAR mapping)
  bam_files = paste0(star_out_dir, samples, "_1.fastq_Aligned.sortedByCoord.out.bam") # to account for STAR output format
  
  dir.create(fastqc_out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(cutadapt_out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(star_out_dir, showWarnings = FALSE, recursive = TRUE)
}

# functions
{
  # run fastqc for samples
  run_fastqc = function(files) {
    log1 = file.path(fastqc_out_dir, paste0(basename(files[[1]]), ".log"))
    log2 = file.path(fastqc_out_dir, paste0(basename(files[[2]]), ".log"))
    cmd1 = paste("fastqc", files[[1]], "-o", fastqc_out_dir, "--nogroup", ">", log1, "2>&1")
    cmd2 = paste("fastqc", files[[2]], "-o", fastqc_out_dir, "--nogroup", ">", log2, "2>&1")
    cat("Running command for file 1:", cmd1, "\n")
    system(cmd1)
    cat("Running command for file 2:", cmd2, "\n")
    system(cmd2)
  }
  
  # run cutadapt for samples
  run_cutadapt = function(files) {
    log = file.path(cutadapt_out_dir, paste0(basename(files[[1]]), "_cutadapt.log"))
    output_file1 = file.path(cutadapt_out_dir, basename(files[[1]]))
    output_file2 = file.path(cutadapt_out_dir, basename(files[[2]]))
    if(forward_adapter != "" & reverse_adapter != ""){
      cmd = paste("cutadapt -q 30 -m 30 -a ", forward_adapter," -A ", reverse_adapter, " -o", shQuote(output_file1), "-p", shQuote(output_file2), shQuote(files[[1]]), shQuote(files[[2]]), ">", shQuote(log), "2>&1")
    } else {
      cmd = paste("cutadapt -q 30 -m 30 ", forward_adapter,"", reverse_adapter, " -o", shQuote(output_file1), "-p", shQuote(output_file2), shQuote(files[[1]]), shQuote(files[[2]]), ">", shQuote(log), "2>&1")
    }
    cat("Running command:", cmd, "\n")
    system(cmd)
  }
  
  # mapping function (parameters are set to output only uniquely mapped alignments)
  run_star = function(files, ext) {
    input_file1 = file.path(cutadapt_out_dir, basename(files[[1]]))
    input_file2 = file.path(cutadapt_out_dir, basename(files[[2]]))
    output_prefix = file.path(star_out_dir, paste0(basename(files[[1]]), "_"))
    log = file.path(star_out_dir, paste0(basename(files[[1]]), "_star.log"))
    
    if(ext == "fastq") {
      cmd = paste("/usr/users/adevada/bin/STAR_2.7.11b/Linux_x86_64_static/STAR", 
                  "--runThreadN 1",
                  "--genomeDir", shQuote(star_index_file_path),
                  "--readFilesIn", shQuote(input_file1), shQuote(input_file2),
                  "--outFileNamePrefix", shQuote(output_prefix),
                  "--outSAMmultNmax 1", 
                  "--outSAMtype BAM SortedByCoordinate",
                  ">", shQuote(log), "2>&1")
    } else {
      cmd = paste("/usr/users/adevada/bin/STAR_2.7.11b/Linux_x86_64_static/STAR", 
                  "--runThreadN 1",
                  "--readFilesCommand zcat",
                  "--genomeDir", shQuote(star_index_file_path),
                  "--readFilesIn", shQuote(input_file1), shQuote(input_file2),
                  "--outFileNamePrefix", shQuote(output_prefix),
                  "--outSAMmultNmax 1", 
                  "--outSAMtype BAM SortedByCoordinate",
                  ">", shQuote(log), "2>&1")
    }
    cat("Running command:", cmd, "\n")
    system(cmd)
  }
  
  # indexing bam files
  index_bam = function(files) {
    bam_file = file.path(star_out_dir, paste0(basename(files[[1]]), "_Aligned.sortedByCoord.out.bam"))
    index_log = file.path(star_out_dir, paste0(basename(files[[1]]), "_index.log"))
    
    cat("BAM file:", bam_file, "\n")
    cat("Index log file:", index_log, "\n")
    
    cmd = paste("samtools index", shQuote(bam_file), ">", shQuote(index_log), "2>&1")
    
    cat("Running command:", cmd, "\n")
    
    result = try(system(cmd, intern = TRUE), silent = TRUE)
    
    if (inherits(result, "try-error")) {
      cat("Error encountered while indexing BAM. Check the log:", index_log, "\n")
      return(result)
    }
    
    return(result)
  }
  
  # function to generate counts dataframe
  count_reads_strict = function(bam_files, granges_annotation) {
    
    # load BAM files into a BamFileList
    bam_file_list = BamFileList(bam_files, yieldSize = 1000000)
    
    # use summarizeOverlaps to count reads overlapping the GRanges annotation
    counts = summarizeOverlaps(
      features = granges_annotation,
      reads = bam_file_list,
      mode = "IntersectionStrict", # strict counting mode
      singleEnd = FALSE, 
      ignore.strand = FALSE,  
      fragments = TRUE  
    )
    
    # extract counts matrix 
    counts_matrix = assay(counts)
    annotation_df = as.data.frame(granges_annotation)
    counts_df = as.data.frame(counts_matrix)
    
    final_df = cbind(annotation_df, counts_df)
    
    return(final_df)
  }
}

# preprocessing
{
  # fastqc for raw data
  fastqc_results = mclapply(fastq_names, run_fastqc, mc.cores = detectCores())
  
  # cutadapt QC for raw data
  forward_adapter = ""
  reverse_adapter = ""
  cutadapt_results = mclapply(fastq_names, run_cutadapt, mc.cores = detectCores())
  
  # fastqc after cutadapt
  fastqc_out_dir = "/usr/users/adevada/Projects/cs/output/fastqc_QC/"
  dir.create(fastqc_out_dir, showWarnings = FALSE, recursive = TRUE)
  fastq_names = lapply(samples, function(id) { list(paste0(cutadapt_out_dir, id, r1_r2_sep[1], ".", extension), paste0(cutadapt_out_dir, id, r1_r2_sep[2], ".", extension)) })
  fastqc_results = mclapply(fastq_names, run_fastqc, mc.cores = detectCores())
  
  # STAR mapping
  star_index_file_path = hg38_ercc_index
  dir.create(star_out_dir, showWarnings = FALSE, recursive = TRUE)
  for (files in fastq_names) { run_star(files, ext = extension) }
  
  # indexing bam files
  for (files in fastq_names) { index_bam(files) }
  
  # counts generation
  outfile = "/usr/users/adevada/Projects/cs/output/genes_spikeins_counts.RData"
  final_counts = count_reads_strict(bam_files, anno)
  save(final_counts, file = outfile)
}

# analysis
{
  data = get(load("arjunserver/Projects/cs/output/genes_spikeins_counts.RData"))
  data[which(data$type == "spikein"), "gene_biotype"] = "spikein"
  
  # distribution of counts 
  {
    columns_of_interest = basename(bam_files)
    long_counts = pivot_longer(data, 
                               cols = all_of(columns_of_interest), 
                               names_to = "sample", 
                               values_to = "count")  
    long_counts = long_counts[-which(long_counts$count == 0), ] # removing all zero counts
    long_counts$log_count = log(long_counts$count + 1)
    
    # plot the log counts distribution for all samples
    ggplot(long_counts, aes(x = log_count, color = sample)) +
      geom_density(size = 1) +
      theme_classic() +
      labs(title = "Density distribution of counts",
           x = "Log2(count + 1)",
           y = "Density") +
      theme(
        legend.title = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(color = "black", size = 14),
        axis.title.y = element_text(color = "black", size = 14),
        plot.title = element_text(color = "black", size = 16, hjust = 0.5),
      )
  }
  
  # distribution of counts in different RNA biotypes
  {
    biotypes_of_interest = c("protein_coding", "lncRNA", "tRNA", "rRNA", "miRNA", "pseudogene", "spikein")
    
    filtered_counts = long_counts %>%
      filter(gene_biotype %in% biotypes_of_interest)
    
    biotype_summary = filtered_counts %>%
      group_by(sample, gene_biotype) %>%
      summarize(total_count = sum(count), .groups = "drop")
    
    # barplot showing total counts for each biotype per sample
    ggplot(biotype_summary, aes(x = sample, y = total_count, fill = gene_biotype)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      theme_minimal() +
      labs(title = "Total Counts of Gene Biotypes Across Samples",
           x = "Sample",
           y = "Total Count",
           fill = "Gene Biotype") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
      )
  }
  
  # spike in counts and normalization factors
  {
    spikein = colSums(data[which(data$type == "spikein"), columns_of_interest])
    
    # barplot showing spike in read numebrs
    df = data.frame(index = seq_along(spikein), value = spikein) 
    ggplot(df, aes(x = factor(index), y = spikein)) +
      geom_bar(stat = "identity", fill = "skyblue") +
      theme_minimal() +
      labs(title = "Simple Barplot",
           x = "Index",
           y = "Value")
    
    spikein_norm = spikein/spikein[1]
    spikein_norm = sf
    norm_data = data
    # normalized counts
    norm_data[, columns_of_interest[1]] = norm_data[, columns_of_interest[1]]/spikein_norm[1]
    norm_data[, columns_of_interest[2]] = norm_data[, columns_of_interest[2]]/spikein_norm[2]
    norm_data[, columns_of_interest[3]] = norm_data[, columns_of_interest[3]]/spikein_norm[3]
    norm_data[, columns_of_interest[4]] = norm_data[, columns_of_interest[4]]/spikein_norm[4]
    norm_data[, columns_of_interest[5]] = norm_data[, columns_of_interest[5]]/spikein_norm[5]
    norm_data[, columns_of_interest[6]] = norm_data[, columns_of_interest[6]]/spikein_norm[6]
  }
  
  # normalized FOXP3 counts
  {
    ctrl = columns_of_interest[1:3]
    trt = columns_of_interest[4:6]
    
    foxp3_counts = norm_data[which(norm_data$gene == "FOXP3"), ]
    foxp3_ctrl = as.numeric(foxp3_counts[, ctrl])
    foxp3_trt = as.numeric(foxp3_counts[, trt])
    
    # Create a data frame with counts for both groups
    combined_df = data.frame(
      sample = c(ctrl, trt),
      count = c(foxp3_ctrl, foxp3_trt),
      group = rep(c("Control", "Treatment"), times = c(length(ctrl), length(trt)))
    )
    
    ggplot(combined_df, aes(x = group, y = count, fill = group)) +
      geom_boxplot() +
      theme_minimal() +
      labs(title = "FOXp3 Gene Counts in Control and Treatment Groups",
           x = "Group",
           y = "FOXp3 Count") +
      theme(legend.position = "none") +
      stat_compare_means(method = "t.test", label = "p.format")
  }
  
  # differential expression and GO analysis
  {
    counts = data[-which(data$type == "spikein"), columns_of_interest]
    rownames(counts) = data[-which(data$type == "spikein"), ]$gene_id
    counts = counts[-which(rowSums(counts) == 0), ]
    count_data = counts
    
    sample_info = data.frame(
      row.names = colnames(count_data),
      condition = factor(rep(c("Control", "Treated"), each = 3)),
      donor = factor(rep(1:3, times = 2))
    )
    
    dds = DESeqDataSetFromMatrix(countData = count_data,
                                 colData = sample_info,
                                 design = ~ donor + condition)
    
    sizeFactors(dds) = spikein_norm
    dds = DESeq(dds)
    results_deseq2 = results(dds)
    filtered_results = results_deseq2[!is.na(results_deseq2$padj), ]
    significant_results = filtered_results[filtered_results$padj < 0.05, ]
    
    res_ordered = results_deseq2[order(results_deseq2$padj), ]
    
    # MA plot
    plotMA(results_deseq2, ylim = c(-5, 5), main = "MA Plot")
    
    top20_up = significant_results[order(significant_results$log2FoldChange, decreasing = TRUE), ][1:20, ]
    top20_down = significant_results[order(significant_results$log2FoldChange), ][1:20, ]
    top_genes = c(rownames(top20_up), rownames(top20_down))
    
    # volcano plot
    lab = rownames(count_data)
    EnhancedVolcano(results_deseq2,
                    lab = lab,
                    selectLab = "FOXP3",
                    labSize = 3,
                    x = 'log2FoldChange',
                    y = 'pvalue',
                    title = 'Differential Expression',
                    subtitle = "",
                    pCutoff = 0.05,
                    FCcutoff = 1)
    
    EnhancedVolcano(results_deseq2,
                    lab = lab,
                    selectLab = top_genes,
                    labSize = 3,
                    x = 'log2FoldChange',
                    y = 'pvalue',
                    title = 'Differential Expression',
                    subtitle = "",
                    pCutoff = 0.05,
                    FCcutoff = 1)
    
    # extract upregulated and downregulated genes
    upregulated_genes = rownames(res_ordered)[res_ordered$log2FoldChange > 1 & res_ordered$padj < 0.05]
    upregulated_genes = upregulated_genes[which(!is.na(upregulated_genes))]
    which(upregulated_genes == "FOXP3")
    downregulated_genes = rownames(res_ordered)[res_ordered$log2FoldChange < -1 & res_ordered$padj < 0.05]
    downregulated_genes = downregulated_genes[which(!is.na(downregulated_genes))]
    
    
    # perform GO analysis
    upregulated_go = enrichGO(gene = upregulated_genes,
                              OrgDb = org.Hs.eg.db,
                              keyType = "SYMBOL",
                              ont = "BP",
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.2)
    
    downregulated_go = enrichGO(gene = downregulated_genes,
                                OrgDb = org.Hs.eg.db,
                                keyType = "SYMBOL",
                                ont = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.2)
    
    # plot GO results for upregulated genes
    barplot(upregulated_go, showCategory = 15, title = "GO Enrichment for Upregulated Genes")
    
    # plot GO results for downregulated genes
    barplot(downregulated_go, showCategory = 15, title = "GO Enrichment for Downregulated Genes")
    
    result = as.data.frame(results_deseq2)
    # optionally, write results to CSV
    write.csv(as.data.frame(res_ordered), file = "deseq2_results.csv")
  }
}

# other plots
{
  {
    # mapping quality
    bam_file = "arjunserver/Projects/cs/output/star_qc/SRR5223570_1.fastq_Aligned.sortedByCoord.out.bam"
    bam_data = scanBam(bam_file, param = ScanBamParam(what = "mapq"))
    mapq_scores = bam_data[[1]]$mapq
    mapq_df = data.frame(mapping_quality = mapq_scores)
    ggplot(mapq_df, aes(x = mapping_quality)) +
      geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
      theme_minimal() +
      labs(title = "Mapping quality : D3_trt",
           x = "Mapping Quality Score",
           y = "Frequency") +
      theme(axis.text = element_text(color = "black"),
            axis.title = element_text(color = "black"),
            plot.title = element_text(color = "black", size = 16, hjust = 0.5))
  }
  
  # some stats, hard coded
  samples = c("SRR5223500", "SRR5223522", "SRR5223543", "SRR5223505", "SRR5223547", "SRR5223570")
  samples.alias = c("D1_ctrl", "D2_ctrl", "D3_ctrl", "D1_trt", "D2_trt", "D3_trt") 
  sample.colors = alpha(c("grey", "grey", "grey", "orange", "orange", "orange"), 0.6)
  raw.read.numbers = c(24981181, 17279017, 19040570, 20979328, 18503753, 4423985) # number of reads in original fastq
  removed.by.cutadapt = c(85.6, 82.3, 80.9, 87.1, 85.2, 81.6) # percent of reads remaining after cutadapt QC
  mapping.percent = c(91.29, 89.37, 89.63, 88.65, 91.03, 91.89) # percent of uniquely mapped reads
  spikein.percent = 100 * spikein/c(19520387, 12712259, 13814437, 16193012, 14346425, 3315672) # percent reads aligning to ercc spike in genome
  
  data = data.frame(
    sample = samples.alias,
    value = spikein.percent, 
    color = sample.colors
  )
  
  # Plot the barplot using ggplot2
  ggplot(data, aes(x = sample, y = value, fill = sample)) +
    geom_bar(stat = "identity", aes(fill = color), color = "black") +
    scale_fill_identity() +
    theme_minimal() +
    labs(title = "Spike-in percent",
         x = "Sample",
         y = "Percent spike-in reads") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 12),
          axis.title.x = element_text(color = "black", size = 14),
          axis.title.y = element_text(color = "black", size = 14),
          plot.title = element_text(color = "black", size = 16, hjust = 0.5),
          legend.text = element_text(color = "black"),
          legend.title = element_text(color = "black"))
  
}
