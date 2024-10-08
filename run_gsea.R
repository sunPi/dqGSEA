# Load Functions
prepareCounts <- function(cts, verbose){
  if(verbose){
    print("Preparing RNA-seq counts for DEG and GSEA...")
  }
  
  # # Remove duplicated values
  # dup.marker <- all(duplicated(cts[,1]))
  # 
  # if(dup.marker){
  #   if(verbose){
  #     print("Found duplicates, removing...")
  #   }
  
  # Move gene names to row names
  if(verbose){
    print("Moving first column to row names...")
  }
  cts <- column_to_rownames(cts, var = "gene_symbol")
  
  # } else{
  #   if(verbose){
  #     print("Found no duplicates!")
  #   }
  # }
  # Round values just in case
  if(verbose){
    print("Rounding the integers...")
  }
  cts <- round(cts)
  
  return(cts)
}
getColData    <- function(exp.data, cts, reference){
  
  # cts <- cts[rowSums(cts == 0) == 0, ] # Remove rows with 0 values
  # if(length(exp.data[[1]]) != length(colnames(cts))){
  #   print("Detecting a difference in sample size between count data and experiment design data, conforming...")
  #   exp.data <- dplyr::filter(exp.data, exp.data[[1]] %in% colnames(cts))
  # }
  
  coldata <- exp.data[order(exp.data[[1]]), ]
  names(coldata)[2] <- "condition"
  coldata$condition <- relevel(factor(coldata$condition), ref = reference)
  
  return(coldata)
  
}
RankGenes <- function(rnaseq.counts, experiment.data, reference, verbose){
  
  # Prepare RNA-seq expression data
  cts <- prepareCounts(rnaseq.counts, verbose)
  
  # Create a 'coldata' object from the experiment design
  coldata <- getColData(exp.data = experiment.data, cts = cts, reference = reference)
  
  # Ensure the order of samples in coldata matches the counts data
  all(colnames(cts) == coldata[[1]])
  
  # Create a DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
  
  if(verbose){
    print("Running DESeq2...")
  }
  
  # Perform the DESeq analysis
  dds <- DESeq(dds)
  
  # Extract the results
  # Contrast is always built such as - "name of the column in coldata" "test" "reference"
  
  test <- levels(coldata$condition)[levels(coldata$condition) != reference]
  ctr <- c("condition", test, reference)
  
  if(verbose){
    print(paste0("Extracting list of DEG in ", test, "..."))
  }
  
  gene.list <- results(dds, tidy = T, contrast = ctr)
  names(gene.list)[1] <- "gene_symbol"
  
  # Filter out the significant genes based on p.value and order/sort them
  sig <- dplyr::filter(gene.list, gene.list$padj < 0.05)
  sig <- sig[order(sig$padj),]
  sig <- sig[order(sig$log2FoldChange, decreasing = T), ]
  sig <- sig %>% as.data.frame() %>% arrange(desc(log2FoldChange), desc(padj))
  top.10.genes <- head(sig[1:10, ])
  bot.10.genes <- tail(sig[1:10, ])
  
  if(nrow(sig) > 0){
    any.sig <- T
    print(paste0("DESeq2 analysis found ", nrow(sig), " significant genes!"))
    
  } else{
    any.sig <- F
    print("DESeq2 found 0 significant genes!")
  }
  
  p <- EnhancedVolcano(gene.list,
                       lab = gene.list[[1]],
                       x = 'log2FoldChange',
                       y = 'pvalue') + 
    theme_prism() +
    labs(title = NULL, subtitle = NULL)
  
  
  return(list(gene.list   = gene.list,
              significant = sig,
              top.10      = top.10.genes,
              bot.10      = bot.10.genes,
              design      = list(ref = reference,
                                 test = test),
              any.sig     = any.sig,
              volcano.p   = p)
  )
}
RunfGSEA <- function(gene.list, pathways, design, verbose){
  res <- gene.list %>% 
    dplyr::select(gene_symbol, log2FoldChange) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(gene_symbol) %>% 
    summarize(stat=mean(log2FoldChange))
  
  
  ranks <- deframe(res)
  head(ranks, 20)
  
  # Load the pathways into a named list
  pathways <- gmtPathways(pathways)
  # names(pathways) <- gsub("HALLMARK_", "", names(pathways))
  
  fgseaRes <- fgsea(pathways=pathways, stats=ranks)
  
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  fgseaResTidy <- fgseaResTidy[fgseaResTidy$padj < 0.1,]
  fgseaResTidy$Direction <- recode(factor(fgseaResTidy$NES > 0), "TRUE" = "Upregulated", "FALSE" = "Downregulated")
  
  custom_colors <- c("Upregulated" = "#077f97", "Downregulated" = "#400257")
  
  # title <- design$test
  
  p <- ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=Direction)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score", title = design$test) + # title=title
    theme(text = element_text(size = 8), 
          title = element_text(size = 18), 
          plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = custom_colors) +
    theme_prism()
  
  return(list(enrichments = fgseaResTidy,
              plot        = p)) 
}

# Load Libraries
pkgs <- c("DESeq2", "tidyverse", "dplyr", "fgsea", "mesocore", "ggprism", "EnhancedVolcano", "docopt")
suppressMessages(mesocore::handleRequirements(pkgs))

#---- Header ----
"Mesothelioma AI Pipeline - Differential Expression using DESeq2

Usage: run_gsea.R [options]

Options:
  -h --help                     Show this screen.
  -d --data_path=<PATH>         Path to the counts data set (RNA-seq data).
  -e --experiment_design=<PATH> Path to the table with experiment design.
  -o --outfolder=<PATH>         Folder where the results are saved.
  -g --gmt_file=<PATH>          Path to the gene signature in the '.gmt' format.
  -p --p_value=<DOUBLE>         FDR cutoff level [Default: 0.1].
  -r --reference=<STRING>       Specify the name of the reference state of the experiment.
  -s --deseq_exp=<PATH>         Specify a pre-ranked list of genes for GSEA.
  -v --verbose=<BOOLEAN>        If set to 1, prints messages verbously.
  -V --version
"-> doc

#---- Arguments ----
arguments <- docopt(doc, quoted_args = TRUE, help = TRUE)
print(arguments)

cts.path  <- normalizePath(arguments$data_path)
exp.path  <- arguments$experiment_design
outfolder <- normalizePath(arguments$outfolder, mustWork = F)
pathways  <- arguments$gmt_file
p.value   <- as.double(arguments$p_value)
reference <- arguments$reference
deseq.exp <- arguments$deseq_exp
verbose   <- as.integer(arguments$verbose)
serialise <- T

if(as.integer(arguments$verbose) == 1){
  verbose <- T
} else{
  verbose <- F
}

#----- Main -----
rnaseq.counts   <- readFile(cts.path)
experiment.data <- readFile(exp.path)

dname           <- pathways %>% basename()
outfolder       <- here(outfolder, dname)

dir.create(outfolder, F, T)

print("Ranking Genes...")
DEG             <- RankGenes(rnaseq.counts, experiment.data, reference, verbose)

print("Running fGSEA...")
fGSEA           <- RunfGSEA(DEG$gene.list, pathways, DEG$design, verbose)

if(serialise){
  print(paste("Serialising results to...", outfolder))
  writexl::write_xlsx(list(significant = DEG$significant,
                           top10       = DEG$top.10,
                           bot10       = DEG$bot.10),
                      here(outfolder,paste0(DEG$design$test, '_DEG.xlsx'))
  )
  
  png(here(outfolder,paste0(DEG$design$test, '_volcano.png')),
      width = 2400, height = 1800, res = 300)
  print(DEG$volcano.p)
  dev.off()
  
  writexl::write_xlsx(list(enriched.p  = fGSEA$enrichments),
                      here(outfolder,paste0(DEG$design$test, '_pathways.xlsx'))
  )
  
  png(here(outfolder,paste0(DEG$design$test, '_patwhays.png')),
      width = 3200, height = 1800, res = 300)
  print(fGSEA$plot)
  dev.off()
}

