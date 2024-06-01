# DESeq2 Functions
removeDuplicates <- function(x){
  x <- x[!duplicated(rownames(x)),]
  return(x)
}
convertGeneType <- function(gene.list, to.type = "symbol"){
  #https://bioinformatics.stackexchange.com/questions/5229/converting-gene-symbol-to-ensembl-id-in-r
  library(EnsDb.Hsapiens.v86)
  # 1. Convert from ensembl.gene to gene.symbol
  if(to.type == "symbol"){
    converted <- ensembldb::select(EnsDb.Hsapiens.v86, keys= gene.list, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
    return(converted)
  }
  if(to.type == "ensemble"){
    # 2. Convert from gene.symbol to ensembl.gene
    converted <- ensembldb::select(EnsDb.Hsapiens.v86, keys= gene.list, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
    return(converted)
  } else{
    stop("Specified the wrong encode destination type, please use either 'symbol' or 'ensemble'!")
  }
}
checkGeneEncoding <- function(gene.id, verbose){
  # Function to test if a string matches the gene name pattern
  isHugo <- function(name) {
    # gene.pattern <- "^[A-Z0-9-]+$"
    gene_pattern <- "^[A-Z]+[0-9]+([-A-Z0-9]*)?$"
    return(grepl(gene_pattern, name))
  }
  
  if (any(startsWith(gene.id, "ENSG"))) {
    if(verbose){
      print("Converting from ENSEMBL id's to HUGO symbols...")
    }
    gene.id <- convertGeneType(gene.id, to.type = "symbol")
    return(gene.id)
    
  } else{
    # stop("Your gene id's are not encoded in either HUGO or ENSEMBLE type. This function supports 
    #      the conversion between these two only!")
    warning("Your gene id's are not encoded in either HUGO or ENSEMBLE type. This function supports 
         the conversion between these two only! THIS FUNCTION IS NOT YET COMPLETE, MIGHT BE OK...CHECK...")
    return(list("SYMBOL" = gene.id))
  }
}
getColData <- function(exp.data, cts){
  if(length(exp.data[[1]]) != length(colnames(cts))){
    print("Detecting a difference in sample size between count data and experiment design data, conforming...")
    exp.data <- dplyr::filter(exp.data, exp.data[[1]] %in% colnames(cts))
  }
  
  coldata <- data.frame(row.names = colnames(cts),
                        "condition" = factor(exp.data[[2]]))
  return(coldata)
}
getExpPerms <- function(coldata){
  condition <- as.vector(unique(coldata$condition))
  
  all_perms <- combinat::permn(condition)
  names(all_perms) <- lapply(all_perms, paste0, collapse = "_vs_")
  
  return(all_perms)
}
getDESeqRes <- function(dds, contrast, alpha.value, p.value){
  
  # Get genes that are differentially expressed
  if(!is.null(contrast)){
    res <- results(dds, alpha = alpha.value, contrast = contrast)
  } else{
    res <- results(dds, alpha = alpha.value)
  }
  
  res$symbol <- row.names(res)
  
  # Filter for significance
  sigs <- na.omit(res)
  sigs <- sigs[sigs$padj < p.value,]
  
  sigs$symbol <- rownames(sigs)
  
  # MA plot
  p.ma.data <- DESeq2::plotMA(res, returnData = T)
  
  # Volcano plot
  p.volc <- EnhancedVolcano::EnhancedVolcano(res, x = "log2FoldChange", y = "padj", lab = res$symbol)
  
  if(nrow(sigs) > 0){
    return(list(significant = sigs,
                res         = as.data.frame(res),
                vplot       = p.volc,
                p.ma.data   = p.ma.data,
                any.sig     = T))
    
  } else{
    return(list(significant = sigs,
                res         = as.data.frame(res),
                vplot       = p.volc,
                p.ma.data   = p.ma.data,
                any.sig     = F))
  }
}
runDESeq <- function(cts, coldata){
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ condition)
  
  
  dds <- dds %>% DESeq()
  
}
getRankedGeneList <- function(data.path, exp.path, outfolder,
                              alpha.value = 0.1, p.value = 0.05, 
                              do.contrasts = T, 
                              verbose = T, serialise = F){
  cts            <- readFile(data.path)
  exp.data       <- readFile(exp.path)
  
  if(any(is.na(cts))){
    if(verbose){
      print("Detected NA values, removing...")
    }
    cts <- na.omit(cts)
  }
  
  # Remove duplicated values
  cts <- cts[!duplicated(cts$gene_name),]
  
  # Assign gene_id as rownames
  row.names(cts) <- cts$gene_name
  cts <- cts[,-1]
  
  # Create a 'coldata' object from the experiment design
  coldata <- getColData(exp.data, cts)
  
  # Create a DESeq2 matrix, filter the matrix for low counts (=cts.filter)
  # and run the DESeq() algorithm
  dds     <- runDESeq(cts, coldata)
  
  # Calculate all possible permutations of the experiment 
  # (e.g. a vs. b ~ b vs. a)
  all_perms <- getExpPerms(coldata)
  
  contrasts <- list()
  res <- list()
  dname <- basename(dirname(exp.path))
  
  
  # Extract the results for all possible contrasts of the experiments
  if(do.contrasts){
    for(i in seq_along(all_perms)){
      pname <- names(all_perms)[i]
      dir.create(here::here(outfolder, dname, pname), showWarnings = F, recursive = T)
      
      print(paste0("Looking for DEGs in ", pname))
      contrasts[[pname]] <- c(names(coldata), paste0(all_perms[[pname]]))
      res[[pname]]       <- getDESeqRes(dds,
                                        contrast = contrasts[[pname]], 
                                        alpha.value = alpha.value,
                                        p.value = p.value)
      if(!res[[pname]]$any.sig){
        print(paste0("After correcting for multiple testing, there were no significant genes left for contrast ", pname))
      } else{
        if(verbose){
          print(paste0("Found ",nrow(res[[pname]]$significant), " significant gene(s) in ", pname))
        }
        png(here::here(outfolder, dname, pname, paste0('MAplot.png')))
        plotMA(res[[pname]]$p.ma.data, ylim=c(-2,2))
        dev.off()
        
        png(here::here(outfolder, dname, pname, paste0('Vplot.png')))
        print(res[[pname]]$vplot)
        dev.off()
        
      }
    }
  } else{
    res <- getDESeqRes(dds, 
                       alpha.value = alpha.value,
                       median.filter = median.filter)
  }
  
  
  dir.create(here::here(outfolder, dname, pname), showWarnings = F, recursive = T)
  # writeLines(dname, here::here(outfolder, ".config.txt"))
  if(serialise){
    saveRDS(res, here::here(outfolder, dname, "diffseq_exp.RDS"))
  }
  
  return(res)
}
# fGSEA Functions
getEnrichedPathways <- function(experiment, pathways, nperm, p.value, dname, 
                                outfolder, verbose){
  runfGSEA <- function(gene.list, pathways, nperm, p.value = 0.05,
                       title = NULL, outfolder = NULL, verbose = F)
  {
    res <- gene.list %>% 
      dplyr::select(symbol, stat) %>% 
      na.omit() %>% 
      distinct() %>% 
      group_by(symbol) %>% 
      summarize(stat=mean(stat))
    
    ranks <- deframe(res)
    head(ranks, 20)
    
    # Load the pathways into a named list
    enriched.p <- gmtPathways(pathways)
    
    # Look at them all if you want (uncomment)
    # pathways.hallmark
    
    if(verbose){
      # Show the first few pathways, and within those, show only the first few genes. 
      print("First few significant pathways and their genes repectively:")
      enriched.p %>% 
        head() %>% 
        lapply(head) %>% 
        print()
    }
    
    if(!is.null(nperm)){
      # Run the fGSEA with n permutations
      print(paste0("Running fGSEA with ", nperm, " permutations..."))
      fgseaRes <- fgsea(pathways=enriched.p, stats=ranks, nPermSimple=nperm)
    } else {
      # Run the fGSEA without permutations for a faster result
      fgseaRes <- fgsea(pathways=enriched.p, stats=ranks)
    }
    
    fgseaResTidy <- fgseaRes %>%
      as_tibble() %>%
      arrange(desc(NES))
    
    # Show in a nice table:
    # fgseaResTidy %>% 
    #   dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
    #   arrange(padj) %>% 
    #   DT::datatable()
    
    # Plot the normalized enrichment scores and color the bars indicating if the pathway was significant
    # unlist(strsplit(names(deseq.exp)[1], "_vs_"))[1]
    fgseaResTidy <- filter(fgseaResTidy, fgseaResTidy$padj < 0.05)
    
    fgseaResTidy$Direction <- recode(factor(fgseaResTidy$NES > 0), "TRUE" = "Upregulated", "FALSE" = "Downregulated")
    
    p <- ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
      geom_col(aes(fill=Direction)) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title=title) + 
      theme(axis.text.y = element_text(size = 4))
    # p + theme_prism() # HERE IS AN ISSUE WITH THE PACKAGE CURRENTLY SO CHECK BACK LATER FOR UPDATES 26.02.2024
    
    attr(p, "pname") <- tools::file_path_sans_ext(basename(pathways))
    
    # Extract genes from each fGSEA enriched pathway
    enriched.p %>% 
      enframe("pathway", "symbol") %>% 
      unnest() %>% 
      inner_join(res, by="symbol")
    
    # if(!missing(outfolder)){
    #   pname <- tools::file_path_sans_ext(basename(gmt.hallmark))
    #   pdf(here(outfolder,pname))
    #   print(p)
    #   dev.off()
    # }
    
    print("Done!")
    return(list(fgseaRes          = fgseaResTidy, 
                enriched.pathways = enriched.p,
                plot              = p))
  }
  
  if(is.list(deseq.exp)){
    res <- list()
    for(i in seq_along(deseq.exp)){
      title <- paste0("Enriched Pathways: ",names(deseq.exp)[i])
      
      exp.name <- names(deseq.exp)[i]
      print(paste0("Running fGSEA on ", exp.name, "..."))
      
      res.dir <- here::here(outfolder, dname, exp.name)
      print(res.dir)
      dir.create(res.dir, showWarnings = F, recursive = T)
      
      res[[exp.name]] <- runfGSEA(gene.list = deseq.exp[[exp.name]]$res, 
                                  pathways = gmt.file, outfolder = outfolder,
                                  nperm = nperm, title = title,
                                  verbose = verbose)
      
      # pdf(here(res.dir, paste0(attributes(res[[exp.name]]$plot)$pname,'.pdf')))
      # print(res[[exp.name]]$plot)
      # dev.off()
      
      ggsave(here(res.dir, paste0(attributes(res[[exp.name]]$plot)$pname,'.pdf')),
             res[[exp.name]]$plot, 
             width = 10, height = 6)
      
      saveRDS(res[[exp.name]], here(res.dir, paste0(attributes(res[[exp.name]]$plot)$pname, '_results.RDS')))
      
    }
  } else{
    res <- runfGSEA(gene.list = deseq.exp[[exp.name]]$res, 
                    pathways = gmt.file, outfolder = outfolder,
                    nperm = nperm, title = title,
                    verbose = verbose)
  }
  
  return(res)
}

pkgs <- c("here", "DESeq2", "org.Hs.eg.db", "EnhancedVolcano",
          "docopt", "ggprism",  "mesocore", "dplyr", "fgsea",
          "tibble", "tidyr")

suppressMessages(mesocore::handleRequirements(pkgs))

#---- Header ----
"Mesothelioma AI Pipeline - Differential Expression using DESeq2

Usage: run_deseq.R [options]

Options:
  -h --help                     Show this screen.
  -d --data_path=<PATH>         Path to the counts data set (RNA-seq data).
  -e --experiment_design=<PATH> Path to the table with experiment design.
  -o --outfolder=<PATH>         Folder where the results are saved.
  -a --alpha_value=<DOUBLE>     Set the threshold for the FDR value cutoff [Default: 0.1].
  -p --p_value=<DOUBLE>         Significance level for the DEG genes [Default: 0.05].
  -g --gmt_file=<PATH>          Path to the gene signature in the '.gmt' format.
  -n --nperm=<INT>              Number of permutations to run the GSEA on.
  -v --verbose=<BOOLEAN>        If set to 1, prints messages verbously.
  -V --version
"-> doc

#---- Arguments ----
arguments <- docopt(doc, quoted_args = TRUE, help = TRUE)
print(arguments)

data.path     <- arguments$data_path
exp.path      <- arguments$experiment_design
gmt.file      <- arguments$gmt_file
outfolder     <- arguments$outfolder
alpha.value   <- arguments$alpha_value
p.value       <- arguments$p_value
nperm         <- as.integer(arguments$nperm)
verbose       <- as.integer(arguments$verbose)

if(as.integer(arguments$verbose) == 1){
  verbose <- T
} else{
  verbose <- F
}

# 1. Generate a ranked list of genes
deseq.exp <- getRankedGeneList(data.path, exp.path, outfolder)

dname <- basename(dirname(exp.path))

# 2. Run Enrichment Analysis using fGSEA
enrichments <- getEnrichedPathways(experiment = deseq.exp, 
                                   pathways = gmt.file, 
                                   outfolder = outfolder,
                                   verbose = verbose,
                                   dname   = dname,
                                   nperm, 
                                   p.value)

#----- Main -----
# main <- function(data.path, exp.path, alpha.value, cts.filter,
#                  p.value, median.filter, contrasts,
#                  outfolder, do.contrasts){ 
# }
# main(data.path = data.path, exp.path = exp.path,
#      alpha.value = alpha.value,
#      p.value = p.value, do.contrasts = T, 
#      verbose = T, serialise = F
#      outfolder = outfolder, do.contrasts = T)
