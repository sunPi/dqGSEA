# DESeq2 Functions
# convertGeneType   <- function(gene.list, to.type = "symbol"){
#   #https://bioinformatics.stackexchange.com/questions/5229/converting-gene-symbol-to-ensembl-id-in-r
#   library(EnsDb.Hsapiens.v86)
#   # 1. Convert from ensembl.gene to gene.symbol
#   if(to.type == "symbol"){
#     converted <- ensembldb::select(EnsDb.Hsapiens.v86, keys= gene.list, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
#     return(converted)
#   }
#   if(to.type == "ensemble"){
#     # 2. Convert from gene.symbol to ensembl.gene
#     converted <- ensembldb::select(EnsDb.Hsapiens.v86, keys= gene.list, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
#     return(converted)
#   } else{
#     stop("Specified the wrong encode destination type, please use either 'symbol' or 'ensemble'!")
#   }
# }
# checkGeneEncoding <- function(gene.id, verbose){
#   # Function to test if a string matches the gene name pattern
#   isHugo <- function(name) {
#     # gene.pattern <- "^[A-Z0-9-]+$"
#     gene_pattern <- "^[A-Z]+[0-9]+([-A-Z0-9]*)?$"
#     return(grepl(gene_pattern, name))
#   }
#   
#   if (any(startsWith(gene.id, "ENSG"))) {
#     if(verbose){
#       print("Converting from ENSEMBL id's to HUGO symbols...")
#     }
#     gene.id <- convertGeneType(gene.id, to.type = "symbol")
#     return(gene.id)
#     
#   } else{
#     # stop("Your gene id's are not encoded in either HUGO or ENSEMBLE type. This function supports 
#     #      the conversion between these two only!")
#     warning("Your gene id's are not encoded in either HUGO or ENSEMBLE type. This function supports 
#          the conversion between these two only! THIS FUNCTION IS NOT YET COMPLETE, MIGHT BE OK...CHECK...")
#     return(list("SYMBOL" = gene.id))
#   }
# }
getExpPerms      <- function(coldata){
  condition <- as.vector(unique(coldata$condition))

  all_perms <- combinat::permn(condition)
  names(all_perms) <- lapply(all_perms, paste0, collapse = "_vs_")

  return(all_perms)
}
removeDuplicates <- function(x){
  x <- x[!duplicated(rownames(x)),]
  return(x)
}
getColData       <- function(exp.data, cts){
  
  # cts <- cts[rowSums(cts == 0) == 0, ] # Remove rows with 0 values
  
  if(length(exp.data[[1]]) != length(colnames(cts))){
    print("Detecting a difference in sample size between count data and experiment design data, conforming...")
    exp.data <- dplyr::filter(exp.data, exp.data[[1]] %in% colnames(cts))
  }
  
  coldata <- data.frame(row.names = colnames(cts),
                        "condition" = factor(exp.data[[2]]))
  return(coldata)
}
runDESeq         <- function(cts, coldata){
  dds <- DESeqDataSetFromMatrix(countData=cts, 
                                colData=coldata, 
                                design= ~ condition) # Set tidy = T if data in long format
  dds <- DESeq(dds)
  
  return(dds)
}
getDESeqRes      <- function(dds, coldata, contrast, p.value, outfolder, serialize){
  # Get genes that are differentially expressed
  if(!is.null(contrast)){
    res <- results(dds, contrast = contrast) %>% na.omit() %>% as.data.frame()
  } else{
    res <- results(dds) %>% na.omit() %>% as.data.frame()
  }
  
  title <- paste(levels(coldata$condition), collapse = " vs. ")
  
  # Filter out the significant genes based on p.value and order/sort them
  sig <- res[res$padj < p.value,]
  sig <- sig[order(sig$padj),]
  sig <- sig[order(sig$log2FoldChange, decreasing = T), ]
  sig <- sig %>% as.data.frame() %>% arrange(desc(log2FoldChange), desc(padj))
  top.10.genes <- rownames(sig[1:10,])
  
  # Prepare a dataframe ready to be serialized
  out.df <- list("res" = res,
                 "sig" = sig,
                 "top10" = as.data.frame(top.10.genes))
  
  out.df$res <- rownames_to_column(out.df$res, var = "gene_symbol")
  out.df$sig <- rownames_to_column(out.df$sig, var = "gene_symbol")
  
  if(serialize){
    writexl::write_xlsx(out.df, 
                        here(outfolder,paste0(gsub(" ", "_", title), '_deg.xlsx')))
  }
  
  if(nrow(sig) > 0){
    out.df$any.sig <- T
    print(paste0("DESeq2 analysis found ", nrow(sig), " significant genes!"))
    
  } else{
    out.df$any.sig <- F
    print("DESeq2 found 0 significant genes!")
  }
  
  # Plot the results via the Volcano Plot
  
  p <- EnhancedVolcano(res,
                       lab = rownames(res),
                       x = 'log2FoldChange',
                       y = 'pvalue',
                       title = title) + theme_prism()

  if(serialize){
    png(here(outfolder,paste0(gsub(" ", "_", title), '_volcano.png')),
        width = 2400, height = 1800, res = 300)
    print(p)
    dev.off()
  }
  
  out.df <- list(significant = out.df$sig,
                 res         = out.df$res,
                 vplot       = p,
                 any.sig     = out.df$any.sig)
  
  return(out.df)

}
getRankedGeneList <- function(cts.path, exp.path, outfolder,
                              p.value, do.contrasts, 
                              verbose, serialize, dname){
  cts <- readFile(cts.path)
  exp <- readFile(exp.path)
  # outfolder <- dirname(cts.path)
  
  # Gets the experiment name, assuming the folder in which the counts file is in is also the experiment name
  exp.name <- cts.path %>% 
    dirname() %>% 
    basename()

  # Remove duplicated values
  dup.marker <- all(duplicated(cts[,1]))
  if(dup.marker){
    print("Found duplicates, removing...")
    cts <- cts[!duplicated(cts[,1]),]
    rownames(cts) <- cts$gene_symbol
    cts <- cts[ ,-1]
    
  } else{
    print("Found no duplicates!")
  }
  
  # Assign gene_id as rownames
  row.names(cts) <- cts[,1]
  cts <- cts[,-1]
  
  # Round values just in case
  cts <- round(cts)
  
  # Create a 'coldata' object from the experiment design
  coldata <- getColData(exp, cts)
  
  # Create a DESeq2 matrix, filter the matrix for low counts (=cts.filter)
  # and run the DESeq() algorithm
  dds     <- runDESeq(cts, coldata)

  # Extract the results for all possible contrasts of the experiments
  if(do.contrasts){
    # Calculate all possible permutations of the experiment 
    # (e.g. a vs. b ~ b vs. a)
    all_perms <- getExpPerms(coldata)
    
    contrasts <- list()
    res <- list()

    for(i in seq_along(all_perms)){
      pname <- names(all_perms)[i]
      # dir.create(here::here(outfolder, dname, pname), showWarnings = F, recursive = T)
      
      print(paste0("Looking for DEGs in ", pname))
      res.dir <- here::here(outfolder, dname, pname)
      
      if(serialize){
        dir.create(res.dir, F, T)
        saveRDS(res, here::here(res.dir, "diffseq_exp.RDS"))
      }
      
      contrasts[[pname]] <- c(names(coldata), paste0(all_perms[[pname]]))
      res[[pname]]       <- getDESeqRes(dds = dds, coldata = coldata,
                                        contrast = contrasts[[pname]],
                                        outfolder = res.dir,
                                        p.value = p.value, 
                                        serialize = serialize)

        # png(here::here(outfolder, dname, pname, paste0('MAplot.png')))
        # DESeq2::plotMA(res[[pname]]$p.ma.data, ylim=c(-2,2))
        # dev.off()
        
        # png(here::here(outfolder, dname, pname, paste0('Vplot.png')))
        # print(res[[pname]]$vplot)
        # dev.off()
        
      }
    } else{
      dir.create(outfolder, F, T)
      res <- getDESeqRes(dds = dds, coldata = coldata, p.value = p.value, outfolder = outfolder, 
                         contrast = NULL, serialize = serialize)
    }
  
  
  # dir.create(here::here(outfolder, dname, pname), showWarnings = F, recursive = T)
  # writeLines(dname, here::here(outfolder, ".config.txt"))
  
  
  return(res)
}

# fGSEA Functions
getEnrichedPathways <- function(deseq.exp, pathways, nperm, p.value, 
                                outfolder, verbose, serialize, dname){
  runfGSEA <- function(gene.list, pathways, nperm, p.value = 0.05,
                       title = NULL, outfolder = NULL, verbose = F)
  {
    res <- gene.list %>% 
      dplyr::select(gene_symbol, stat) %>% 
      na.omit() %>% 
      distinct() %>% 
      group_by(gene_symbol) %>% 
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
      theme(text = element_text(size = 8))
    # p + theme_prism() # HERE IS AN ISSUE WITH THE PACKAGE CURRENTLY SO CHECK BACK LATER FOR UPDATES 26.02.2024
    
    attr(p, "pname") <- tools::file_path_sans_ext(basename(pathways))
    
    # Extract genes from each fGSEA enriched pathway
    enriched.p %>% 
      enframe("pathway", "gene_symbol") %>% 
      unnest() %>% 
      inner_join(res, by="gene_symbol")
    
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
      
      res[[exp.name]] <- runfGSEA(gene.list = deseq.exp[[exp.name]]$res, 
                                  pathways = gmt.file, outfolder = res.dir,
                                  nperm = nperm, title = title,
                                  verbose = verbose)
      
      if(serialize){
        res.dir <- here::here(outfolder, dname, exp.name)
        
        if(verbose){
          print(res.dir)
        }
        
        dir.create(res.dir, showWarnings = F, recursive = T)
        
        ggsave(here(res.dir, paste0(attributes(res[[exp.name]]$plot)$pname,'.png')),
               res[[exp.name]]$plot, 
               width = 5, height = 3)
        
        saveRDS(res[[exp.name]], here(res.dir, paste0(attributes(res[[exp.name]]$plot)$pname, '_results.RDS')))
        
        # pdf(here(res.dir, paste0(attributes(res[[exp.name]]$plot)$pname,'.pdf')))
        # print(res[[exp.name]]$plot)
        # dev.off()
      }
    }
  } else{
    res <- runfGSEA(gene.list = deseq.exp[[exp.name]]$res, 
                    pathways = gmt.file, outfolder = outfolder,
                    nperm = nperm, title = title,
                    verbose = verbose)
    
    if(serialize){
      dname <- tools::file_path_sans_ext(basename(gmt.file))
      res.dir <- here::here(outfolder, dname)
      if(verbose){
        print(res.dir)
      }
      
      dir.create(res.dir, showWarnings = F, recursive = T)
      
      ggsave(here(res.dir, paste0(attributes(res$plot)$pname,'.png')),
             res$plot, 
             width = 5, height = 3)
      
      saveRDS(res, here(res.dir, paste0(attributes(res$plot)$pname, '_results.RDS')))
      
      # pdf(here(res.dir, paste0(attributes(res[[exp.name]]$plot)$pname,'.pdf')))
      # print(res[[exp.name]]$plot)
      # dev.off()
    }
  }
  
  return(res)
}
# Main Function
main <- function(cts.path, exp.path, gmt.file,
                 outfolder, p.value, nperm, 
                 verbose, do.contrasts = T, serialize = T, deseq.exp = NULL){
  
  # dname <- basename(dirname(exp.path))
  dname <- tools::file_path_sans_ext(basename(gmt.file))
  # 1. Generate a ranked list of genes
  if(is.null(deseq.exp)){
    deseq.exp <- getRankedGeneList(cts.path     = cts.path, 
                                   exp.path     = exp.path, 
                                   outfolder    = outfolder, 
                                   p.value      = p.value, 
                                   verbose      = verbose, 
                                   do.contrasts = do.contrasts, 
                                   serialize    = serialize, 
                                   dname        = dname)
  }
  
  # 2. Run Enrichment Analysis using fGSEA
  enrichments <- getEnrichedPathways(deseq.exp = deseq.exp, 
                                     pathways  = gmt.file, 
                                     outfolder = outfolder,
                                     verbose   = verbose,
                                     nperm     = nperm, 
                                     p.value   = p.value,
                                     serialize = serialize,
                                     dname     = dname )
  
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
  -p --p_value=<DOUBLE>         Significance level for the DEG genes [Default: 0.05].
  -g --gmt_file=<PATH>          Path to the gene signature in the '.gmt' format.
  -n --nperm=<INT>              Number of permutations to run the GSEA on.
  -s --deseq_exp=<PATH>         Specify a ranked list of genes for GSEA.
  -v --verbose=<BOOLEAN>        If set to 1, prints messages verbously.
  -V --version
"-> doc

#---- Arguments ----
arguments <- docopt(doc, quoted_args = TRUE, help = TRUE)
print(arguments)

cts.path  <- normalizePath(arguments$data_path)
exp.path  <- arguments$experiment_design
gmt.file  <- arguments$gmt_file
outfolder <- normalizePath(arguments$outfolder)
p.value   <- as.double(arguments$p_value)
nperm     <- as.integer(arguments$nperm)
deseq.exp <- arguments$deseq_exp
verbose   <- as.integer(arguments$verbose)

if(as.integer(arguments$verbose) == 1){
  verbose <- T
} else{
  verbose <- F
}


#----- Main -----
# do.contrasts <- T
# serialize    <- T
# verbose      <- T
# p.value      <- 0.05
# nperm        <- 1000
# cts.path     <- "~/Documents/Cancer_Studies_PhD/Studies/Study_Biphasic/datasets/UvM/Cunniff_counts.csv"
# exp.path     <- "~/Documents/Cancer_Studies_PhD/Studies/Study_Biphasic/datasets/UvM/Cunniff_expdata.xlsx"
# outfolder    <- "~/Documents/Cancer_Studies_PhD/Studies/Study_Biphasic/results/UvM"
# gmt.file     <- "/home/jr453/bioinf-tools/db/molsigs/msigdb/h.all.v2023.2.Hs.symbols.gmt"
if(arguments$deseq_exp == "NULL"){
  deseq.exp <- NULL
}
if(!is.null(deseq.exp)){
  main(cts.path, exp.path, gmt.file,
       outfolder, p.value, nperm, 
       verbose, do.contrasts = T, serialize = T, deseq.exp = deseq.exp)
} else{
  main(cts.path, exp.path, gmt.file,
       outfolder, p.value, nperm, 
       verbose, do.contrasts = T, serialize = T, deseq.exp = NULL)
}


  