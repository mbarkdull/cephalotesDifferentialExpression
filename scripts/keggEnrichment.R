library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(data.table)
library(splitstackshape)
library(snakecase)

# Write a function to do KEGG enrichment for a selected differential expression contrast:
testKEGGenrichment <- function(selectedContrast) {
  #### Read in all of our genes and their KEGG annotations (from eggnog-mapper): ####
  keggAnnotations <- read_delim("keggAnnotations.emapper.annotations",
                                delim = "\t",
                                skip = 4) %>%
    filter(!is.na(seed_ortholog))
  
  # Get just the Cephalotes varians gene names and the KO annotations:
  genesToKO <- select(keggAnnotations,
                      c("#query", "KEGG_ko"))
  
  # Strip off the "ko" before each term:
  genesToKO$KEGG_ko <- gsub(pattern = "ko:", 
                            replacement = "", 
                            as.character(genesToKO$KEGG_ko))
  
  # Split the instances where a single gene has multiple KO annotations:
  genesToKO <- concat.split.multiple(genesToKO, 
                                     split.cols = "KEGG_ko", 
                                     seps = ",", 
                                     direction = "long") %>%
    distinct()  %>%
    dplyr::filter(KEGG_ko != "-")
  
  # Strip off the "-R*" after gene names:
  genesToKO$`#query` <- gsub(pattern = "-R.*", 
                             replacement = "", 
                             as.character(genesToKO$`#query`))
  
  # Change the column names in the terms-to-gene dictionary:
  colnames(genesToKO) <- c("gene", "term")
  
  # Re-order columns, because enricher wants terms first then genes:
  genesToKO <- genesToKO %>%
    dplyr::select("term",
                  "gene")
  
  #### Select just our genes of interest (i.e. differentially expressed genes:) ####
  pupalDEGenes <- read_csv(file = "allDifferentialExpressionResultsAndFunctions.csv") %>%
    filter(contrast == selectedContrast) %>%
    filter(padj <= 0.05)
  pupalDEGenes <- unique(pupalDEGenes$gene_name)
  
  #### Create a list of the KEGG terms I am interested in ####
  interestingKEGGTerms <- genesToKO %>%
    dplyr::filter(gene %in% pupalDEGenes) %>%
    unlist() %>%
    as.vector()
  
  #### Create a list of all of the KEGG terms in my gene set: ####
  background_genes <- genesToKO$gene
  backgroundKEGGTerms <- genesToKO %>%
    dplyr::filter(gene %in% background_genes) %>%
    unlist() %>%
    as.vector()
  
  #### Test for KEGG term enrichment ####
  enrichment_kegg <- enrichKEGG(interestingKEGGTerms,
                                organism = "ko",
                                keyType = "kegg",
                                pvalueCutoff = 0.05,
                                pAdjustMethod = "BH",
                                universe = backgroundKEGGTerms,
                                minGSSize = 10,
                                maxGSSize = 500,
                                qvalueCutoff = 0.05,
                                use_internal_data = FALSE)
  
  #save the enrichment result
  contrast <- snakecase::to_upper_camel_case(selectedContrast)
  resultsFile <- paste(contrast,
                       "_EnrichmentKEGGResults.csv",
                       sep = "")
  write.csv(file = resultsFile,                 
            x = enrichment_kegg@result)
  
  if (any(enrichment_kegg@result$p.adjust <= 0.05)){
    plotData <- enrichment_kegg@result
    plotData$geneRatioDecimal <- sapply(plotData$GeneRatio, function(x) eval(parse(text=x)))
    
    dotPlot <- ggplot(data = filter(plotData,
                                    `p.adjust` <= 0.05 &
                                      category != "Human Diseases")) + 
      geom_point(mapping = aes(x = geneRatioDecimal, 
                               y = fct_reorder(Description, 
                                               geneRatioDecimal),
                               size = geneRatioDecimal, 
                               color = `p.adjust`)) +
      theme_bw(base_size = 12) +
      scale_color_gradient2(high = "#e8d056", 
                            mid = "#0455b3") +
      guides(colour = guide_colorbar(reverse = T)) +
      ylab(NULL) +
      xlab("Fold-enrichment") +
      ggtitle(selectedContrast) +
      theme(axis.text = element_text(color = "black")) + 
      guides(colour = guide_colourbar(order = 1), 
             size = guide_legend(order = 2)) +
      labs(color = "p-value (FDR-adjusted)",
           size = "Fold-enrichment")
    
    return(dotPlot)
    
  }
  
}
possiblyTestKEGGenrichment <- possibly(testKEGGenrichment,
                                       otherwise = "Error.")

# List all the contrasts 
contrasts <- read_csv(file = "allDifferentialExpressionResultsAndFunctions.csv")
contrasts <- unique(contrasts$contrast)
contrasts <- contrasts[!is.na(contrasts)]

allKEGGenrichments <- purrr::map(contrasts,
                                 possiblyTestKEGGenrichment)
allKEGGplots <- patchwork::wrap_plots(allKEGGenrichments, 
                                      ncol = 2) + 
  patchwork::plot_annotation(title = 'Enriched KEGG pathways across differential expression contrasts',
                             theme = theme(plot.title = element_text(size = 19)))
plot(allKEGGplots)

# Write a function to make scales and axes consistent across plots:
makePlotsConsistent <- function(plot){
  # Get the total number of plots that are combined together:
  num_plots <- length(plot)
  
  # Fix x limits and bubble sizes:
  # Get the minimum and maximum values of geneRatioDecimal for each plot:
  xLimits <- lapply(1:num_plots, function(x) ggplot_build(plot[[x]])$layout$panel_scales_x[[1]]$range$range)
  # Get the minimum and maximum x values for each plot:
  minX <- min(unlist(xLimits))
  maxX <- max(unlist(xLimits))
  
  # Fix color scale:
  # Get the minimum and maximum values of p.adjust for each plot:
  colorMin <- function(x) {
    minimum <- min(plot[[x]][["data"]][["p.adjust"]])
    return(minimum)
  }
  colorMins <- purrr::map(1:num_plots,
                          colorMin)
  colorMins <- min(unlist(colorMins))
  colorMax <- function(x) {
    maximum <- max(plot[[x]][["data"]][["p.adjust"]])
    return(maximum)
  }
  colorMaxs <- purrr::map(1:num_plots,
                          colorMax)
  colorMaxs <- max(unlist(colorMaxs))
  
  plot & 
    xlim(minX - 0.01, 
         maxX + 0.01) & 
    scale_color_gradient(high = "#e8d056", 
                         low = "#0455b3",
                         limits = c(colorMins, 
                                    colorMaxs)) & 
    scale_size(limits = c(minX,
                         maxX)) 
}

# Plot all results with consistent scales and a single legend:
makePlotsConsistent(allKEGGplots) + 
  patchwork::plot_layout(guides = "collect",
                         axes = "collect")

ggsave(filename = "allKEGGresults.png",
       width = 16, 
       height = 10, 
       units = "in")
