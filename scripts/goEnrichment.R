library(tidyverse)
library(splitstackshape)
library(snakecase)
library(topGO)

#### Get data prepped for GO term enrichment: ####
# Read in the annotations from eggnog:
eggnogAnnotations <- read_delim("./04_keggAnnotations/results/keggAnnotations.emapper.annotations",
                                delim = "\t",
                                skip = 4) %>%
  filter(!is.na(seed_ortholog))

# Subset so we just have gene name and GO domain IDs:
longAnnotations <- dplyr::select(eggnogAnnotations,
                                 c("#query", "GOs"))

# Reshape into a long dataframe:
longAnnotations <- cSplit(longAnnotations, 
                          splitCols = "GOs",
                          sep = ",",
                          direction = "long")

# Take out genes without GO terms and get distinct rows:
longAnnotations <- filter(longAnnotations,
                          GOs != "-") %>%
  distinct()

# Rename the column #query to geneName
longAnnotations <- longAnnotations %>%
  dplyr::rename(geneName = `#query`)

# Convert transcripts to genes:
longAnnotations$geneName <- gsub(pattern = "-R.*", 
                                 replacement = "", 
                                 as.character(longAnnotations$geneName))

# Create list with element for each gene, containing vectors with all terms for each gene
wideListAnnotations <- tapply(longAnnotations$GOs, longAnnotations$geneName, function(x)x)

# Read in the differential expression information:
differentialExpression <- read_csv(file = "./finalResults/allDifferentialExpressionResultsAndFunctions.csv")

#### Do GO term enrichment: ####
# Write a function to do GO term enrichment for each DE contrast:
goTermsByContrast <- function(selectedContrast) {
  contrast <- snakecase::to_upper_camel_case(selectedContrast)
  resultsFile <- paste("./finalResults/",
                       contrast,
                       "_EnrichmentGOResults.tsv",
                       sep = "")
  if (file.exists(resultsFile)) {
    print("GO analysis has already been conducted.")
    enrichedGOTerms <- read_delim(file = resultsFile,
                                  delim = "\t",
                                  quote = "none")
    # Plotting:
    enrichedGOTerms$ratio <- enrichedGOTerms$Significant/enrichedGOTerms$Annotated
    enrichedGOTerms$raw.p.value <- as.numeric(enrichedGOTerms$raw.p.value)
    
    fullDotPlot <- ggplot(data = enrichedGOTerms) + 
      geom_point(mapping = aes(x = ratio, 
                               y = fct_reorder(Term, 
                                               ratio),
                               size = ratio, 
                               color = `raw.p.value`)) +
      theme_bw(base_size = 12) +
      scale_color_gradient2(high = "#F1A208", 
                            mid = "#3e8a04") +
      guides(colour = guide_colorbar(reverse = T)) +
      ylab(NULL) +
      xlab("Ratio of enrichment") +
      ggtitle(selectedContrast) +
      theme(axis.text = element_text(color = "black")) + 
      guides(colour = guide_colourbar(order = 1), 
             size = guide_legend(order = 2)) +
      labs(color = "p-value (FDR-adjusted)",
           size = "Ratio of enrichment")
    
    selectedEnrichedGOTerms <- enrichedGOTerms %>%
      dplyr::arrange(`raw.p.value`) %>%
      head(n = 15)
    
    dotPlot <- ggplot(data = selectedEnrichedGOTerms) + 
      geom_point(mapping = aes(x = ratio, 
                               y = fct_reorder(Term, 
                                               ratio),
                               size = ratio, 
                               color = `raw.p.value`)) +
      theme_bw(base_size = 12) +
      scale_color_gradient2(high = "#F1A208", 
                            mid = "#3e8a04") +
      guides(colour = guide_colorbar(reverse = T)) +
      ylab(NULL) +
      xlab("Ratio of enrichment") +
      ggtitle(selectedContrast) +
      theme(axis.text = element_text(color = "black")) + 
      guides(colour = guide_colourbar(order = 1), 
             size = guide_legend(order = 2)) +
      labs(color = "p-value (FDR-adjusted)",
           size = "Ratio of enrichment")
    
    return(dotPlot)
  } else {
    # Get significance info from the differential expression results:
    significanceInfo <- differentialExpression %>%
      filter(contrast == selectedContrast) %>%
      dplyr::select(gene_name, 
                    padj, 
                    log2FoldChange) %>%
      distinct()
    
    # Create a formatted list of significance information:
    geneList <- ifelse(significanceInfo$padj <= 0.05, 
                       1, 
                       0)
    # Give geneList names:
    names(geneList) <- significanceInfo$gene_name
    
    # Analyze Biological Process terms:
    GOdataBP <- new("topGOdata",
                    ontology = "BP",
                    allGenes = geneList,
                    geneSelectionFun = function(x)(x == 1),
                    annot = annFUN.gene2GO, 
                    gene2GO = wideListAnnotations)
    # Run Fisher's exact test to check for enrichment:
    resultsFisherBP <- runTest(GOdataBP, 
                               algorithm = "elim", 
                               statistic = "fisher")
    resultsFisherBP
    resultsFisherBPTable <- GenTable(GOdataBP, 
                                     raw.p.value = resultsFisherBP, 
                                     topNodes = length(resultsFisherBP@score),
                                     numChar = 120)
    
    # Analyze Molecular Function terms:
    GOdataMF <- new("topGOdata",
                    ontology = "MF",
                    allGenes = geneList,
                    geneSelectionFun = function(x)(x == 1),
                    annot = annFUN.gene2GO, 
                    gene2GO = wideListAnnotations)
    # Run Fisher's exact test to check for enrichment:
    resultFisherMF <- runTest(GOdataMF, 
                              algorithm = "elim", 
                              statistic = "fisher")
    resultFisherMF
    resultsFisherMFTable <- GenTable(GOdataMF, 
                                     raw.p.value = resultFisherMF, 
                                     topNodes = length(resultFisherMF@score),
                                     numChar = 120)
    
    # Analyze Cellular Component terms:
    GOdataCC <- new("topGOdata",
                    ontology = "CC",
                    allGenes = geneList,
                    geneSelectionFun = function(x)(x == 1),
                    annot = annFUN.gene2GO, 
                    gene2GO = wideListAnnotations)
    # Run Fisher's exact test to check for enrichment:
    resultFisherCC <- runTest(GOdataCC, 
                              algorithm = "elim", 
                              statistic = "fisher")
    resultFisherCC
    resultsFisherCCTable <- GenTable(GOdataCC, 
                                     raw.p.value = resultFisherCC, 
                                     topNodes = length(resultFisherCC@score),
                                     numChar = 120)
    
    # Combine all of the results:
    enrichedGOTerms <- rbind(resultsFisherBPTable, 
                             resultsFisherMFTable, 
                             resultsFisherCCTable) %>%
      dplyr::filter(raw.p.value <= 0.01)
    
    # Export the results to a csv:
    readr::write_delim(enrichedGOTerms,
                       file = resultsFile,
                       delim = "\t",
                       quote = "none")
    # Plotting:
    enrichedGOTerms$ratio <- enrichedGOTerms$Significant/enrichedGOTerms$Annotated
    enrichedGOTerms$raw.p.value <- as.numeric(enrichedGOTerms$raw.p.value)
    
    fullDotPlot <- ggplot(data = enrichedGOTerms) + 
      geom_point(mapping = aes(x = ratio, 
                               y = fct_reorder(Term, 
                                               ratio),
                               size = ratio, 
                               color = `raw.p.value`)) +
      theme_bw(base_size = 12) +
      scale_color_gradient2(high = "#F1A208", 
                            mid = "#3e8a04") +
      guides(colour = guide_colorbar(reverse = T)) +
      ylab(NULL) +
      xlab("Ratio of enrichment") +
      ggtitle(selectedContrast) +
      theme(axis.text = element_text(color = "black")) + 
      guides(colour = guide_colourbar(order = 1), 
             size = guide_legend(order = 2)) +
      labs(color = "p-value (FDR-adjusted)",
           size = "Ratio of enrichment")
    
    selectedEnrichedGOTerms <- enrichedGOTerms %>%
      dplyr::arrange(`raw.p.value`) %>%
      head(n = 15)
    
    dotPlot <- ggplot(data = selectedEnrichedGOTerms) + 
      geom_point(mapping = aes(x = ratio, 
                               y = fct_reorder(Term, 
                                               ratio),
                               size = ratio, 
                               color = `raw.p.value`)) +
      theme_bw(base_size = 12) +
      scale_color_gradient2(high = "#F1A208", 
                            mid = "#3e8a04") +
      guides(colour = guide_colorbar(reverse = T)) +
      ylab(NULL) +
      xlab("Ratio of enrichment") +
      ggtitle(selectedContrast) +
      theme(axis.text = element_text(color = "black")) + 
      guides(colour = guide_colourbar(order = 1), 
             size = guide_legend(order = 2)) +
      labs(color = "p-value (FDR-adjusted)",
           size = "Ratio of enrichment")
    
    return(dotPlot)
  }
}
possiblyGoTermsByContrast <- possibly(goTermsByContrast,
                                      otherwise = "error.")
contrasts <- unique(differentialExpression$contrast)
contrasts <- contrasts[!is.na(contrasts)]

goTermPlots <- purrr::map(contrasts,
                          possiblyGoTermsByContrast)

allGOplots <- patchwork::wrap_plots(goTermPlots, 
                                    ncol = 2) + 
  patchwork::plot_annotation(title = 'Enriched GO terms across differential expression contrasts',
                             theme = theme(plot.title = element_text(size = 19)))
plot(allGOplots)

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
    minimum <- min(plot[[x]][["data"]][["raw.p.value"]])
    return(minimum)
  }
  colorMins <- purrr::map(1:num_plots,
                          colorMin)
  colorMins <- min(unlist(colorMins))
  colorMax <- function(x) {
    maximum <- max(plot[[x]][["data"]][["raw.p.value"]])
    return(maximum)
  }
  colorMaxs <- purrr::map(1:num_plots,
                          colorMax)
  colorMaxs <- max(unlist(colorMaxs))
  
  plot & 
    xlim(minX - 0.01, 
         maxX + 0.01) & 
    scale_color_gradient2(high = "#F1A208", 
                         mid = "#e8dfa5",
                         low = "#3e8a04",
                         limits = c(colorMins, 
                                    colorMaxs),
                         midpoint = (colorMins + colorMaxs)/2) & 
    scale_size(limits = c(minX,
                          maxX)) 
}

# Plot all results with consistent scales and a single legend:
goPlotsList <- makePlotsConsistent(allGOplots) + 
  patchwork::plot_layout(guides = "collect",
                         axes = "collect") 

goPlotsList

dir.create("./images/")
ggsave(filename = "./images/allGOresults.png",
       width = 18, 
       height = 10, 
       units = "in",
       dpi = 900)

# Read in all of the GO term enrichment results as dataframes for ease of viewing:
goTermFiles <- list.files(path = "./finalResults",
                          pattern = "*_EnrichmentGOResults.tsv",
                          full.names = TRUE)

goTermResults <- read_tsv(goTermFiles,
                          id = "contrast") 
goTermResults$contrast <- str_split_i(goTermResults$contrast,
                                      pattern = "/",
                                      i = 3)
goTermResults$contrast <- str_split_i(goTermResults$contrast,
                                      pattern = "_",
                                      i = 1)

goTermResults$contrast <- gsub("([a-z])([A-Z])", 
                               "\\1 \\L\\2", 
                               goTermResults$contrast, 
                               perl = TRUE)

goTermResultsTable <- goTermResults %>%
  dplyr::select("contrast",
                "GO.ID",
                "Term",
                "Annotated",
                "Significant",
                "raw.p.value") %>%
  gt::gt()
gtsave(goTermResultsTable, 
       "allGOTerms.docx")


numberPupae <- length(which(goTermResults$contrast == "./finalResults/WorkerPupaeVsSoldierPupae_EnrichmentGOResults.tsv"))
numberAdults <- length(which(goTermResults$contrast == "./finalResults/AdultWorkersVsAdultSoldiers_EnrichmentGOResults.tsv"))
numberSoldiers <- length(which(goTermResults$contrast == "./finalResults/AdultSoldiersVsSoldierPupae_EnrichmentGOResults.tsv"))
numberWorkers <- length(which(goTermResults$contrast == "./finalResults/AdultWorkersVsWorkerPupae_EnrichmentGOResults.tsv"))









