library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(data.table)
library(splitstackshape)
library(snakecase)

#### Get a list of the contrasts to explore: ####
contrasts <- read_csv(file = "./finalResults/allDifferentialExpressionResultsAndFunctions.csv")
contrasts <- unique(contrasts$contrast)
contrasts <- contrasts[!is.na(contrasts)]

#### Write a function that does KEGG enrichment: #### 
doKeggEnrichment <- function(selectedContrast) {
  #### Read in all of our genes and their KEGG annotations (from eggnog-mapper): ####
  keggAnnotations <- read_delim("./04_keggAnnotations/results/keggAnnotations.emapper.annotations",
                                delim = "\t",
                                skip = 4) %>%
    filter(!is.na(seed_ortholog))
  
  # Get just the Cephalotes varians gene names and the KO annotations:
  genesToKO <- dplyr::select(keggAnnotations,
                             c("#query", "KEGG_ko"))
  
  # Strip off the "ko" before each term:
  genesToKO$KEGG_ko <- gsub(pattern = "ko:", 
                            replacement = "", 
                            as.character(genesToKO$KEGG_ko))
  
  # Split the instances where a single gene has multiple KO annotations:
  genesToKO <- splitstackshape::concat.split.multiple(genesToKO, 
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
  
  # Identify just our genes of interest (i.e. differentially expressed genes:) 
  DEGenes <- read_csv(file = "./finalResults/allDifferentialExpressionResultsAndFunctions.csv") %>%
    filter(contrast == selectedContrast) %>%
    filter(padj <= 0.05)
  DEGenes <- unique(DEGenes$gene_name)
  
  #### Create a list of the KEGG terms I am interested in ####
  interestingKEGGTerms <- genesToKO %>%
    dplyr::filter(gene %in% DEGenes) %>%
    unlist() %>%
    as.vector()
  
  #### Create a list of all of the KEGG terms in my gene set: ####
  background_genes <- genesToKO$gene
  backgroundKEGGTerms <- genesToKO %>%
    dplyr::filter(gene %in% background_genes) %>%
    unlist() %>%
    as.vector()
  
  #### Test for KEGG term enrichment ####
  enrichment_kegg <- clusterProfiler::enrichKEGG(interestingKEGGTerms,
                                                 organism = "ko",
                                                 keyType = "kegg",
                                                 pvalueCutoff = 0.05,
                                                 pAdjustMethod = "BH",
                                                 universe = backgroundKEGGTerms,
                                                 minGSSize = 10,
                                                 maxGSSize = 500,
                                                 qvalueCutoff = 0.05,
                                                 use_internal_data = FALSE)
  enrichment_kegg@result$contrast <- selectedContrast
  return(enrichment_kegg)
}

allKEGGs <- purrr::map(contrasts,
                       doKeggEnrichment)
extractKEGGS <- function(index) {
  dataframe <- allKEGGs[[index]]@result
  dataframe$index <- index
  return(dataframe)
}

keggDataframes <- purrr::map(1:4,
                             extractKEGGS)
keggDataframes <- as.data.frame(do.call(rbind, keggDataframes)) 

keggDataframes$p.adjust <- as.numeric(keggDataframes$p.adjust)
keggDataframes <- keggDataframes %>%
  filter(p.adjust <= 0.05) 

ggplot(data = keggDataframes) + 
  geom_bar(mapping = aes(x = subcategory)) + 
  facet_wrap(~ contrast) + 
  theme(axis.text = element_text(angle = 90))


repeatedKEGGpathways <- keggDataframes %>% 
  group_by(ID) %>% 
  filter(n() > 1) %>% 
  filter(category != "Human Diseases") %>%
  select(c("index",
           "contrast",
           "category",
           "subcategory",
           "ID",
           "Description",
           "p.adjust"))

repeatedKEGGpathways <- repeatedKEGGpathways %>% 
  group_by(ID) %>% 
  mutate(combinedIndex = paste0(index,
                                collapse = ", ")) %>%
  select(-c(index, contrast, p.adjust)) %>%
  distinct() %>%
  ungroup()

ontogenyRepeats <- repeatedKEGGpathways %>%
  filter(grepl("3", repeatedKEGGpathways$combinedIndex) & grepl("4", repeatedKEGGpathways$combinedIndex))

morphRepeats <- repeatedKEGGpathways %>%
  filter(grepl("1", repeatedKEGGpathways$combinedIndex) & grepl("2", repeatedKEGGpathways$combinedIndex))

#### Write a function to run the KEGG enrichment and generate a dot plot: ####
generateDotPlot <- function(selectedContrast) {
  enrichment_kegg <- doKeggEnrichment(selectedContrast)
  
  # Save the enrichment result
  dir.create("./finalResults/")
  contrast <- snakecase::to_upper_camel_case(selectedContrast)
  plotData <- enrichment_kegg@result
  plotData$geneRatioDecimal <- sapply(plotData$GeneRatio, function(x) eval(parse(text=x)))
  significantResults <- plotData %>%
    dplyr::filter(`p.adjust` <= 0.05)
  
  resultsFile <- paste("./finalResults/",
                       contrast,
                       "_EnrichmentKEGGResults.csv",
                       sep = "")
  write.csv(file = resultsFile,                 
            x = significantResults)
  
  # Generate a dotplot:
  if (any(enrichment_kegg@result$p.adjust <= 0.05)){
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
possiblyGenerateDotPlot <- possibly(generateDotPlot,
                                    otherwise = "error")

allDotPlots <- purrr::map(contrasts,
                          possiblyGenerateDotPlot)

allDotPlots <- patchwork::wrap_plots(allDotPlots, 
                                     ncol = 2) + 
  patchwork::plot_annotation(title = 'Enriched KEGG pathways across differential expression contrasts',
                             theme = theme(plot.title = element_text(size = 19)))
plot(allDotPlots)

# Plot results together:
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
makePlotsConsistent(allDotPlots) + 
  patchwork::plot_layout(guides = "collect",
                         axes = "collect")

dir.create("./images/")
ggsave(filename = "./images/allKEGGresults.png",
       width = 16, 
       height = 10, 
       units = "in")

#### Write a function to run the KEGG enrichment and generate a map plot: ####
generateMapPlot <- function(selectedContrast) {
  enrichment_kegg <- doKeggEnrichment(selectedContrast)
  # Generate an enrichment map:
  # Get the pairwise similarity matrix for the KEGG terms:
  pairwiseSimilarity <- pairwise_termsim(enrichment_kegg,
                                         showCategory = 200)
  # Get a vector of terms that are not human disease terms, so those are not displayed:
  diseaseTerms <- select(enrichment_kegg@result,
                         c(Description, category)) %>%
    filter(category == "Human Diseases")
  diseaseTerms <- diseaseTerms$Description
  
  termsToDisplay <- rownames(pairwiseSimilarity@termsim)
  termsToDisplay <- setdiff(termsToDisplay, diseaseTerms)
  
  # Make the plot:
  mapPlot <- emapplot(pairwiseSimilarity,
                      showCategory = termsToDisplay, 
                      cex.params = list(line = 0.25,
                                        category_node = 1,
                                        category_label = 0.5)) +
    set_enrichplot_color(colors = c("#0455b3",
                                    "#e8d056"),
                         type = "fill") +
    ggtitle(selectedContrast) +
    labs(color = "p-value (FDR-adjusted)",
         size = "Number of genes")
  plot(mapPlot)
  return(mapPlot)
}

possiblyGenerateMapPlot <- possibly(generateMapPlot,
                                    otherwise = "error")

allMapPlots <- purrr::map(contrasts,
                          possiblyGenerateMapPlot)

allMapPlotsPatchwork <- patchwork::wrap_plots(allMapPlots, 
                                              ncol = 2) + 
  patchwork::plot_annotation(title = 'Enriched KEGG pathways across differential expression contrasts',
                             theme = theme(plot.title = element_text(size = 19))) 
plot(allMapPlotsPatchwork)

# Plot results together:
# Write a function to make scales and axes consistent across plots: 
makeMapPlotsConsistent <- function(plot){
  # Get the total number of plots that are combined together:
  num_plots <- length(plot)
  
  # Fix x limits: ####
  # Get the minimum and maximum values of geneRatioDecimal for each plot:
  xLimits <- lapply(1:num_plots, function(x) ggplot_build(plot[[x]])$layout$panel_scales_x[[1]]$range$range)
  # Get the minimum and maximum x values for each plot:
  minX <- min(unlist(xLimits))
  maxX <- max(unlist(xLimits))
  
  # Fix color scale: ####
  # Get the minimum and maximum values of color for each plot:
  colorMin <- function(x) {
    minimum <- min(plot[[x]][["data"]][["color"]])
    return(minimum)
  }
  colorMins <- purrr::map(1:num_plots,
                          colorMin)
  colorMins <- min(unlist(colorMins))
  colorMax <- function(x) {
    maximum <- max(plot[[x]][["data"]][["color"]])
    return(maximum)
  }
  colorMaxs <- purrr::map(1:num_plots,
                          colorMax)
  colorMaxs <- max(unlist(colorMaxs))
  
  # Fix the bubble sizes: ####
  # Get the minimum and maximum values of size for each plot:
  sizeMin <- function(x) {
    minimum <- min(plot[[x]][["data"]][["size"]])
    return(minimum)
  }
  sizeMins <- purrr::map(1:num_plots,
                         sizeMin)
  sizeMins <- min(unlist(sizeMins))
  sizeMax <- function(x) {
    maximum <- max(plot[[x]][["data"]][["size"]])
    return(maximum)
  }
  sizeMaxs <- purrr::map(1:num_plots,
                         sizeMax)
  sizeMaxs <- max(unlist(sizeMaxs))
  
  plot & 
    xlim(minX - 0.01, 
         maxX + 0.01) &
    set_enrichplot_color(colors = c("#0455b3",
                                    "#e8d056"),
                         limits = c(colorMins, 
                                    colorMaxs),
                         type = "fill") & 
    scale_size(limits = c(sizeMins,
                          sizeMaxs)) 
}

# Plot all results with consistent scales and a single legend:
makeMapPlotsConsistent(allMapPlotsPatchwork) + 
  patchwork::plot_layout(guides = "collect",
                         axes = "collect") & 
  theme(panel.border = element_rect(colour = "black", 
                                    fill = NA))

# Save the plot:
ggsave(filename = "./images/allMapPlotsPatchwork.png",
       width = 16, 
       height = 10, 
       units = "in")

#### Generate a ridgeline plot: ####

generateRidgelinePlot <- function(selectedContrast) {
  enrichment_kegg <- doKeggEnrichment(selectedContrast)
  
  # Get a key that links signficant KEGG pathways to KEGG terms:
  signficantKEGGPathToTerm <- enrichment_kegg@result %>%
    filter(p.adjust < 0.05) 
  signficantKEGGPathToTerm <- splitstackshape::concat.split.multiple(signficantKEGGPathToTerm, 
                                                    split.cols = "geneID", 
                                                    seps = "/", 
                                                    direction = "long")
  
  
  # Get a key that links KEGG terms to genes:
  keggTermToGene <- read_delim("./04_keggAnnotations/results/keggAnnotations.emapper.annotations",
                               delim = "\t",
                               skip = 4) %>%
    filter(!is.na(seed_ortholog))
  # Get just the Cephalotes varians gene names and the KO annotations:
  keggTermToGene <- dplyr::select(keggTermToGene,
                                  c("#query", "KEGG_ko"))
  # Strip off the "ko" before each term:
  keggTermToGene$KEGG_ko <- gsub(pattern = "ko:", 
                                 replacement = "", 
                                 as.character(keggTermToGene$KEGG_ko))
  # Split the instances where a single gene has multiple KO annotations:
  keggTermToGene <- splitstackshape::concat.split.multiple(keggTermToGene, 
                                          split.cols = "KEGG_ko", 
                                          seps = ",", 
                                          direction = "long") %>%
    distinct()  %>%
    dplyr::filter(KEGG_ko != "-")
  # Strip off the "-R*" after gene names:
  keggTermToGene$`#query` <- gsub(pattern = "-R.*", 
                                  replacement = "", 
                                  as.character(keggTermToGene$`#query`))
  # Change the column names in the terms-to-gene dictionary:
  colnames(keggTermToGene) <- c("gene", "term")
  
  # Get key linking genes to their expression: 
  genesToExpression <- read_csv(file = "./finalResults/allDifferentialExpressionResultsAndFunctions.csv") %>%
    filter(contrast == selectedContrast) %>%
    select(-c(seq.name)) %>%
    distinct()
  
  # Join those three keys together:
  ridgelineData <- full_join(keggTermToGene,
                             signficantKEGGPathToTerm,
                             by = c("term" = "geneID"),
                             relationship = "many-to-many")
  ridgelineData <- full_join(ridgelineData,
                             genesToExpression,
                             by = c("gene" = "gene_name")) 
  ridgelineData <- ridgelineData %>%
    select("gene",
           "term",
           "category",
           "Description",
           "p.adjust",
           "log2FoldChange",
           "padj")
  colnames(ridgelineData) <- c("Gene",
                               "Term",
                               "category",
                               "KEGG pathway",
                               "Adjusted p-value",
                               "Fold change",
                               "padj")
  ridgelineData <- ridgelineData %>%
    dplyr::arrange(`KEGG pathway`)
  ridgelineData$`KEGG pathway` <- as.factor(ridgelineData$`KEGG pathway`)
  
  # Make the ridgeline plot:
  ggplot(data = filter(ridgelineData,
                       !is.na(`KEGG pathway`),
                       category != "Human Diseases"),
         mapping = aes(x = `Fold change`, 
                       y = `KEGG pathway`,
                       fill = `Adjusted p-value`)) +
    ggridges::geom_density_ridges(linewidth = 0.15) +
    scale_fill_gradientn(colors = c("#0455b3", 
                                    "#e8d056")) +
    theme_bw() +
    labs(title = selectedContrast)
}
possiblyGenerateRidgelinePlot <- possibly(generateRidgelinePlot,
                                          otherwise = "error")
  
listOfRidgelinePlots <- purrr::map(contrasts,
                                   possiblyGenerateRidgelinePlot)

allRidgelinePlotsPatchwork <- patchwork::wrap_plots(listOfRidgelinePlots, 
                                                    ncol = 2) + 
  patchwork::plot_annotation(title = 'Enriched KEGG pathways across differential expression contrasts',
                             theme = theme(plot.title = element_text(size = 19))) 
plot(allRidgelinePlotsPatchwork)

# Plot results together:
# Write a function to make scales and axes consistent across plots: 
makeRidgelinePlotsConsistent <- function(plot){
  # Get the total number of plots that are combined together:
  num_plots <- length(plot)
  
  # Fix x limits: ####
  # Get the minimum and maximum values of geneRatioDecimal for each plot:
  xLimits <- lapply(1:num_plots, function(x) ggplot_build(plot[[x]])$layout$panel_scales_x[[1]]$range$range)
  # Get the minimum and maximum x values for each plot:
  minX <- min(unlist(xLimits))
  maxX <- max(unlist(xLimits))
  
  # Fix color scale: ####
  # Get the minimum and maximum values of color for each plot:
  colorMin <- function(x) {
    minimum <- min(plot[[1]][["data"]][["Adjusted p-value"]])
    return(minimum)
  }
  colorMins <- purrr::map(1:num_plots,
                          colorMin)
  colorMins <- min(unlist(colorMins))
  colorMax <- function(x) {
    maximum <- max(plot[[1]][["data"]][["Adjusted p-value"]])
    return(maximum)
  }
  colorMaxs <- purrr::map(1:num_plots,
                          colorMax)
  colorMaxs <- max(unlist(colorMaxs))
  
  plot & 
    xlim(minX - 0.01, 
         maxX + 0.01) & 
    scale_fill_gradient(high = "#e8d056", 
                         low = "#0455b3",
                         limits = c(colorMins, 
                                    colorMaxs))
}

# Plot all results with consistent scales and a single legend:
makeRidgelinePlotsConsistent(allRidgelinePlotsPatchwork) + 
  patchwork::plot_layout(guides = "collect",
                         axes = "collect")

# Save the plot:
ggsave(filename = "./images/allRidgelinePlotsPatchwork.png",
       width = 16, 
       height = 10, 
       units = "in")








#### Explore results numerically: ####
# Read in all results as a dataframe:
keggTermFiles <- list.files(path = "./finalResults",
                            pattern = "*_EnrichmentKEGGResults.csv",
                            full.names = TRUE)

keggTermResults <- read_csv(keggTermFiles,
                            id = "contrast") %>%
  dplyr::filter(category != "Human Diseases") %>%
  distinct()

numberPupae <- length(which(keggTermResults$contrast == "./finalResults/WorkerPupaeVsSoldierPupae_EnrichmentKEGGResults.csv"))
numberAdults <- length(which(keggTermResults$contrast == "./finalResults/AdultWorkersVsAdultSoldiers_EnrichmentKEGGResults.csv"))
numberSoldiers <- length(which(keggTermResults$contrast == "./finalResults/AdultSoldiersVsSoldierPupae_EnrichmentKEGGResults.csv"))
numberWorkers <- length(which(keggTermResults$contrast == "./finalResults/AdultWorkersVsWorkerPupae_EnrichmentKEGGResults.csv"))

#### Explore specific pathways: ####
# List all of the Hippo pathway components from KEGG:
hippoComponents <- "K04382/K04662/K21283/K06269/K06630/K16197/K04671/K05692/K05689/K16686/K04663/K16620/K16621/K02375/K04491/K19631/K09448/K04503/K10151/K10152/K16175/K04676/K00182/K00408/K00444/K00572/K08959/K08960/K04500/K06069/K18952"
hippoComponents <- strsplit(hippoComponents, "/")[[1]]

# Write a function to find and extract matching genes from the differential expression results:
keggAnnotations <- read_delim("./04_keggAnnotations/results/keggAnnotations.emapper.annotations",
                              delim = "\t",
                              skip = 4) %>%
  filter(!is.na(seed_ortholog))

checkHippoExpression <- function(ko) {
  hippoComponent <- keggAnnotations %>% 
    filter(str_detect(KEGG_ko, 
                      ko))
  cvarGene <- hippoComponent$`#query`[1] %>%
    str_split_i(pattern = "-",
                i = 1)
  allResultsCombined <- read_csv("./finalResults/allDifferentialExpressionResultsAndFunctions.csv")
  expressionOfGene <- allResultsCombined %>% 
    filter(str_detect(gene_name, 
                      cvarGene) &
             contrast == "Worker pupae vs. soldier pupae")
  return(expressionOfGene)
  
}
possiblyCheckHippoExpression <- possibly(checkHippoExpression, otherwise = "error")

hippoGenes <- purrr::map(hippoComponents,
                         possiblyCheckHippoExpression)
hippoGenes <- as.data.frame(do.call(rbind, hippoGenes))

# Add in Hyphy results:
hyphyResults <- read_tsv("./hyphyResultsForDE.tsv")
hyphyResults$HOGmembers <- gsub(pattern = "_R.*$",
                                replacement = "",
                                x = hyphyResults$HOGmembers)
hippoGenes <- left_join(hippoGenes,
                        hyphyResults,
                        by = c("gene_name" = "HOGmembers"))

getTextSelectiveRegimes <- function(dataframe) {
  dataframe <- dataframe %>% mutate(positiveSelectionOn =
                                      case_when(as.numeric(as.character(testPvalueFDR)) <= 0.05 & 
                                                  as.numeric(as.character(backgroundPvalueFDR)) > 0.05 &
                                                  as.numeric(as.character(differencePvalueFDR)) <= 0.05 ~ "ForegroundOnly",
                                                
                                                as.numeric(as.character(testPvalueFDR)) <= 0.05 & 
                                                  as.numeric(as.character(backgroundPvalueFDR)) <= 0.05 &
                                                  as.numeric(as.character(differencePvalueFDR)) <= 0.05 ~ "SelectionOnBothButDifferent",
                                                
                                                as.numeric(as.character(testPvalueFDR)) <= 0.05 & 
                                                  as.numeric(as.character(backgroundPvalueFDR)) <= 0.05 &
                                                  as.numeric(as.character(differencePvalueFDR)) > 0.05 ~ "SelectionOnBothButNoSignificantDifference",
                                                
                                                as.numeric(as.character(testPvalueFDR)) <= 0.05 & 
                                                  as.numeric(as.character(backgroundPvalueFDR)) > 0.05 &
                                                  as.numeric(as.character(differencePvalueFDR)) > 0.05 ~ "EvidenceOfSelectionAssociatedWithTraitButNS",
                                                
                                                as.numeric(as.character(testPvalueFDR)) > 0.05 & 
                                                  as.numeric(as.character(backgroundPvalueFDR)) <= 0.05 &
                                                  as.numeric(as.character(differencePvalueFDR)) <= 0.05 ~ "BackgroundOnly",
                                                
                                                as.numeric(as.character(testPvalueFDR)) > 0.05 & 
                                                  as.numeric(as.character(backgroundPvalueFDR)) <= 0.05 &
                                                  as.numeric(as.character(differencePvalueFDR)) > 0.05 ~ "EvidenceOfSelectionAssociatedWithLackOfTraitButNS",
                                                
                                                as.numeric(as.character(testPvalueFDR)) > 0.05 & 
                                                  as.numeric(as.character(backgroundPvalueFDR)) > 0.05 &
                                                  as.numeric(as.character(differencePvalueFDR)) <= 0.05 ~ "NoEvidenceOfSelection",
                                                
                                                as.numeric(as.character(testPvalueFDR)) > 0.05 & 
                                                  as.numeric(as.character(backgroundPvalueFDR)) > 0.05 &
                                                  as.numeric(as.character(differencePvalueFDR)) > 0.05 ~ "NoEvidenceOfSelection"))
  dataframe <- dataframe %>%
    mutate(selectiveIntensity = case_when(pValueFDR < 0.05 & kValue < 1 ~ "selectionRelaxed",
                                          pValueFDR < 0.05 & kValue > 1 ~ "selectionIntensified",
                                          pValueFDR >= 0.05 ~ "noShift"))
  return(dataframe)
}

hippoGenes <- getTextSelectiveRegimes(hippoGenes)


uniqueHippoGenes <- hippoGenes %>%
  select(-c("seq.name",
            "inputFile.x",
            "testPvalue",
            "backgroundPvalue",
            "differencePvalue",
            "inputFile.y",
            "gene_name.y",
            "pValue")) %>%
  distinct()









test <- doKeggEnrichment("Worker pupae vs. soldier pupae")
test2 <- test@result

# List all of the Hippo pathway components from KEGG:
dopaminergicComponents <- "K07882/K07299/K07593/K04131/K01539/K18211/K04632/K02677/K19662/K19663/K04936/K12407/K04634/K04635/K04526/K08049/K01540"
dopaminergicComponents <- strsplit(dopaminergicComponents, "/")[[1]]

# Write a function to find and extract matching genes from the differential expression results:
keggAnnotations <- read_delim("./04_keggAnnotations/results/keggAnnotations.emapper.annotations",
                              delim = "\t",
                              skip = 4) %>%
  filter(!is.na(seed_ortholog))

checkHippoExpression <- function(ko) {
  hippoComponent <- keggAnnotations %>% 
    filter(str_detect(KEGG_ko, 
                      ko))
  cvarGene <- hippoComponent$`#query`[1] %>%
    str_split_i(pattern = "-",
                i = 1)
  allResultsCombined <- read_csv("./finalResults/allDifferentialExpressionResultsAndFunctions.csv")
  expressionOfGene <- allResultsCombined %>% 
    filter(str_detect(gene_name, 
                      cvarGene) &
             contrast == "Worker pupae vs. soldier pupae")
  return(expressionOfGene)
  
}
possiblyCheckHippoExpression <- possibly(checkHippoExpression, otherwise = "error")

dopaGenes <- purrr::map(dopaminergicComponents,
                         possiblyCheckHippoExpression)
dopaGenes <- as.data.frame(do.call(rbind, dopaGenes))

# Add in Hyphy results:
hyphyResults <- read_tsv("./hyphyResultsForDE.tsv")
hyphyResults$HOGmembers <- gsub(pattern = "_R.*$",
                                replacement = "",
                                x = hyphyResults$HOGmembers)
dopaGenes <- left_join(dopaGenes,
                        hyphyResults,
                        by = c("gene_name" = "HOGmembers"))

getTextSelectiveRegimes <- function(dataframe) {
  dataframe <- dataframe %>% mutate(positiveSelectionOn =
                                      case_when(as.numeric(as.character(testPvalueFDR)) <= 0.05 & 
                                                  as.numeric(as.character(backgroundPvalueFDR)) > 0.05 &
                                                  as.numeric(as.character(differencePvalueFDR)) <= 0.05 ~ "ForegroundOnly",
                                                
                                                as.numeric(as.character(testPvalueFDR)) <= 0.05 & 
                                                  as.numeric(as.character(backgroundPvalueFDR)) <= 0.05 &
                                                  as.numeric(as.character(differencePvalueFDR)) <= 0.05 ~ "SelectionOnBothButDifferent",
                                                
                                                as.numeric(as.character(testPvalueFDR)) <= 0.05 & 
                                                  as.numeric(as.character(backgroundPvalueFDR)) <= 0.05 &
                                                  as.numeric(as.character(differencePvalueFDR)) > 0.05 ~ "SelectionOnBothButNoSignificantDifference",
                                                
                                                as.numeric(as.character(testPvalueFDR)) <= 0.05 & 
                                                  as.numeric(as.character(backgroundPvalueFDR)) > 0.05 &
                                                  as.numeric(as.character(differencePvalueFDR)) > 0.05 ~ "EvidenceOfSelectionAssociatedWithTraitButNS",
                                                
                                                as.numeric(as.character(testPvalueFDR)) > 0.05 & 
                                                  as.numeric(as.character(backgroundPvalueFDR)) <= 0.05 &
                                                  as.numeric(as.character(differencePvalueFDR)) <= 0.05 ~ "BackgroundOnly",
                                                
                                                as.numeric(as.character(testPvalueFDR)) > 0.05 & 
                                                  as.numeric(as.character(backgroundPvalueFDR)) <= 0.05 &
                                                  as.numeric(as.character(differencePvalueFDR)) > 0.05 ~ "EvidenceOfSelectionAssociatedWithLackOfTraitButNS",
                                                
                                                as.numeric(as.character(testPvalueFDR)) > 0.05 & 
                                                  as.numeric(as.character(backgroundPvalueFDR)) > 0.05 &
                                                  as.numeric(as.character(differencePvalueFDR)) <= 0.05 ~ "NoEvidenceOfSelection",
                                                
                                                as.numeric(as.character(testPvalueFDR)) > 0.05 & 
                                                  as.numeric(as.character(backgroundPvalueFDR)) > 0.05 &
                                                  as.numeric(as.character(differencePvalueFDR)) > 0.05 ~ "NoEvidenceOfSelection"))
  dataframe <- dataframe %>%
    mutate(selectiveIntensity = case_when(pValueFDR < 0.05 & kValue < 1 ~ "selectionRelaxed",
                                          pValueFDR < 0.05 & kValue > 1 ~ "selectionIntensified",
                                          pValueFDR >= 0.05 ~ "noShift"))
  return(dataframe)
}

dopaGenes <- getTextSelectiveRegimes(dopaGenes)


uniqueDopaGenes <- dopaGenes %>%
  select(-c("seq.name",
            "inputFile.x",
            "testPvalue",
            "backgroundPvalue",
            "differencePvalue",
            "inputFile.y",
            "gene_name.y",
            "pValue")) %>%
  distinct()


ggplot(data = uniqueDopaGenes) +
  geom_bar(mapping = aes(x = positiveSelectionOn))
