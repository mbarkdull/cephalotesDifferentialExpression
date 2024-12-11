library(ggplot2)
library(ComplexUpset)
library(gt)
library(tidyverse)

#### Read in and combine the HYPHY results: ####
hyphyResults <- read_tsv("./hyphyResultsForDE.tsv")

# Read in the differential expression results:
differentialExpression <- read_csv(file = "./finalResults/allDifferentialExpressionResultsAndFunctions.csv")

# Combine with the Hyphy results:
allResultsCombined <- full_join(hyphyResults,
                                differentialExpression,
                                by = c("gene_name" = "gene_name"))

#### Generate the upset plots ####
generateUpsets <- function(specificContrast) {
  # To construct proper labels later, get the condition in which genes are up- and downregulated in this contrast:
  upregulatedExpression <- str_split_i(string = specificContrast,
                                       pattern = " vs. ",
                                       i = 2) %>%
    tolower()
  upregulatedExpression <- paste("Upregulated in",
                                 upregulatedExpression)
  
  downregulatedExpression <- str_split_i(string = specificContrast,
                                       pattern = " vs. ",
                                       i = 1) %>%
    tolower()
  downregulatedExpression <- paste("Upregulated in",
                                   downregulatedExpression)
  
  upsetInputData <- allResultsCombined %>%
    filter(contrast == specificContrast) %>%
   # filter(!is.na(kValue)) %>%
   # filter(!is.na(testPvalueFDR)) %>%
   # filter(!is.na(padj)) %>%
    mutate(`Upregulated` = case_when(log2FoldChange > 0 & 
                                          padj < 0.05 ~ TRUE,
                                        TRUE ~ FALSE)) %>%
    mutate(`Downregulated` = case_when(log2FoldChange < 0 & 
                                        padj < 0.05 ~ TRUE,
                                      TRUE ~ FALSE)) %>%
    mutate(`Positive selection` = case_when(testPvalueFDR < 0.05 &      
                                              backgroundPvalueFDR >= 0.05 &
                                              differencePvalueFDR < 0.05 ~ TRUE,
                                            TRUE ~ FALSE)) %>%
    mutate(`Intensified selection` = case_when(kValue > 1 &      
                                                 pValueFDR < 0.05 ~ TRUE,
                                               TRUE ~ FALSE)) %>%
    mutate(`Relaxed selection` = case_when(kValue < 1 &      
                                             pValueFDR < 0.05 ~ TRUE,
                                           TRUE ~ FALSE)) %>%
    dplyr::select(c("gene_name",
             "Upregulated",  
             "Downregulated",
             "Positive selection",
             "Intensified selection",
             "Relaxed selection")) %>%
    distinct()
  rownames(upsetInputData) <- upsetInputData$gene_name
  upsetInputData <- upsetInputData %>%
    dplyr::select(-c("gene_name"))
  
  categories <- c("Upregulated",  
                  "Downregulated",
                  "Positive selection",
                  "Intensified selection",
                  "Relaxed selection")
  
  theme_border <- theme_gray() + 
    theme(plot.background = element_rect(fill = NA, 
                                         colour = 'grey', 
                                         size = 1))
  
  hyphyAndDEUpset <- upset(upsetInputData, 
                           categories,
                           width_ratio = 0.2,
                           max_degree = 2,
                           mode = "inclusive_intersection",
                           intersections = list(c('Downregulated', 'Positive selection'),
                                                c('Upregulated', 'Positive selection'),
                                                c('Downregulated', 'Relaxed selection'),
                                                c('Upregulated', 'Relaxed selection'),
                                                c('Downregulated', 'Intensified selection'),
                                                c('Upregulated', 'Intensified selection')),
                           sort_intersections = FALSE,
                           sort_sets = FALSE,
                           queries = list(upset_query(intersect=c('Downregulated', 
                                                                  'Positive selection'),
                                                      color = '#bed0a5',
                                                      fill = '#bed0a5',
                                                      only_components = c('intersections_matrix', 'Intersection size')),
                                          upset_query(intersect=c('Upregulated', 
                                                                  'Positive selection'),
                                                      color = '#bed0a5',
                                                      fill = '#bed0a5',
                                                      only_components = c('intersections_matrix', 'Intersection size')),
                                          upset_query(intersect=c('Downregulated', 
                                                                  'Intensified selection'),
                                                      color = '#05a77d',
                                                      fill = '#05a77d',
                                                      only_components = c('intersections_matrix', 'Intersection size')),
                                          upset_query(intersect=c('Upregulated', 
                                                                  'Intensified selection'),
                                                      color = '#05a77d',
                                                      fill = '#05a77d',
                                                      only_components = c('intersections_matrix', 'Intersection size')),
                                          upset_query(intersect=c('Downregulated', 
                                                                  'Relaxed selection'),
                                                      color = '#c5e8a8',
                                                      fill = '#c5e8a8',
                                                      only_components = c('intersections_matrix', 'Intersection size')),
                                          upset_query(intersect=c('Upregulated', 
                                                                  'Relaxed selection'),
                                                      color = '#c5e8a8',
                                                      fill = '#c5e8a8',
                                                      only_components = c('intersections_matrix', 'Intersection size')),
                                          upset_query(set = 'Downregulated',
                                                      fill = '#b3b3b3'),
                                          upset_query(set = 'Upregulated',
                                                      fill = '#737272'),
                                          upset_query(set = 'Relaxed selection',
                                                      fill = '#c5e8a8'),
                                          upset_query(set = 'Intensified selection',
                                                      fill = '#05a77d'),
                                          upset_query(set = 'Positive selection',
                                                      fill = '#bed0a5')),
                           labeller = ggplot2::as_labeller(c('Downregulated' = downregulatedExpression,
                                                             'Upregulated' = upregulatedExpression,
                                                             "Relaxed selection" = "Relaxed selection",
                                                             "Intensified selection" = "Intensified selection",
                                                             "Positive selection" = "Positive selection")),
                           set_sizes = (upset_set_size() + 
                                          theme(axis.text.x = element_text(angle = 45)))) + 
    patchwork::plot_annotation(theme = theme_border)
  
  # Changing text colors doesn't work for one of the plots, but for future reference, the code is here:
  #,
  #base_annotations = list('Intersection size' = intersection_size(text_colors = c(on_background = 'black', 
                                                                                  #on_bar = 'black')))
  
  
  plot <- plot(hyphyAndDEUpset) +
    ggtitle(specificContrast) +
    theme(plot.title = element_text(size = 10, 
                                    hjust = 0.5)) 

  plot
  
  return(plot)
}
possiblyGenerateUpsets <- purrr::possibly(generateUpsets,
                                          otherwise = "Error")

contrasts <- unique(allResultsCombined$contrast) %>%
  na.omit() %>%
  as.character()

allPlots <- purrr::map(contrasts,
                       possiblyGenerateUpsets)


hyphyAndDEUpset <- cowplot::plot_grid(allPlots[[1]], allPlots[[2]], allPlots[[3]], allPlots[[4]],
                                      ncol = 2,
                                      scale = 0.9,
                                      labels = c('A.', 'B.', 'C.', 'D.'), 
                                      label_size = 14,
                                      label_x = 0.05,
                                      label_y = 0.95)

hyphyAndDEUpset

ggsave(filename = "./images/upsetPlot.png",
       plot = hyphyAndDEUpset,
       height = 10, 
       width = 16,
       units = "in",
       dpi = 600)

morphConstrasts <- cowplot::plot_grid(allPlots[[1]], allPlots[[2]],
                                      ncol = 2,
                                      scale = 0.9,
                                      labels = c('A.', 'B.'), 
                                      label_size = 14,
                                      label_x = 0.05,
                                      label_y = 0.95)

morphConstrasts

ggsave(filename = "./images/upsetPlotMorphConstrastsOnly.png",
       plot = morphConstrasts,
       height = 15, 
       width = 35,
       units = "cm",
       dpi = 1800)

#### Do hypergeometric tests: ####

hypergeometricTests <- function(specificContrast) {
  upsetInputData <- allResultsCombined %>%
    filter(contrast == specificContrast) %>%
    # filter(!is.na(kValue)) %>%
    # filter(!is.na(testPvalueFDR)) %>%
    # filter(!is.na(padj)) %>%
    mutate(`Upregulated` = case_when(log2FoldChange > 0 & 
                                       padj < 0.05 ~ TRUE,
                                     TRUE ~ FALSE)) %>%
    mutate(`Downregulated` = case_when(log2FoldChange < 0 & 
                                         padj < 0.05 ~ TRUE,
                                       TRUE ~ FALSE)) %>%
    mutate(`Positive selection` = case_when(testPvalueFDR < 0.05 &      
                                              backgroundPvalueFDR >= 0.05 &
                                              differencePvalueFDR < 0.05 ~ TRUE,
                                            TRUE ~ FALSE)) %>%
    mutate(`Intensified selection` = case_when(kValue > 1 &      
                                                 pValueFDR < 0.05 ~ TRUE,
                                               TRUE ~ FALSE)) %>%
    mutate(`Relaxed selection` = case_when(kValue < 1 &      
                                             pValueFDR < 0.05 ~ TRUE,
                                           TRUE ~ FALSE)) %>%
    select(c("gene_name",
             "Upregulated",  
             "Downregulated",
             "Positive selection",
             "Intensified selection",
             "Relaxed selection")) %>%
    distinct()
  rownames(upsetInputData) <- upsetInputData$gene_name
  upsetInputData <- upsetInputData %>%
    select(-c("gene_name"))
  
  #### Upregulated and positive ####
  totalUpAndPositive <- length(which(upsetInputData$Upregulated == TRUE &
                                       upsetInputData$`Positive selection` == TRUE))
  totalPositive <- length(which(upsetInputData$`Positive selection` == TRUE))
  totalCount <- length(unique(rownames(upsetInputData)))
  totalUpregulated <- length(which(upsetInputData$Upregulated == TRUE))
  
  pValuePositiveUp <- phyper(totalUpAndPositive-1,
                             totalPositive, 
                             totalCount - totalPositive, 
                             totalUpregulated, 
                             lower.tail = FALSE, 
                             log.p = FALSE) 
  
  #### Upregulated and relaxed ####
  totalUpAndRelaxed <- length(which(upsetInputData$Upregulated == TRUE &
                                      upsetInputData$`Relaxed selection` == TRUE))
  totalRelaxed <- length(which(upsetInputData$`Relaxed selection` == TRUE))
  totalCount <- length(unique(rownames(upsetInputData)))
  totalUpregulated <- length(which(upsetInputData$Upregulated == TRUE))
  
  pValueRelaxedUp <- phyper(totalUpAndRelaxed-1,
                            totalRelaxed, 
                            totalCount - totalRelaxed, 
                            totalUpregulated, 
                            lower.tail = FALSE, 
                            log.p = FALSE)
  #### Upregulated and intensified ####
  totalUpAndIntensified <- length(which(upsetInputData$Upregulated == TRUE &
                                          upsetInputData$`Intensified selection` == TRUE))
  totalIntensified <- length(which(upsetInputData$`Intensified selection` == TRUE))
  totalCount <- length(unique(rownames(upsetInputData)))
  totalUpregulated <- length(which(upsetInputData$Upregulated == TRUE))
  
  pValueIntensifiedUp <- phyper(totalUpAndIntensified-1,
                                totalIntensified, 
                                totalCount - totalIntensified, 
                                totalUpregulated, 
                                lower.tail = FALSE, 
                                log.p = FALSE)
  #### Downregulated and positive ####
  totalDownAndPositive <- length(which(upsetInputData$Downregulated == TRUE &
                                         upsetInputData$`Positive selection` == TRUE))
  totalPositive <- length(which(upsetInputData$`Positive selection` == TRUE))
  totalCount <- length(unique(rownames(upsetInputData)))
  totalDownregulated <- length(which(upsetInputData$Downregulated == TRUE))
  
  pValuePositiveDown <- phyper(totalDownAndPositive-1,
                               totalPositive, 
                               totalCount - totalPositive, 
                               totalDownregulated, 
                               lower.tail = FALSE, 
                               log.p = FALSE) 
  #### Downregulated and relaxed ####
  totalDownAndRelaxed <- length(which(upsetInputData$Downregulated == TRUE &
                                        upsetInputData$`Relaxed selection` == TRUE))
  totalRelaxed <- length(which(upsetInputData$`Relaxed selection` == TRUE))
  totalCount <- length(unique(rownames(upsetInputData)))
  totalDownregulated <- length(which(upsetInputData$Downregulated == TRUE))
  
  pValueRelaxedDown <- phyper(totalDownAndRelaxed-1,
                              totalRelaxed, 
                              totalCount - totalRelaxed, 
                              totalDownregulated, 
                              lower.tail = FALSE, 
                              log.p = FALSE)
  #### Downregulated and intensified ####
  totalDownAndIntensified <- length(which(upsetInputData$Downregulated == TRUE &
                                            upsetInputData$`Intensified selection` == TRUE))
  totalIntensified <- length(which(upsetInputData$`Intensified selection` == TRUE))
  totalCount <- length(unique(rownames(upsetInputData)))
  totalDownregulated <- length(which(upsetInputData$Downregulated == TRUE))
  
  pValueIntensifiedDown <- phyper(totalDownAndIntensified-1,
                                  totalIntensified, 
                                  totalCount - totalIntensified, 
                                  totalDownregulated, 
                                  lower.tail = FALSE, 
                                  log.p = FALSE)
  
  allpValues <- c(pValuePositiveUp,
                  pValueRelaxedUp,
                  pValueIntensifiedUp,
                  pValuePositiveDown,
                  pValueRelaxedDown,
                  pValueIntensifiedDown,
                  specificContrast)
  names(allpValues) <- c("pValuePositiveUp",
                         "pValueRelaxedUp",
                         "pValueIntensifiedUp",
                         "pValuePositiveDown",
                         "pValueRelaxedDown",
                         "pValueIntensifiedDown",
                         "specificContrast")
  return(allpValues)
}

contrasts <- unique(allResultsCombined$contrast) %>%
  na.omit() %>%
  as.character()

allPvalues <- purrr::map(contrasts,
                         hypergeometricTests)
allPvalues <- as.data.frame(do.call(rbind, allPvalues))

allPvalues <- allPvalues %>%
  pivot_longer(cols = -c(specificContrast),
               names_to = c("Overlap"))

# Do FDR corrections:
allPvalues <- allPvalues %>%
  mutate(pValueFDR = p.adjust(value, method='BH')) %>%
  select(-c(value)) %>%
  arrange(pValueFDR)

allPvalues$upregulated <- str_split_i(allPvalues$specificContrast,
                                      pattern = " vs. ",
                                      i = 2) %>%
  tolower()
allPvalues$downregulated <- str_split_i(allPvalues$specificContrast,
                                      pattern = " vs. ",
                                      i = 1) %>%
  tolower()

allPvalues <- allPvalues %>%
  mutate(expressionPattern = case_when(
    grepl("Down", Overlap) ~ paste("Upregulated in", downregulated, "relative to", upregulated),
    grepl("Up", Overlap) ~ paste("Upregulated in", upregulated, "relative to", downregulated)))
allPvalues <- allPvalues %>%
  mutate(selectionType = case_when(
    grepl("Relaxed", Overlap) ~ "Relaxed selection",
    grepl("Intensified", Overlap) ~ "Intensified selection",
    grepl("Positive", Overlap) ~ "Positive selection"))
allPvalues <- allPvalues %>% 
  select(c(expressionPattern,
           selectionType,
           pValueFDR))


# Export a table of results:
exportTable <- gt(allPvalues)
exportTable
gtsave(exportTable, 
       "overlapPValues.docx")


#### Bootstrapping? ####
#numberOfRelaxedDevelopmentalContrasts <- length(which((allPvalues$specificContrast == "Worker pupae vs. soldier pupae" | 
#                                                         allPvalues$specificContrast == "Adult workers vs. adult soldiers") &
#                                                        (allPvalues$Overlap == "pValueRelaxedUp" | 
#                                                           allPvalues$Overlap == "pValueRelaxedDown") &
#                                                        allPvalues$pValueFDR < 0.05))
#
#bootstrapRelaxDevelopment <- function(dataframe) {
#  bootstrap <- dataframe
#  bootstrap$pValueFDR <- sample(bootstrap$pValueFDR,
#                            length(bootstrap$pValueFDR))
#  bootstrapResult <- length(which((bootstrap$specificContrast == "Worker pupae vs. soldier pupae" | 
#                                     bootstrap$specificContrast == "Adult workers vs. adult soldiers") &
#                                    (bootstrap$Overlap == "pValueRelaxedUp" | 
#                                       bootstrap$Overlap == "pValueRelaxedDown") &
#                                    bootstrap$pValueFDR < 0.05))
#  return(bootstrapResult)
#}
#
#bootstraps <- replicate(1000000,
#                        bootstrapRelaxDevelopment(dataframe = allPvalues), 
#                        simplify = TRUE) %>%
#  as.vector()
#length(which(bootstraps > numberOfRelaxedDevelopmentalContrasts))/1000000
