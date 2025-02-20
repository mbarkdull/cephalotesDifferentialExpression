library(tidyverse)
library(scales)
library(googlesheets4)
library(splitstackshape)

#### Plot contrast-vs-contrast as suggested by Mike ####
# Prep the differential expression data:
differentialExpression <- read_csv(file = "./finalResults/allDifferentialExpressionResultsAndFunctions.csv") %>%
  dplyr::select(c(gene_name, log2FoldChange, padj, contrast)) %>%
  distinct() %>%
  filter(!is.na(contrast)) %>%
  pivot_wider(names_from = "contrast", values_from = c("log2FoldChange", "padj")) %>%
  mutate(significantAnywhere = case_when(`padj_Pupal workers vs. pupal soldiers` <= 0.05 |   
                                           `padj_Adult workers vs. adult soldiers` <= 0.05 |
                                           `padj_Adult soldiers vs. pupal soldiers` <= 0.05 |
                                           `padj_Adult workers vs. pupal workers`<= 0.05 ~ "Yes",
                                         TRUE ~ "No"),
         significantWithinLifeStage = case_when(`padj_Pupal workers vs. pupal soldiers` <= 0.05 |   
                                                  `padj_Adult workers vs. adult soldiers`<= 0.05 ~ "Yes",
                                                TRUE ~ "No"),
         significantWithinMorph = case_when(`padj_Adult workers vs. pupal workers` <= 0.05 |
                                              `padj_Adult soldiers vs. pupal soldiers` <= 0.05 ~ "Yes",
                                            TRUE ~ "No"))

differentialExpression$significantAnywhere <- factor(differentialExpression$significantAnywhere,
                                                     levels = c("Yes",
                                                                "No"))
differentialExpression$significantWithinLifeStage <- factor(differentialExpression$significantWithinLifeStage,
                                                            levels = c("Yes",
                                                                       "No"))
differentialExpression$significantWithinMorph <- factor(differentialExpression$significantWithinMorph,
                                                        levels = c("Yes",
                                                                   "No"))

# Get the list of Hippo genes:
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
                      cvarGene))
  return(expressionOfGene)
  
}
possiblyCheckHippoExpression <- possibly(checkHippoExpression, otherwise = "error")

hippoGenes <- purrr::map(hippoComponents,
                         possiblyCheckHippoExpression)
hippoGenes <- as.data.frame(do.call(rbind, hippoGenes))

# Add a column indicating Hippo genes in the differential expression data:
differentialExpression <- differentialExpression %>%
  dplyr::mutate(hippo = case_when(gene_name %in% hippoGenes$gene_name ~ "Yes",
                                  TRUE ~ "No"))

differentialExpression$hippo <- factor(differentialExpression$hippo,
                                       levels = c("Yes",
                                                  "No"))

# Make a color scheme
colorScheme <- c(`Yes` = "#1167a8D1",
                 `No` = "#b6bdc280")

# Contrasts within a single life stage, across castes:
lifeStageCorrelation <- cor.test(formula = ~ `log2FoldChange_Pupal workers vs. pupal soldiers` + `log2FoldChange_Adult workers vs. adult soldiers`,
                                 data = differentialExpression,
                                 method = "pearson")
lifeStageCorrelation

hippoColorScheme <- c(`Yes` = "red",
                      `No` = "#b6bdc200")

withinLifeStage <- ggplot(data = differentialExpression) +
  geom_point(mapping = aes(x = `log2FoldChange_Pupal workers vs. pupal soldiers`,
                           y = `log2FoldChange_Adult workers vs. adult soldiers`,
                           color = significantWithinLifeStage,
                           text = gene_name),
             size = 0.75) +
  geom_point(data = filter(differentialExpression,
                           significantWithinLifeStage == "Yes"),
             mapping = aes(x = `log2FoldChange_Pupal workers vs. pupal soldiers`,
                           y = `log2FoldChange_Adult workers vs. adult soldiers`,
                           color = significantWithinLifeStage),
             size = 0.75) +
  scale_color_manual(values = colorScheme) +
  labs(x = "Expression in worker pupae relative to soldier pupae",
       y = "Expression in worker adults relative to soldier adults") +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        legend.title = element_text(size = 8)) + 
  labs(fill = 'Signficant in at\nleast one of the\ntwo contrasts?',
       color = 'Gene in the\nHippo pathway') 

withinLifeStage

plotly::ggplotly(withinLifeStage)

# Contrasts within a single caste, over development:
casteCorrelation <- cor.test(formula = ~ `log2FoldChange_Adult soldiers vs. pupal soldiers` + `log2FoldChange_Adult workers vs. pupal workers`,
                             data = differentialExpression,
                             method = "pearson")
casteCorrelation

withinACaste <- ggplot(data = differentialExpression) +
  geom_point(mapping = aes(x = `log2FoldChange_Adult soldiers vs. pupal soldiers`,
                           y = `log2FoldChange_Adult workers vs. pupal workers`,
                           color = significantWithinMorph,
                           text = gene_name),
             size = 0.75) +
  geom_point(data = filter(differentialExpression,
                           significantWithinMorph == "Yes"),
             mapping = aes(x = `log2FoldChange_Adult soldiers vs. pupal soldiers`,
                           y = `log2FoldChange_Adult workers vs. pupal workers`,
                           color = significantWithinMorph),
             size = 0.75) +
  scale_color_manual(values = colorScheme) +
  labs(x = "Expression in soldier adults relative to soldier pupae",
       y = "Expression in worker adults relative to worker pupae") +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        legend.title = element_text(size = 8)) + 
  labs(fill = 'Signficant in at\nleast one of the\ntwo contrasts?') 

withinACaste

plotly::ggplotly(withinACaste)

# Arrange with patchwork:
library(patchwork)

# Get legends:
correlationLegend <- ggpubr::get_legend(withinLifeStage) %>%
  ggpubr::as_ggplot()

# Plot only the correlation plots:
(withinACaste &
    theme(legend.position = "none",
          plot.margin = unit(c(0, 30, 0, 0), "pt")) | withinLifeStage &
    theme(legend.position = "none",
          plot.margin = unit(c(0, 0, 0, 30), "pt")) |
    (correlationLegend)) + 
  patchwork::plot_layout(widths = c(1, 1, 0.45)) 

# Save a png:
ggsave("./images/logFoldScatter.png", 
       width = 24, 
       height = 15*0.66, 
       units = "cm", 
       dpi = 1200)

# Save an svg:
ggsave("./images/logFoldScatter.svg", 
       width = 30, 
       height = 20, 
       units = "cm", 
       dpi = 1200)

#### Generate volcano plots: ####
listOfContrasts <- colnames(differentialExpression) 
listOfContrasts <- listOfContrasts[grepl("vs", listOfContrasts)]
listOfContrasts <- str_replace(listOfContrasts,
                               pattern = ".*_",
                               replacement = "") %>%
  unique()

plotPairwiseContrasts <- function(contrast) {
  log2Column <- paste("log2FoldChange_",
                      contrast,
                      sep = "")
  padjColumn <- paste("padj_",
                      contrast,
                      sep = "")
  
  # Select only columns for the particular contrast:
  filteredDataframe <- differentialExpression %>%
    select(c(gene_name,
             log2Column,
             padjColumn,
             hippo))
  
  colnames(filteredDataframe) <- c("gene_name",
                                   "log2FoldChange",
                                   "padj",
                                   "hippo")
  
  # Filter the data to remove rows with no p-value:
  filteredDataframe <- filter(filteredDataframe,
                              !is.na(padj))
  
  # Add a column indicating the color for a gene, based on the adjusted p-value:
  filteredDataframe <- filteredDataframe %>% 
    mutate(color = case_when(padj <= 0.01 & log2FoldChange <= 0 ~ "p ≤ 0.01, upregulated in pupal workers",
                             padj <= 0.05 & log2FoldChange <= 0 ~ "p ≤ 0.05, upregulated in pupal workers",
                             padj > 0.05 & log2FoldChange <= 0 ~ "p > 0.05",
                             padj <= 0.01 & log2FoldChange > 0 ~ "p ≤ 0.01, upregulated in pupal soldiers",
                             padj <= 0.05 & log2FoldChange > 0 ~ "p ≤ 0.05, upregulated in pupal soldiers",
                             padj > 0.05 & log2FoldChange > 0 ~ "p > 0.05"))
  
  # Make a color scheme
  significanceColorScheme <- c("p > 0.05" = "#e8e6e6",
                               "p ≤ 0.05, upregulated in pupal workers" = "#a7cfc3",
                               "p ≤ 0.05, upregulated in pupal soldiers" = "#edce8c",
                               "p ≤ 0.01, upregulated in pupal workers" = "#00996e",
                               "p ≤ 0.01, upregulated in pupal soldiers" = "#f0a000")
  
  # Get the range limits for the x-axis, which should be the most extreme large or small value of the fold-change, so that the axis can be symmetrical around zero. 
  xLimRange <- c(abs(min(filteredDataframe$log2FoldChange,
                         na.rm = TRUE)), 
                 max(filteredDataframe$log2FoldChange,
                     na.rm = TRUE)) %>%
    max()
  
  plotTitle <- contrast
  plotTitle <- gsub(pattern = " vs.",
                    replacement = "\nvs.",
                    plotTitle) %>%
    tolower()
  
  ggplot() +
    geom_point(data = filteredDataframe,
               mapping = aes(x = log2FoldChange, 
                             y = padj,
                             color = color),
               size = 0.5, 
               alpha = 1) +
    scale_color_manual(values = significanceColorScheme) +
    
    geom_point(data = filter(filteredDataframe,
                             hippo == "Yes"),
               mapping = aes(x = log2FoldChange, 
                             y = padj,
                             fill = color),
               size = 1.5, 
               alpha = 1,
               pch = 21,
               color = "black",
               stroke = 1) +
    scale_fill_manual(values = significanceColorScheme,
                      guide="none") +
    
    scale_y_continuous(trans = compose_trans("log10", 
                                             "reverse"),
                       labels = label_log()) + 
    xlim(floor(-xLimRange),
         ceiling(xLimRange)) +
    labs(title = paste("Gene expression in",
                       tolower(contrast)), 
         x = "Log<sub>2</sub> fold change", 
         y = "-log<sub>10</sub>(adjusted p-value)",
         color = "Signficance of differential expression") +
    theme_bw() +
    theme(axis.title = ggtext::element_markdown(size = 8),
          plot.title = element_text(size = 9),
          legend.title = element_text(size = 8)) + 
    guides(colour = guide_legend(override.aes = list(size = 3)))
}

pupaPlot <- plotPairwiseContrasts(contrast = "Pupal workers vs. pupal soldiers")

pupaPlot

ggsave("./images/pupalVolcanoHippo.png",
       width = 7,
       height = 3, 
       units = "in",
       dpi = 1200)


