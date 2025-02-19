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
                           fill = significantWithinLifeStage,
                           color = hippo,
                           text = gene_name),
             size = 0.75,
             pch = 21) +
  geom_point(data = filter(differentialExpression,
                           significantWithinLifeStage == "Yes"),
             mapping = aes(x = `log2FoldChange_Pupal workers vs. pupal soldiers`,
                           y = `log2FoldChange_Adult workers vs. adult soldiers`,
                           fill = significantWithinLifeStage,
                           color = hippo),
             size = 0.75,
             pch = 21) +
  geom_point(data = filter(differentialExpression,
                           hippo == "Yes"),
             mapping = aes(x = `log2FoldChange_Pupal workers vs. pupal soldiers`,
                           y = `log2FoldChange_Adult workers vs. adult soldiers`,
                           fill = significantWithinLifeStage,
                           color = hippo),
             size = 2,
             pch = 21) +
  scale_fill_manual(values = colorScheme) +
  scale_color_manual(values = hippoColorScheme) +
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
                           fill = significantWithinMorph,
                           text = gene_name),
             size = 0.75,
             pch = 21,
             color = "#00000000") +
  geom_point(data = filter(differentialExpression,
                           significantWithinMorph == "Yes"),
             mapping = aes(x = `log2FoldChange_Adult soldiers vs. pupal soldiers`,
                           y = `log2FoldChange_Adult workers vs. pupal workers`,
                           fill = significantWithinMorph),
             size = 0.75,
             pch = 21,
             color = "#00000000") +
  scale_fill_manual(values = colorScheme) +
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
ggsave("./images/logFoldScatterHippo.png", 
       width = 24, 
       height = 15*0.66, 
       units = "cm", 
       dpi = 1200)

# Save an svg:
ggsave("./images/logFoldScatterHippo.svg", 
       width = 30, 
       height = 20, 
       units = "cm", 
       dpi = 1200)

