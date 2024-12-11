library(gganimate)
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(data.table)
library(splitstackshape)
library(snakecase)

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
  return(enrichment_kegg)
}

selectedContrast <- "Worker pupae vs. soldier pupae"

enrichment_kegg <- doKeggEnrichment(selectedContrast)

# Get a key that links signficant KEGG pathways to KEGG terms:
signficantKEGGPathToTerm <- enrichment_kegg@result %>%
  filter(p.adjust < 0.05) 
signficantKEGGPathToTerm <- concat.split.multiple(signficantKEGGPathToTerm, 
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
keggTermToGene <- concat.split.multiple(keggTermToGene, 
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


ridgelineData$`KEGG pathway` <- as.factor(ridgelineData$`KEGG pathway`)
ridgelineData <- ridgelineData %>%
  dplyr::arrange(desc(`KEGG pathway`))
ridgelineData <- na.omit(ridgelineData)

filteredRidgelineData <- ridgelineData %>%
  filter(!is.na(`KEGG pathway`),
         category != "Human Diseases") %>%
  mutate(`KEGG pathway` = fct_rev(as_factor(`KEGG pathway`))) 

pupaRidgeline <- ggplot(data = filteredRidgelineData,
                        mapping = aes(x = `Fold change`, 
                                      y = `KEGG pathway`,
                                      group = `KEGG pathway`,
                                      fill = `Adjusted p-value`)) +
  ggridges::geom_density_ridges() +
  scale_fill_gradientn(colors = c("#0455b3", 
                                  "#e8d056")) +
  theme_bw() + 
  transition_states(`KEGG pathway`,
                    transition_length = 1, 
                    state_length = 1,
                    wrap = FALSE) +
  enter_fade() +
  shadow_mark()

animate(pupaRidgeline, 
        height = 6,
        width = 10, 
        units = "in", 
        res = 400,
        rewind = FALSE,
        renderer = gifski_renderer(loop = FALSE))

anim_save("pupaRidgeline.gif")

staticPlot <- ggplot(data = filteredRidgelineData,
                     mapping = aes(x = `Fold change`, 
                                   y = `KEGG pathway`,
                                   group = `KEGG pathway`,
                                   fill = `Adjusted p-value`)) +
  ggridges::geom_density_ridges() +
  scale_fill_gradientn(colors = c("#0455b3", 
                                  "#e8d056")) +
  theme_bw()

ggsave(filename = "./images/staticPupalRidgeline.png", 
       plot = staticPlot,
       height = 6,
       width = 10, 
       units = "in", 
       dpi = 400)
