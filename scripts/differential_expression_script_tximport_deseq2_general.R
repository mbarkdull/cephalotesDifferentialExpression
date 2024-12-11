# User supplies:
  # path to Salmon outputs
  # link to google sheet with sample metadata
  # genome annotation 
  # study design for DESeq to set up contrasts. 

library(tidyverse)
library(scales)
library(DESeq2)
library(tximport)
library(googlesheets4)
library(rtracklayer)
library(splitstackshape)

#### Prep input data for DeSeq ####
# List all of the Salmon quantification output files:
salmon_outputs <- list.files(path = "./02_alignments/salmon/", 
                             full.names = TRUE,
                             recursive = TRUE,
                             pattern = "*sf")

# Read in the genome annotation to get a table linking transcripts to genes:
# Filter to only mRNA features:
annotation <- ape::read.gff(file = "./CVAR_OGS_v1.0.gff3") %>%
  filter(type == "mRNA")
# Split the attributes column to separate out ID and Parent:
annotation$ID <- stringr::str_extract(annotation$attributes,
                                      regex("ID[^;]+", 
                                            ignore_case = T)) %>%
  str_split_i(pattern = "=",
              i = 2)
annotation$Parent <- stringr::str_extract(annotation$attributes,
                                          regex("Parent[^;]+", 
                                                ignore_case = T)) %>%
  str_split_i(pattern = "=",
              i = 2)
# Create the dictionary linking transcripts to genes:
transcriptsToGenes <- dplyr::select(annotation, ID)
transcriptsToGenes$gene_name <- annotation$Parent

# Use tximport to import transcript-level abundance estimates across all samples, from the salmon outputs:
transcriptAbundances <- tximport(files = salmon_outputs, 
                                 type = "salmon", 
                                 tx2gene = transcriptsToGenes, 
                                 countsFromAbundance = "lengthScaledTPM")

# Make a table of sample metadata, used to set up the contrasts in the differential expression analysis:
sampleData <- read_sheet("https://docs.google.com/spreadsheets/d/1FEim9DwRUdOW5ZG0c6VzdkCdOIWDtCSekWrA_kN1AKk/edit?usp=sharing")
sampleData <- data.frame(sampleData, 
                         row.names = colnames(transcriptAbundances$counts))

#### Test for differential expression, using DESeq to run a Wald test: ####
# Combine all relevant information into a DESeq dataset object:
# Using a complex design formula:
dataForDESeq <- DESeqDataSetFromTximport(transcriptAbundances, 
                                         colData = sampleData, 
                                         design = ~ stage + caste + stage:caste + caste:stage)

# Using a simple design formula to do simple, pairwise comparisons:
dataForPairwise <- dataForDESeq
# Create a new condition that is the combination of caste and stage:
dataForPairwise$group <- factor(paste0(dataForPairwise$caste, 
                                       dataForPairwise$stage))
# Update the design to only consider the effects of this group variable:
design(dataForPairwise) <- ~ group
# Run the analysis:
pairwiseAnalysis <- DESeq(dataForPairwise)

# Use contrasts to pull out each interesting pairwise comparisons:
pupaOnlyResults <- results(pairwiseAnalysis, 
                           contrast = c("group", 
                                        "soldierpupa", 
                                        "workerpupa"))
pupaOnlyResultsDataframe <- data.frame(pupaOnlyResults)
pupaOnlyResultsDataframe$contrast <- "Pupal workers vs. pupal soldiers"

adultOnlyResults <- results(pairwiseAnalysis, 
                            contrast = c("group", 
                                         "soldieradult", 
                                         "workeradult"))
adultOnlyResultsDataframe <- data.frame(adultOnlyResults)
adultOnlyResultsDataframe$contrast <- "Adult workers vs. adult soldiers"

soldierOnlyResults <- results(pairwiseAnalysis, 
                              contrast = c("group", 
                                           "soldierpupa", 
                                           "soldieradult"))
soldierOnlyResultsDataframe <- data.frame(soldierOnlyResults)
soldierOnlyResultsDataframe$contrast <- "Adult soldiers vs. pupal soldiers"

workerOnlyResults <- results(pairwiseAnalysis, 
                             contrast = c("group", 
                                          "workerpupa", 
                                          "workeradult"))
workerOnlyResultsDataframe <- data.frame(workerOnlyResults)
workerOnlyResultsDataframe$contrast <- "Adult workers vs. pupal workers"

#### See, for each differentially expressed gene, if it is unique to that contrast or not: ####
pupalDEGenes <- filter(pupaOnlyResultsDataframe,
                       padj < 0.05)
adultDEGenes <- filter(adultOnlyResultsDataframe,
                       padj < 0.05)
workerDEGenes <- filter(workerOnlyResultsDataframe,
                        padj < 0.05)
soldierDEGenes <- filter(soldierOnlyResultsDataframe,
                         padj < 0.05)

# Pupal genes:
pupalDEGenes$notUnique <- row.names(pupalDEGenes) %in%
  c(row.names(adultDEGenes),
    row.names(workerDEGenes),
    row.names(soldierDEGenes))
length(which(pupalDEGenes$notUnique == FALSE))
uniquePupalDEGenes <- filter(pupalDEGenes,
                             notUnique == FALSE)
length(which(uniquePupalDEGenes$log2FoldChange < 0))
length(which(uniquePupalDEGenes$log2FoldChange > 0))

# Adult genes:
adultDEGenes$notUnique <- row.names(adultDEGenes) %in%
  c(row.names(pupalDEGenes),
    row.names(workerDEGenes),
    row.names(soldierDEGenes))
length(which(adultDEGenes$notUnique == FALSE))
uniqueAdultDEGenes <- filter(adultDEGenes,
                             notUnique == FALSE)
length(which(uniqueAdultDEGenes$log2FoldChange < 0))
length(which(uniqueAdultDEGenes$log2FoldChange > 0))

# Worker genes:
workerDEGenes$notUnique <- row.names(workerDEGenes) %in%
  c(row.names(pupalDEGenes),
    row.names(adultDEGenes),
    row.names(soldierDEGenes))
length(which(workerDEGenes$notUnique == FALSE))
uniqueWorkerDEGenes <- filter(workerDEGenes,
                             notUnique == FALSE)
length(which(uniqueWorkerDEGenes$log2FoldChange < 0))
length(which(uniqueWorkerDEGenes$log2FoldChange > 0))

# Soldier genes:
soldierDEGenes$notUnique <- row.names(soldierDEGenes) %in%
  c(row.names(pupalDEGenes),
    row.names(adultDEGenes),
    row.names(workerDEGenes))
length(which(soldierDEGenes$notUnique == FALSE))
uniqueSoldierDEGenes <- filter(soldierDEGenes,
                             notUnique == FALSE)
length(which(uniqueSoldierDEGenes$log2FoldChange < 0))
length(which(uniqueSoldierDEGenes$log2FoldChange > 0))

# Worker genes not differentially expressed in soldiers:
workerDEGenes$alsoInSoldier <- row.names(workerDEGenes) %in%
  c(row.names(soldierDEGenes))
length(which(workerDEGenes$alsoInSoldier == FALSE))
workerDEGenesNotSolider <- filter(workerDEGenes,
                                  alsoInSoldier == FALSE)
length(which(workerDEGenesNotSolider$log2FoldChange < 0))
length(which(workerDEGenesNotSolider$log2FoldChange > 0))

# Soldier genes not differentially expressed in workers:
soldierDEGenes$alsoInWorker <- row.names(soldierDEGenes) %in%
  c(row.names(workerDEGenes))
length(which(soldierDEGenes$alsoInWorker == FALSE))
soldierDEGenesNotWorker <- filter(soldierDEGenes,
                                  alsoInWorker == FALSE)
length(which(soldierDEGenesNotWorker$log2FoldChange < 0))
length(which(soldierDEGenesNotWorker$log2FoldChange > 0))

#### Generate volcano plots for the pairwise comparisons: ####
# Write a function to do a single volcano plot:
plotPairwiseContrasts <- function(dataframe) {
  # Filter the data to remove rows with no p-value:
  filteredDataframe <- filter(dataframe,
                              !is.na(padj))
  
  # Add a column indicating the color for a gene, based on the adjusted p-value:
  filteredDataframe <- filteredDataframe %>% 
    mutate(color = case_when(padj <= 0.01 ~ "p <= 0.01",
                             padj <= 0.05  ~ "p <= 0.05",
                             padj > 0.05 ~ "p > 0.05"))
  
  # Make a color scheme
  colorScheme <- c(`p <= 0.01` = "#06A77D",
                   `p <= 0.05` = "#1378ed",
                   `p > 0.05` = "#D5C67A")
  
  # Get the range limits for the x-axis, which should be the most extreme large or small value of the fold-change, so that the axis can be symmetrical around zero. 
  xLimRange <- c(abs(min(filteredDataframe$log2FoldChange,
                         na.rm = TRUE)), 
                 max(filteredDataframe$log2FoldChange,
                     na.rm = TRUE)) %>%
    max()
  ggplot(data = filteredDataframe) +
    geom_point(mapping = aes(x = log2FoldChange, 
                             y = padj,
                             color = color),
               size = 0.25, 
               alpha = 0.5) +
    scale_y_continuous(trans = compose_trans("log10", 
                                             "reverse"),
                       labels = label_log()) + 
    xlim(-xLimRange - 5,
         xLimRange + 5) +
    #geom_hline(yintercept = 0.05, 
    #           linewidth = 1, 
    #           colour = "#FF3721", 
    #           linetype = "dashed") + 
    #geom_hline(yintercept = 0.01, 
    #           linewidth = 1, 
    #           colour = "#FFA500", 
    #           linetype = "dashed") + 
    scale_color_manual(values = colorScheme) +
    labs(title = unique(dataframe$contrast), 
         x = "Log<sub>2</sub> fold change", 
         y = "-log<sub>10</sub>(adjusted p-value)",
         color = "Signficance of\ndifferential expression") +
    theme_bw() +
    theme(axis.title = ggtext::element_markdown(size = 8),
          plot.title = element_text(size = 9),
          #axis.text.y = element_text(angle = 90),   #### Skipping because I can't figure out how to change the number of breaks
          legend.title = element_text(size = 8))
}

# List all of the pairwise comparisons:
resultsList <- list(pupaOnlyResultsDataframe,
                    adultOnlyResultsDataframe,
                    soldierOnlyResultsDataframe,
                    workerOnlyResultsDataframe)

# Use purrr to generate a volcano plot for each comparison:
volcanoPlotList <- purrr::map(resultsList, 
                              plotPairwiseContrasts)

# Write a function to make scales and axes consistent across plots: 
makePlotsConsistent <- function(plot){
  # Get the total number of plots that are combined together:
  num_plots <- length(plot)
  
  # Fix x limits:
  xLimits <- lapply(1:num_plots, function(x) ggplot_build(plot[[x]])$layout$panel_scales_x[[1]]$range$range)
  # Get the minimum and maximum x values for each plot:
  minX <- min(unlist(xLimits))
  maxX <- max(unlist(xLimits))
  xLimRange <- c(abs(minX), 
                 maxX) %>%
    max()
  
  # Fix y limits:
  yLimits <- lapply(1:num_plots, function(x) ggplot_build(plot[[x]])$layout$panel_scales_y[[1]]$range$range)
  # Get the minimum and maximum x values for each plot:
  minY <- min(unlist(yLimits))
  maxY <- max(unlist(yLimits))
  
  plot & 
    xlim(-xLimRange - 0.01, 
         xLimRange + 0.01) &
    scale_y_continuous(trans = compose_trans("log10", 
                                             "reverse"),
                       labels = label_log())
}

# Plot them all together using patchwork:
allVolcanoes <- patchwork::wrap_plots(volcanoPlotList, 
                                      ncol = 2) + 
  patchwork::plot_annotation(title = 'Differential gene expression across caste and life stage contrasts',
                             theme = theme(plot.title = element_text(size = 19)))

allVolcanoes <- makePlotsConsistent(allVolcanoes) + 
  patchwork::plot_layout(guides = "collect",
                         axes = "collect")

allVolcanoes

ggsave(filename = "./images/volcanoPlots.png",
       width = 12,
       height = 6,
       units = "in", 
       dpi = 900)

#### Link the cvar genes to functions from other species based on OrthoFinder clustering ####
# Read in the orthogroup membership file:
orthogroupMembers <- read_delim(file = "03_Orthofinder/fasta/OrthoFinder/Results_Aug16/Orthogroups/Orthogroups.tsv",
                                delim = "\t")
# Split the column that lists genes from Cephalotes varians in a long direction (so each Cephalotes varians gene has one row linking it to other genes):
cvarToFunctions <- splitstackshape::concat.split.multiple(orthogroupMembers, 
                                         split.cols = "cvar_proteins", 
                                         seps = ",", 
                                         direction = "long") 
# Filter any orthogroups that don't have a Cephalotes varians protein, since these are not useful:
cvarToFunctions <- filter(cvarToFunctions, 
                          !is.na(cvar_proteins))

# Pivot this long, so there is only one column with genes from other species:
cvarToFunctionsLong <- cvarToFunctions %>% 
  pivot_longer(cols = c("acep_proteins",
                        "acol_proteins",
                        "amel_proteins",
                        "dmel_proteins",
                        "tcas_proteins"),
               names_to='species',
               values_to='gene') %>% 
  filter(!is.na(gene))
cvarToFunctionsLong <- splitstackshape::concat.split.multiple(cvarToFunctionsLong, 
                                             split.cols = "gene", 
                                             seps = ",", 
                                             direction = "long") 

# Write a function to read in the protein fasta files, since these actually have the gene functions:
readFastaGetFunction <- function(file) {
  species <- str_split_i(string = file, 
                         pattern = "/",
                         4)
  species <- str_split_i(string = species, 
                         pattern = "\\.",
                         1)
  amelProteins <- phylotools::read.fasta(file = file) %>%
    dplyr::select(seq.name)
  # Get just the locus name to join up with the big dictionary we're building:
  amelProteins$geneName <- str_split_i(string = amelProteins$seq.name, 
                                       pattern = " ", 
                                       i = 1)
  amelProteins$species <- species
  return(amelProteins)
}

# Iterate this over all proteomes:
allProteomes <- list.files(path = "./03_Orthofinder/fasta",
                           pattern = "*fasta", 
                           full.names = TRUE)
possiblyReadFastaGetFunction <- purrr::possibly(readFastaGetFunction, 
                                                otherwise = "error")
allProteomes <- purrr::map(allProteomes, 
                           possiblyReadFastaGetFunction)
allProteomes <- as.data.frame(do.call(rbind, allProteomes))

# Combine the other-species-genes-to-functions table with the cvar-genes-to-other-species-genes table:
allFunctionInformation <- full_join(cvarToFunctionsLong,
                                    allProteomes, 
                                    by = c("species" = "species",
                                           "gene" = "geneName")) %>%
  distinct()

# Combine this table, which is transcript-based, with the gene-to-transcript key:
functionsDictionary <- full_join(transcriptsToGenes,
                                 allFunctionInformation,
                                 by = c("ID" = "cvar_proteins"))

# Combine the functions with the differential expression results: 
combineFunctionsWithResults <- function(resultsDataframe) {
  modifiedResults <- resultsDataframe
  modifiedResults$gene_name <- row.names(modifiedResults)
  modifiedResults <- full_join(modifiedResults, 
                               functionsDictionary,
                               by = c("gene_name" = "gene_name")) %>%
    dplyr::select(c(gene_name,
                    seq.name,
                    baseMean,
                    log2FoldChange,
                    lfcSE,
                    stat,
                    pvalue,
                    padj, 
                    contrast)) %>%
    distinct()
  return(modifiedResults)
}

pupaResultsFunctions <- combineFunctionsWithResults(pupaOnlyResultsDataframe)
pupaResultsFunctions$uniqueToContrast <- pupaResultsFunctions$gene_name %in% 
  row.names(uniquePupalDEGenes)

adultResultsFunctions <- combineFunctionsWithResults(adultOnlyResultsDataframe)
adultResultsFunctions$uniqueToContrast <- adultResultsFunctions$gene_name %in% 
  row.names(uniqueAdultDEGenes)

soldierResultsFunctions <- combineFunctionsWithResults(soldierOnlyResultsDataframe)
soldierResultsFunctions$uniqueToContrast <- soldierResultsFunctions$gene_name %in% 
  row.names(uniqueSoldierDEGenes)

workerResultsFunctions <- combineFunctionsWithResults(workerOnlyResultsDataframe)
workerResultsFunctions$uniqueToContrast <- workerResultsFunctions$gene_name %in% 
  row.names(uniqueWorkerDEGenes)

# Combine all of these dataframes and save them as a csv file:
allResultsForExport <- rbind(pupaResultsFunctions,
                             adultResultsFunctions,
                             soldierResultsFunctions,
                             workerResultsFunctions)

# Export the results as a .csv:
write_csv(allResultsForExport,
          file = "./finalResults/allDifferentialExpressionResultsAndFunctions.csv")

#### Assess how many genes are differentially expressed across all contrasts: ####
allResults <- read_csv("./finalResults/allDifferentialExpressionResultsAndFunctions.csv") %>%
  select(-c("seq.name")) %>%
  distinct()

length(which(allResults$contrast == "Pupal workers vs. pupal soldiers" &
               allResults$log2FoldChange < 0 &
               allResults$padj < 0.05))
length(which(allResults$contrast == "Pupal workers vs. pupal soldiers" &
               allResults$log2FoldChange > 0 &
               allResults$padj < 0.05))

length(which(allResults$contrast == "Adult workers vs. adult soldiers" &
               allResults$log2FoldChange < 0 &
               allResults$padj < 0.05))
length(which(allResults$contrast == "Adult workers vs. adult soldiers" &
               allResults$log2FoldChange > 0 &
               allResults$padj < 0.05))

length(which(allResults$contrast == "Adult soldiers vs. pupal soldiers" &
               allResults$log2FoldChange < 0 &
               allResults$padj < 0.05))
length(which(allResults$contrast == "Adult soldiers vs. pupal soldiers" &
               allResults$log2FoldChange > 0 &
               allResults$padj < 0.05))

length(which(allResults$contrast == "Adult workers vs. pupal workers" &
               allResults$log2FoldChange < 0 &
               allResults$padj < 0.05))
length(which(allResults$contrast == "Adult workers vs. pupal workers" &
               allResults$log2FoldChange > 0 &
               allResults$padj < 0.05))


#### Assess how many of the C. varians genes have matches in another species: ####
genesToOrthologs <- dplyr::select(allResultsForExport, 
                                  c(gene_name,
                                    seq.name)) %>%
  filter(!is.na(gene_name)) %>%
  distinct() %>% 
  group_by(gene_name) %>%
  arrange(seq.name) %>%
  slice_head(n = 5) %>%
  mutate(index = row_number()) %>%
  pivot_wider(names_from = index, 
              values_from = seq.name, 
              names_prefix = 'Names_')

totalGenes <- length(genesToOrthologs$gene_name)
genesWithOrthologs <- length(which(!is.na(genesToOrthologs$Names_1)))
percentGenesWithOrthologs <- genesWithOrthologs/totalGenes











# Run the differential expression analysis:
differentialExpression <- DESeq(dataForDESeq)

#### Generate dot plots for individual genes as a sanity check ####
plotCounts(pairwiseAnalysis,
           gene = "CVAR_10161",
           intgroup = c("group"))

#### Run a PCA of the samples ####
transformedData <- rlog(dataForPairwise, 
                        blind = TRUE)
DESeq2::plotPCA(transformedData,
                intgroup = c("group"),
                ntop = 10000)

#### Get normalized read counts across all samples ####
counts <- estimateSizeFactors(pairwiseAnalysis) 
counts <- counts(counts, 
                 normalized = TRUE) %>%
  as.data.frame()
colnames(counts) <- pairwiseAnalysis@colData$id
counts$gene <- rownames(counts)

averageCounts <- counts %>%
  rowwise() %>%
  mutate(adultSoldierMean = mean(c_across(c("MEB_AS_1",
                                            "MEB_AS_2",
                                            "MEB_AS_3")),
                                 na.rm = TRUE),
         adultWorkerMean = mean(c_across(c("MEB_AW_1",
                                           "MEB_AW_2",
                                           "MEB_AW_3")),
                                na.rm = TRUE),
         pupalSoldierMean = mean(c_across(c("MEB_SP_1",
                                            "MEB_SP_2",
                                            "MEB_SP_3")),
                                 na.rm = TRUE),
         pupalWorkerMean = mean(c_across(c("MEB_WP_3",
                                           "MEB_WP_4",
                                           "MEB_WP_5")),
                                na.rm = TRUE)) %>%
  select(c(gene,
           adultSoldierMean,
           adultWorkerMean,
           pupalSoldierMean,
           pupalWorkerMean))

totalGenes <- length(averageCounts$gene)

adultWorkerSpecific <- filter(averageCounts, 
                              adultSoldierMean == 0 &
                                adultWorkerMean > 0)
adultWorkerSpecificCount <- length(adultWorkerSpecific$gene)
adultWorkerSpecificCount/totalGenes
  
adultSoldierSpecific <- filter(averageCounts, 
                               adultWorkerMean == 0 &
                                adultSoldierMean > 0)
adultSoldierSpecificCount <- length(adultSoldierSpecific$gene)
adultSoldierSpecificCount/totalGenes

pupalWorkerSpecific <- filter(averageCounts, 
                              pupalSoldierMean == 0 &
                                pupalWorkerMean > 0)
pupalWorkerSpecificCount <- length(pupalWorkerSpecific$gene)
pupalWorkerSpecificCount/totalGenes

pupalSoldierSpecific <- filter(averageCounts, 
                               pupalWorkerMean == 0 &
                                pupalSoldierMean > 0)
pupalSoldierSpecificCount <- length(pupalSoldierSpecific$gene)
pupalSoldierSpecificCount/totalGenes

#### Plot contrast-vs-contrast as suggested by Mike ####
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

# Make a color scheme
colorScheme <- c(`Yes` = "#0B95A7",
                 `No` = "#D5C67A")

# Contrasts within a single life stage, across castes:
lifeStageCorrelation <- cor.test(formula = ~ `log2FoldChange_Pupal workers vs. pupal soldiers` + `log2FoldChange_Adult workers vs. adult soldiers`,
                                 data = differentialExpression,
                                 method = "pearson")
lifeStageCorrelation

withinLifeStage <- ggplot(data = differentialExpression) +
  geom_point(mapping = aes(x = `log2FoldChange_Pupal workers vs. pupal soldiers`,
                           y = `log2FoldChange_Adult workers vs. adult soldiers`,
                           color = significantWithinLifeStage,
                           text = gene_name),
             size = 0.75,
             alpha = 0.5) +
  geom_point(data = filter(differentialExpression,
                           significantWithinLifeStage == "Yes"),
             mapping = aes(x = `log2FoldChange_Pupal workers vs. pupal soldiers`,
                           y = `log2FoldChange_Adult workers vs. adult soldiers`,
                           color = significantWithinLifeStage),
             size = 0.75,
             alpha = 0.75) +
  scale_color_manual(values = colorScheme) +
  labs(x = "Expression in worker pupae relative to soldier pupae",
       y = "Expression in worker adults relative to soldier adults") +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        legend.title = element_text(size = 8)) + 
  labs(color = 'Signficant in at\nleast one of the two contrasts') 

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
             size = 0.75,
             alpha = 0.5) +
  geom_point(data = filter(differentialExpression,
                           significantWithinMorph == "Yes"),
             mapping = aes(x = `log2FoldChange_Adult soldiers vs. pupal soldiers`,
                           y = `log2FoldChange_Adult workers vs. pupal workers`,
                           color = significantWithinMorph),
             size = 0.75,
             alpha = 0.75) +
  scale_color_manual(values = colorScheme) +
  labs(x = "Expression in soldier adults relative to soldier pupae",
       y = "Expression in worker adults relative to worker pupae") +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        legend.title = element_text(size = 8)) + 
  labs(color = 'Signficant in at\nleast one of the two contrasts') 

withinACaste

plotly::ggplotly(withinACaste)

# Arrange with patchwork:
library(patchwork)

# Combo method one:
castePlots <- withinACaste / (allVolcanoes[[4]] + allVolcanoes[[3]]) + 
  patchwork::plot_layout(guides = "collect",
                         heights = c(2, 1)) +
  patchwork::plot_annotation("Gene expression within morphs across life stages")

lifestagePlots <- withinLifeStage / (allVolcanoes[[2]] + allVolcanoes[[1]]) + 
  patchwork::plot_layout(guides = "collect",
                         heights = c(2, 1)) &
  patchwork::plot_annotation("Gene expression within life stages across morphs")

# Get legends:
correlationLegend <- ggpubr::get_legend(withinACaste) %>%
  ggpubr::as_ggplot()
volcanoLegend <- ggpubr::get_legend(allVolcanoes[[1]]) %>%
  ggpubr::as_ggplot()

(castePlots &
  theme(legend.position = "none") | lifestagePlots &
  theme(legend.position = "none") |
  (correlationLegend / volcanoLegend)) + 
  patchwork::plot_layout(widths = c(1, 1, 0.45)) 

ggsave("./images/logFoldScatterPlusVolcanoes.png", 
       width = 30, 
       height = 20, 
       units = "cm", 
       dpi = 1200)

# Combo method two:
scatterPlots <- (withinACaste | withinLifeStage) + 
  patchwork::plot_layout(guides = "collect",
                         axes = "collect")
volcanoPlots <- (allVolcanoes[[3]] | allVolcanoes[[4]] | allVolcanoes[[1]] | allVolcanoes[[2]]) + 
  patchwork::plot_layout(guides = "collect",
                         axes = "collect") 

(scatterPlots / volcanoPlots) +
  patchwork::plot_layout(heights = c(2, 1))

# Combo method three, arrange with cowplot:
library(cowplot)
casteVolcanoes <- plot_grid(allVolcanoes[[3]] + theme(legend.position="none"), 
                            allVolcanoes[[4]] + theme(legend.position="none"), 
                            nrow = 1)
withinCastePlots <- plot_grid(withinACaste + theme(legend.position="none"), 
                              casteVolcanoes,
                              ncol = 1) 

lifestageVolcanoes <- plot_grid(allVolcanoes[[1]] + theme(legend.position="none"), 
                                allVolcanoes[[2]] + theme(legend.position="none"), 
                                nrow = 1)
withinLifestagePlots <- plot_grid(withinLifeStage + theme(legend.position="none"), 
                                  lifestageVolcanoes,
                                  ncol = 1) 

plot_grid(withinCastePlots,
          withinLifestagePlots,
          ncol = 2) 



