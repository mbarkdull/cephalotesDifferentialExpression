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

# Make an object that connects transcript names of the salmon files to gene names from the genome annotation:
# Read in the genome annotation to get a table linking transcripts to genes:
annotation <- ape::read.gff(file = "./CVAR_OGS_v1.0.gff3") %>%
  filter(type == "mRNA")
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

transcriptsToGenes <- select(annotation, ID)
transcriptsToGenes$gene_name <- annotation$Parent


# Use tximport to import transcript-level abundance estimates across all samples:
transcriptAbundances <- tximport(files = salmon_outputs, 
                                 type = "salmon", 
                                 tx2gene = transcriptsToGenes, 
                                 countsFromAbundance = "lengthScaledTPM")

# Make a table of sample metadata, used to set up the contrasts in the differential expression analysis:
sampleData <- read_sheet("https://docs.google.com/spreadsheets/d/1FEim9DwRUdOW5ZG0c6VzdkCdOIWDtCSekWrA_kN1AKk/edit?usp=sharing")
sampleData <- data.frame(sampleData, 
                         row.names = colnames(transcriptAbundances$counts))

#### Test for differential expression: ####
# Use DESeq to run a Wald test for differential expression:
# Combine all relevant information into a DESeq dataset object:
dataForDESeq <- DESeqDataSetFromTximport(transcriptAbundances, 
                                         colData = sampleData, 
                                         design = ~ stage + caste + stage:caste + caste:stage)

#### Do simple, pairwise comparisons ####
dataForPairwise <- dataForDESeq
# Create a new condition that is the combination of caste and stage:
dataForPairwise$group <- factor(paste0(dataForPairwise$caste, 
                                       dataForPairwise$stage))
# Update the design to only consider the effects of this group variable:
design(dataForPairwise) <- ~ group
# Run the analysis:
pairwiseAnalysis <- DESeq(dataForPairwise)

# Use contrasts to pull out single pairwise comparisons:
pupaOnlyResults <- results(pairwiseAnalysis, 
                           contrast = c("group", 
                                        "soldierpupa", 
                                        "workerpupa"))
pupaOnlyResultsDataframe <- data.frame(pupaOnlyResults)
pupaOnlyResultsDataframe$contrast <- "Worker pupae vs. soldier pupae"

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
soldierOnlyResultsDataframe$contrast <- "Adult soldiers vs. soldier pupae"

workerOnlyResults <- results(pairwiseAnalysis, 
                             contrast = c("group", 
                                          "workerpupa", 
                                          "workeradult"))
workerOnlyResultsDataframe <- data.frame(workerOnlyResults)
workerOnlyResultsDataframe$contrast <- "Adult workers vs. worker pupae"

# Do volcano plots for pairwise comparisons:
plotPairwiseContrasts <- function(dataframe) {
  # Filter the data to remove rows with no p-value:
  filteredDataframe <- filter(dataframe,
                              !is.na(padj))
  # Get the range limits for the x-axis, which should be the most extreme large or small value of the fold-change, so that the axis can be symetrical around zero. 
  xLimRange <- c(abs(min(filteredDataframe$log2FoldChange,
                         na.rm = TRUE)), 
                 max(filteredDataframe$log2FoldChange,
                     na.rm = TRUE)) %>%
    max()
  ggplot(data = filteredDataframe) +
    geom_point(mapping = aes(x = log2FoldChange, 
                             y = padj),
               size = 0.5, 
               alpha = 0.5) +
    scale_y_continuous(trans = compose_trans("log10", 
                                             "reverse"),
                       labels = label_log()) + 
    xlim(-xLimRange - 5,
         xLimRange + 5) +
    geom_hline(yintercept = 0.05, 
               linewidth = 1, 
               colour = "#FF3721", 
               linetype = "dashed") + 
    geom_hline(yintercept = 0.01, 
               linewidth = 1, 
               colour = "#FFA500", 
               linetype = "dashed") + 
    labs(title = unique(dataframe$contrast), 
         x = paste("Fold expression change between",
                   tolower(unique(dataframe$contrast))), 
         y = "Negative log10-adjusted p-value") +
    theme_bw()
}

resultsList <- list(pupaOnlyResultsDataframe,
                    adultOnlyResultsDataframe,
                    soldierOnlyResultsDataframe,
                    workerOnlyResultsDataframe)
volcanoPlotList <- purrr::map(resultsList, 
                              plotPairwiseContrasts)
patchwork::wrap_plots(volcanoPlotList, 
                      ncol = 2) + 
  patchwork::plot_annotation(title = 'Differential gene expression across caste and life stage contrasts',
                             theme = theme(plot.title = element_text(size = 19)))



#### Create a dictionary linking cvar genes to functions from other species from OrthoFinder ####
orthogroupMembers <- read_delim(file = "03_Orthofinder/fasta/OrthoFinder/Results_Apr30/Orthogroups/Orthogroups.tsv",
                                delim = "\t")
cvarToFunctions <- concat.split.multiple(orthogroupMembers, 
                                         split.cols = "cvar_proteins", 
                                         seps = ",", 
                                         direction = "long") 
cvarToFunctions <- filter(cvarToFunctions, 
                          !is.na(cvar_proteins))
cvarToFunctionsLong <- cvarToFunctions %>% 
  pivot_longer(cols = c("acep_proteins",
                        "acol_proteins",
                        "amel_proteins",
                        "dmel_proteins",
                        "tcas_proteins"),
               names_to='species',
               values_to='gene') %>% 
  filter(!is.na(gene))
cvarToFunctionsLong <- concat.split.multiple(cvarToFunctionsLong, 
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
    select(seq.name)
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

# Combine the functions with the cvar-to-other-species table:
allFunctionInformation <- full_join(cvarToFunctionsLong,
                                    allProteomes, 
                                    by = c("species" = "species",
                                           "gene" = "geneName")) %>%
  distinct()

# Combine this table, which is transcript-based, with the gene-to-transcript key:
functionsDictionary <- full_join(transcriptsToGenes,
                                 allFunctionInformation,
                                 by = c("ID" = "cvar_proteins"))
test <- select(functionsDictionary, -c("Orthogroup"))
test <- pivot_wider(test,
                    names_from = species, 
                    values_from = c(gene, seq.name)) %>%
  dplyr::select(-c("gene_NA",
                   "gene_cvar_proteins",
                   "seq.name_NA",
                   "seq.name_cvar_proteins")) 

##### Combine the functions with the differential expression results: ####
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
adultResultsFunctions <- combineFunctionsWithResults(adultOnlyResultsDataframe)
soldierResultsFunctions <- combineFunctionsWithResults(soldierOnlyResultsDataframe)
workerResultsFunctions <- combineFunctionsWithResults(workerOnlyResultsDataframe)

# Combine all of these dataframes and save them as a csv file:
allResultsForExport <- rbind(pupaResultsFunctions,
                             adultResultsFunctions,
                             soldierResultsFunctions,
                             workerResultsFunctions)

write_csv(allResultsForExport,
          file = "allDifferentialExpressionResultsAndFunctions.csv")




















# Run the differential expression analysis:
differentialExpression <- DESeq(dataForDESeq)

#### Explore results: ####
# List all contrasts:
resultsNames(differentialExpression)
# Extract the results for particular contrast of interest:
resultsCaste <- results(differentialExpression, 
                        name = "caste_worker_vs_soldier" )
mcols(resultsCaste)$description
resultsCasteDataframe <- data.frame(resultsCaste)
significantCaste <- filter(resultsCasteDataframe, 
                           padj <= 0.05)

# Volcano plots with ggplot:
ggplot(data = resultsCasteDataframe) +
  geom_point(mapping = aes(x = log2FoldChange, 
                           y = padj),
             size = 0.5) +
  scale_y_continuous(trans = compose_trans("log10", 
                                           "reverse"),
                     labels = label_log()) + 
  geom_hline(yintercept = 0.05, 
             size = 1, 
             colour = "#FF3721", 
             linetype = "dashed") + 
  geom_hline(yintercept = 0.01, 
             size = 1, 
             colour = "#FFA500", 
             linetype = "dashed") + 
  labs(title = "Genes differentially expressed between workers and soldiers", 
       x = "Fold expression change between workers and soldiers (positive means higher in workers)", 
       y = "Negative log10-adjusted p-value") 






# Extract the results for particular contrast of interest:
resultsStage <- results(differentialExpression, 
                        name = "stage_pupa_vs_adult" )
resultsStageDataframe <- data.frame(resultsStage)
significantStage <- filter(resultsStageDataframe, 
                           padj <= 0.05)

# Volcano plots with ggplot:
ggplot(data = resultsStageDataframe) +
  geom_point(mapping = aes(x = log2FoldChange, 
                           y = padj),
             size = 0.5) +
  scale_y_continuous(trans = compose_trans("log10", 
                                           "reverse"),
                     labels = label_log()) + 
  geom_hline(yintercept = 0.05, 
             size = 1, 
             colour = "#FF3721", 
             linetype = "dashed") + 
  geom_hline(yintercept = 0.01, 
             size = 1, 
             colour = "#FFA500", 
             linetype = "dashed") + 
  labs(title = "Genes differentially expressed between pupae and adults", 
       x = "Fold expression change in pupae vs adults (positive means higher in adults)", 
       y = "Negative log10-adjusted p-value") 


####### Violin plots for individual genes ########
plotCounts(pairwiseAnalysis,
           gene = "CVAR_01260",
           intgroup = c("group"))




###### stage specific #########
resultsCombo <- results(differentialExpression, 
                        name = "stagepupa.casteworker")
resultsComboDataframe <- data.frame(resultsCombo)
signficantPupalCastes <- filter(resultsComboDataframe,
                                padj <= 0.05)

plotCounts(differentialExpression,
           gene = "CVAR_10388",
           intgroup = c("caste",
                        "stage"))

###### PCA of samples ##########
object <- rlog(dataForDESeq, blind = TRUE)
plotPCA(object,
        intgroup = c("caste",
                     "stage"))

###### Add functional info to genes ##########
#### Create a dictionary linking cvar genes to functions from other species from OrthoFinder ####
orthogroupMembers <- read_delim(file = "03_Orthofinder/fasta/OrthoFinder/Results_Apr30/Orthogroups/Orthogroups.tsv",
                                delim = "\t")
cvarToFunctions <- concat.split.multiple(orthogroupMembers, 
                                         split.cols = "cvar_proteins", 
                                         seps = ",", 
                                         direction = "long") 
cvarToFunctions <- filter(cvarToFunctions, 
                          !is.na(cvar_proteins))
cvarToFunctionsLong <- cvarToFunctions %>% 
  pivot_longer(cols = c("acep_proteins",
                        "acol_proteins",
                        "amel_proteins",
                        "dmel_proteins",
                        "tcas_proteins"),
               names_to='species',
               values_to='gene') %>% 
  filter(!is.na(gene))
cvarToFunctionsLong <- concat.split.multiple(cvarToFunctionsLong, 
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
    select(seq.name)
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

# Combine the functions with the cvar-to-other-species table:
allFunctionInformation <- full_join(cvarToFunctionsLong,
                                    allProteomes, 
                                    by = c("species" = "species",
                                           "gene" = "geneName")) %>%
  distinct()

# Combine this table, which is transcript-based, with the gene-to-transcript key:
functionsDictionary <- full_join(transcriptsToGenes,
                                 allFunctionInformation,
                                 by = c("ID" = "cvar_proteins"))
test <- select(functionsDictionary, -c("Orthogroup"))
test <- pivot_wider(test,
                    names_from = species, 
                    values_from = c(gene, seq.name)) %>%
  dplyr::select(-c("gene_NA",
                   "gene_cvar_proteins",
                   "seq.name_NA",
                   "seq.name_cvar_proteins")) 

##### Combine the functions with the differential expression results: ####
resultsCasteWithFunctions <- resultsCasteDataframe
resultsCasteWithFunctions$gene_name <- row.names(resultsCasteWithFunctions)
resultsCasteWithFunctions <- full_join(resultsCasteWithFunctions, 
                                       functionsDictionary,
                                       by = c("gene_name" = "gene_name")) %>%
  dplyr::select(c(gene_name,
                  seq.name,
                  baseMean,
                  log2FoldChange,
                  lfcSE,
                  stat,
                  pvalue,
                  padj)) %>%
  distinct()




