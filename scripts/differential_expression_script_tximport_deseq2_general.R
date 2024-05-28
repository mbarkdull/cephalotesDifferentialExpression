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
transcriptsToGenes <- select(annotation, ID)
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
               size = 0.5, 
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
         x = paste("Fold expression change between",
                   tolower(unique(dataframe$contrast))), 
         y = "Negative log10-adjusted p-value") +
    theme_bw()
}

# List all of the pairwise comparisons:
resultsList <- list(pupaOnlyResultsDataframe,
                    adultOnlyResultsDataframe,
                    soldierOnlyResultsDataframe,
                    workerOnlyResultsDataframe)

# Use purrr to generate a volcano plot for each comparison:
volcanoPlotList <- purrr::map(resultsList, 
                              plotPairwiseContrasts)

# Plot them all together using patchwork:
patchwork::wrap_plots(volcanoPlotList, 
                      ncol = 2) + 
  patchwork::plot_annotation(title = 'Differential gene expression across caste and life stage contrasts',
                             theme = theme(plot.title = element_text(size = 19)))

#### Link the cvar genes to functions from other species based on OrthoFinder clustering ####
# Read in the orthogroup membership file:
orthogroupMembers <- read_delim(file = "03_Orthofinder/fasta/OrthoFinder/Results_Apr30/Orthogroups/Orthogroups.tsv",
                                delim = "\t")
# Split the column that lists genes from Cephalotes varians in a long direction (so each Cephalotes varians gene has one row linking it to other genes):
cvarToFunctions <- concat.split.multiple(orthogroupMembers, 
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
adultResultsFunctions <- combineFunctionsWithResults(adultOnlyResultsDataframe)
soldierResultsFunctions <- combineFunctionsWithResults(soldierOnlyResultsDataframe)
workerResultsFunctions <- combineFunctionsWithResults(workerOnlyResultsDataframe)

# Combine all of these dataframes and save them as a csv file:
allResultsForExport <- rbind(pupaResultsFunctions,
                             adultResultsFunctions,
                             soldierResultsFunctions,
                             workerResultsFunctions)

# Export the results as a .csv:
write_csv(allResultsForExport,
          file = "./finalResults/allDifferentialExpressionResultsAndFunctions.csv")




















# Run the differential expression analysis:
differentialExpression <- DESeq(dataForDESeq)

#### Generate dot plots for individual genes as a sanity check ####
plotCounts(pairwiseAnalysis,
           gene = "CVAR_01260",
           intgroup = c("group"))

#### Run a PCA of the samples ####
transformedData <- rlog(dataForPairwise, 
                        blind = TRUE)
DESeq2::plotPCA(transformedData,
                intgroup = c("group"),
                ntop = 10000)
