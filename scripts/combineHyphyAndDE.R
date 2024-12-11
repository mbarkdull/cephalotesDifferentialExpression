# Combine differential expression results with HYPHY results:
library(splitstackshape)
library(tidyverse)
library(googlesheets4)
library(eulerr)

#### Read in and combine the HYPHY results: ####
bustedPHResults <- read_csv(file = "allBustedPHResults.csv")
relaxPHResults <- read_csv(file = "allRelaxResults.csv")
hyphyResults <- full_join(bustedPHResults,
                          relaxPHResults,
                          by = c("orthogroup" = "orthogroup",
                                 "HOGmembers" = "HOGmembers"))

# Split the HOG members column and filter to get a column with CVAR transcript IDs for each result:
hyphyResults <- cSplit(hyphyResults, 
                       splitCols = "HOGmembers",
                       sep = "|",
                       direction = "long") %>% 
  filter(grepl("^CVAR", 
               HOGmembers)) %>%
  distinct()
# Turn all hyphens to underscores:
hyphyResults$HOGmembers <- gsub(pattern = "\\-",
                                replacement = "_",
                                hyphyResults$HOGmembers)

# Create a dictionary linking transcripts to genes:
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

transcriptsToGenes <- dplyr::select(annotation, ID)
transcriptsToGenes$gene_name <- annotation$Parent
# Turn all hyphens to underscores:
transcriptsToGenes$ID <- gsub(pattern = "\\-",
                              replacement = "_",
                              transcriptsToGenes$ID)

# Join the dictionary to the HYPHY results:
hyphyResults <- full_join(hyphyResults, 
                          transcriptsToGenes, 
                          by = c("HOGmembers" = "ID")) %>%
  filter(!is.na(inputFile.x) &
           !is.na(inputFile.y))

# Export:
write_tsv(hyphyResults,
          file = "./hyphyResultsForDE.tsv")

# Read in the differential expression results:
differentialExpression <- read_csv(file = "./finalResults/allDifferentialExpressionResultsAndFunctions.csv")

# Combine with the Hyphy results:
allResultsCombined <- full_join(hyphyResults,
                                differentialExpression,
                                by = c("gene_name" = "gene_name"))

#### Export top candidates to Google Drive ####
getTopTwenty <- function(specifiedContrast) {
  filterResults <- filter(allResultsCombined, 
                  contrast == specifiedContrast &
                    !is.na(padj)) 
  filterResultsArranged <- arrange(filterResults,
                                   log2FoldChange)
  topTen <- head(unique(filterResultsArranged$gene_name),
                 n = 10L)
  bottomTen <- tail(unique(filterResultsArranged$gene_name),
                    n = 10L)
  candidates <- c(topTen,
                  bottomTen)
  candidates <- filter(filterResultsArranged, 
                       gene_name %in% candidates)
  return(candidates)
}

lifeStageContrasts <- c("Worker pupae vs. soldier pupae",
                        "Adult workers vs. adult soldiers")
topTwenty <- purrr::map(lifeStageContrasts,
                        getTopTwenty)
topTwenty <- as.data.frame(do.call(rbind, topTwenty))

# For the pupal-to-adult contrasts, we only want the top genes that are NOT also differentially expressed in the other caste. 
soldierContrast <- filter(allResultsCombined, 
                          contrast == "Adult soldiers vs. soldier pupae" &
                            padj <= 0.05) 
soldierContrastArranged <- arrange(soldierContrast,
                                   log2FoldChange)

workerContrast <- filter(allResultsCombined, 
                         contrast == "Adult workers vs. worker pupae" &
                           padj <= 0.05) 
workerContrastArranged <- arrange(workerContrast,
                                  log2FoldChange)

soldierContrastExclusive <- subset(soldierContrastArranged, 
                                   !gene_name %in% workerContrastArranged$gene_name)

workerContrastExclusive <- subset(workerContrastArranged, 
                                  !gene_name %in% soldierContrastArranged$gene_name)

getTopTwentyLifestage <- function(filteredCandidates) {
  topTen <- head(unique(filteredCandidates$gene_name),
                 n = 10L)
  bottomTen <- tail(unique(filteredCandidates$gene_name),
                    n = 10L)
  candidates <- c(topTen,
                  bottomTen)
  candidates <- filter(filteredCandidates, 
                       gene_name %in% candidates)
  return(candidates)
}
soldierTopTwenty <- getTopTwentyLifestage(soldierContrastExclusive)
workerTopTwenty <- getTopTwentyLifestage(workerContrastExclusive)
allTopTwenty<- rbind(topTwenty,
                     soldierTopTwenty,
                     workerTopTwenty)


write_sheet(data = allTopTwenty,
            "https://docs.google.com/spreadsheets/d/1X6DWQjeL9CnfLEA38rErFgBRIdG9MnSTVS7swyaaIes/edit?usp=sharing")



#### Generate Euler diagrams: ####
library(eulerr)
library(gridExtra)

# Write a function to get the data for and generate an Euler diagram:
generateEulerDiagrams <- function(specificContrast) {
  # Get just genes that are differentially expressed in the chosen contrast:
  genesDifferentiallyExpressed <- filter(allResultsCombined, 
                                         contrast == specificContrast &
                                           padj < 0.05) %>%
    select(gene_name) %>%
    distinct() 
  genesDifferentiallyExpressed <- genesDifferentiallyExpressed$gene_name
  
  # Get positively selected genes:
  genesPositivelySelected <- filter(allResultsCombined, 
                                      testPvalueFDR < 0.05 &
                                      backgroundPvalueFDR >= 0.05 &
                                      differencePvalueFDR < 0.05) %>%
    select(gene_name) %>%
    distinct() 
  genesPositivelySelected <- genesPositivelySelected$gene_name
  
  # Get genes with a shift in selective intensity:
  genesSelectiveIntensityShift <- filter(allResultsCombined, 
                                           pValueFDR < 0.05) %>%
    select(gene_name) %>%
    distinct() 
  genesSelectiveIntensityShift <- genesSelectiveIntensityShift$gene_name
  
  # List those sets together:
  sets <- list(differentialExpression = genesDifferentiallyExpressed,
               positiveSelection = genesPositivelySelected,
               shiftedIntensity = genesSelectiveIntensityShift)
  
  # Euler plot of just DE and posSel:
  plot(euler(list(differentialExpression = genesDifferentiallyExpressed,
                  positiveSelection = genesPositivelySelected),
             shape = "ellipse"), 
       quantities = TRUE)
  
  # Euler plot of just DE and shifted:
  plot(euler(list(differentialExpression = genesDifferentiallyExpressed,
                  shiftedIntensity = genesSelectiveIntensityShift),
             shape = "ellipse"), 
       quantities = TRUE)
  
  # Euler plot of all three sets:
  diagram <- euler(sets, shape = "ellipse")
  plot(diagram, quantities = TRUE)
  return(diagram)
}

pupae <- generateEulerDiagrams("Worker pupae vs. soldier pupae")
pupae <- plot(pupae, 
              quantities = TRUE,
              main = "Pupal genes",
              fills = list(fill = c("#e9d677", 
                                    "#06a77d", 
                                    "#1378ed"), 
                           alpha = 0.5),
              labels = FALSE)
adults <- generateEulerDiagrams("Adult workers vs. adult soldiers")
adults <- plot(adults, 
               quantities = TRUE,
               main = "Adult genes",
               fills = list(fill = c("#e9d677", 
                                     "#06a77d", 
                                     "#1378ed"), 
                            alpha = 0.5),
               labels = FALSE)
soldiers <- generateEulerDiagrams("Adult soldiers vs. soldier pupae")
soldiers <- plot(soldiers, 
                 quantities = TRUE,
                 main = "Soldier genes",
                 fills = list(fill = c("#e9d677", 
                                       "#06a77d", 
                                       "#1378ed"), 
                              alpha = 0.5),
                 labels = FALSE)
workers <- generateEulerDiagrams("Adult workers vs. worker pupae")
workers <- plot(workers, 
                quantities = TRUE,
                main = "Worker genes",
                fills = list(fill = c("#e9d677", 
                                      "#06a77d", 
                                      "#1378ed"), 
                             alpha = 0.5),
                labels = FALSE) 
eulerPlotsList <- c(pupae,
                    adults,
                    soldiers,
                    workers)
allPlots <- grid.arrange(pupae,
                         adults,
                         soldiers,
                         workers, 
                         nrow = 2)
allPlots

ggsave("./images/eulerPlots.png",
       plot = allPlots,
       width = 12, 
       height = 8, 
       units = "in")








#### Break down pupal results in more detail: ####
minorUpregulated <- filter(allResultsCombined, 
                                       contrast == "Worker pupae vs. soldier pupae" &
                                         padj < 0.05 &
                             log2FoldChange < 1) %>%
  select(gene_name) %>%
  distinct() 
minorUpregulated <- minorUpregulated$gene_name

soldierUpregulated <- filter(allResultsCombined, 
                           contrast == "Worker pupae vs. soldier pupae" &
                             padj < 0.05 &
                             log2FoldChange > 1) %>%
  select(gene_name) %>%
  distinct() 
soldierUpregulated <- soldierUpregulated$gene_name


# Get positively selected genes:
genesPositivelySelected <- filter(allResultsCombined, 
                                  testPvalueFDR < 0.05 &
                                    backgroundPvalueFDR >= 0.05 &
                                    differencePvalueFDR < 0.05) %>%
  select(gene_name) %>%
  distinct() 
genesPositivelySelected <- genesPositivelySelected$gene_name

# Get genes with a shift in selective intensity:
genesIntensified <- filter(allResultsCombined, 
                           pValueFDR < 0.05 &
                             kValue > 1) %>%
  select(gene_name) %>%
  distinct() 
genesIntensified <- genesIntensified$gene_name

genesRelaxed <- filter(allResultsCombined, 
                       pValueFDR < 0.05 &
                         kValue < 1) %>%
  select(gene_name) %>%
  distinct() 
genesRelaxed <- genesRelaxed$gene_name

# List those sets together:
sets <- list(minorUpregulated = minorUpregulated,
             soldierUpregulated = soldierUpregulated,
             positiveSelection = genesPositivelySelected,
             intensified = genesIntensified,
             relaxed = genesRelaxed)

# Euler plot of all sets:
##### THIS IS INACCURATE
diagram <- euler(sets, 
                 shape = "ellipse")
plot(diagram, quantities = TRUE)



numberMinorUpregulated <- length(minorUpregulated)
numberSoldierUpregulated <- length(soldierUpregulated)

numberGenesPositivelySelected <- length(genesPositivelySelected)
numberMinorPositive <- length(intersect(minorUpregulated,
                                        genesPositivelySelected))
numberSoldierPositive <- length(intersect(soldierUpregulated,
                                          genesPositivelySelected))

numberGenesIntensified <- length(genesIntensified)
numberMinorIntensified <- length(intersect(minorUpregulated,
                                           genesIntensified))
numberSoldierIntensified <- length(intersect(soldierUpregulated,
                                             genesIntensified))

numberGenesRelaxed <- length(genesRelaxed)
numberMinorRelaxed <- length(intersect(minorUpregulated,
                                       genesRelaxed))
numberSoldierRelaxed <- length(intersect(soldierUpregulated,
                                         genesRelaxed))
numberGenes <- length(unique(allResultsCombined$gene_name))

# The following is from https://www.badgrammargoodsyntax.com/compbio/2017/12/16/compbio-017-is-your-overlap-significant
# Compare overlap with a hypergeometric test, where:
# q = (size of overlap-1)
# m = number of differential genes in experiment 1
# n = (total number of genes on platform - m)
# k = number of differential genes in experiment 2. 

# Minor and positive
pValueMinorPositive <- phyper(numberMinorPositive-1,
                              numberGenesPositivelySelected, 
                              numberGenes - numberGenesPositivelySelected, 
                              numberMinorUpregulated, 
                      lower.tail = FALSE, 
                      log.p = FALSE) 
# Minor and relaxed
pValueMinorRelaxed <- phyper(numberMinorRelaxed-1,
                              numberGenesRelaxed, 
                              numberGenes - numberGenesRelaxed, 
                              numberMinorUpregulated, 
                              lower.tail = FALSE, 
                              log.p = FALSE) 
# Minor and intensified
pValueMinorIntensified <- phyper(numberMinorIntensified-1,
                             numberGenesIntensified, 
                             numberGenes - numberGenesIntensified, 
                             numberMinorUpregulated, 
                             lower.tail = FALSE, 
                             log.p = FALSE) 

# Soldier and positive
pValueSoldierPositive <- phyper(numberSoldierPositive-1,
                              numberGenesPositivelySelected, 
                              numberGenes - numberGenesPositivelySelected, 
                              numberSoldierUpregulated, 
                              lower.tail = FALSE, 
                              log.p = FALSE) 
# Soldier and relaxed
pValueSoldierRelaxed <- phyper(numberSoldierRelaxed-1,
                             numberGenesRelaxed, 
                             numberGenes - numberGenesRelaxed, 
                             numberSoldierUpregulated, 
                             lower.tail = FALSE, 
                             log.p = FALSE) 
# Soldier and intensified
pValueSoldierIntensified <- phyper(numberSoldierIntensified-1,
                                 numberGenesIntensified, 
                                 numberGenes - numberGenesIntensified, 
                                 numberSoldierUpregulated, 
                                 lower.tail = FALSE, 
                                 log.p = FALSE) 



