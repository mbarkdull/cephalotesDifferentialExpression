# Combine differential expression results with HYPHY results:
library(splitstackshape)
library(tidyverse)
library(googlesheets4)

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

transcriptsToGenes <- select(annotation, ID)
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

