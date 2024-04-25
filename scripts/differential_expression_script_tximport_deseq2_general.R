# User supplies:
  # path to Salmon outputs
  # link to google sheet with sample metadata
  # genome annotation 
  # study design for DESeq to set up contrasts. 

# List all of the Salmon quantification output files:
salmon_outputs <- list.files(path = "/../../local/workdir/noah/JC/tissue_culture/results/salmon_4_27_23", 
                             full.names = TRUE,
                             recursive = TRUE,
                             pattern = "*sf")

# Make an object that connects transcript names of the salmon files to gene names from the genome annotation:
# Will want to change this:
tx2gene <- read.csv(salmon_outputs[1], 
                    sep = "\t" )
tx2gene <- subset(tx2gene, 
                  select = c(Name) )
tx2gene$gene_name <- tx2gene$Name

# Use tximport to import transcript-level abundance estimates across all samples:
txi <- tximport(files = salmon_outputs, 
                type = "salmon", 
                tx2gene = tx2gene, 
                countsFromAbundance = "lengthScaledTPM")

# Make a table of sample metadata, used to set up the contrasts in the differential expression analysis:
sample_type <- read_sheet("https://docs.google.com/spreadsheets/d/1ojUn849kIFZLH2bKsHsW34IZWv-bskkCUG6tlvkgEGQ/edit?usp=sharing")
meta <- data.frame(sample_type, 
                   row.names = colnames(txi$counts))

# Read in the genome annotation:
# Set up a named list used to filter the annotation file when read in:
my_filter <- list(type=c("mRNA"))
# Read in the annotation:
gff <- readGFF("/../../local/workdir/share/JCv2.OGS3.beta/JC_OGS3.beta.gff", 
               filter = my_filter, 
               columns = c("seqid"), 
               tags = c("ID", "Note"))
# Convert to a dataframe:
gene_names <- as.data.frame(gff)

# Use DESeq to run a Wald test for differential expression:
# Combine all relevant information into a DESeq dataset object:
dataForDESeq <- DESeqDataSetFromTximport(txi, 
                                         colData = meta, 
                                         design = ~ stage + treatment + stage:treatment)
# Run the differential expression analysis:
differentialExpression <- DESeq(dataForDESeq)
# List all contrasts:
resultsNames(differentialExpression)
# Extract the results for particular contrast of interest:
testresults_all <- results(differentialExpression, 
                           name = "treatment_t_vs_c" )
testresults_all_df <- data.frame(testresults_all)
