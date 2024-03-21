### Introduction to cBioPortal
### Demo on using the cBioPortal API to access data, then MAFtools to create plots
### First delivered 21/03/2024 via Zoom by Dr Andrew Mason (University of York, UK)
### Session delivered as an Elixir-UK training event

# check if installs are needed
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("cbioportalR", quietly = TRUE))
  BiocManager::install("cbioportalR")

if (!requireNamespace("maftools", quietly = TRUE))
  BiocManager::install("maftools")

# load libraries
library(dplyr)
library(cbioportalR)
library(maftools)

### 1) Accessing the cBioPortal API
### more information on cbioportalR - https://www.karissawhiting.com/cbioportalR/index.html

# connect to cBioPortal (requires internet)
# we will use public data only
# if you eventually use restricted access data, you will need a token
set_cbioportal_db("public")

# working with the API requires an internet connection
# poor connection can give intermittent/non-reproducible errors
# if you're having problems, check the connection
test_cbioportal_db()

# let's look at the available studies
studies <- available_studies()
head(studies)

# check for available bladder cancer studies
studies[grep("Bladder", studies$name),c("studyId","description")]

# we will use TCGA's muscle invasive bladder cancer study
# use the study Id and transpose for better viewing
get_study_info("blca_tcga_pub_2017") |> t()

# we can check the available molecular data
available_profiles(study_id = "blca_tcga_pub_2017")

# you can see that some of these appear as default (including online) and some are hidden
# let's check the IDs for the individual data types so we can work with them
available_profiles(study_id = "blca_tcga_pub_2017") |> pull(molecularProfileId)

# now we have used the API to find what data are available, now we can pull the data we want to our local machines
# this can be done by sample ID or by study ID, the latter can only do one study at a time
# you can pull an entire study with all its data with get_genetics_by_study()
# but we only want to access the mutation data for this demo
# this is the only command in this demo that takes time to run
tcga_blca_muts <- get_mutations_by_study(study_id = "blca_tcga_pub_2017")

# check the available data
colnames(tcga_blca_muts)

# visualise the first row
# format is just like what you can download in a text-based file - almost like a MAF file
head(tcga_blca_muts, 1) |> t()

### 2) Using MAFtools to visualise the data
### more information on MAFtools here - https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html
### and on the MAF file format - https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/

# cbioportalR gets lots of information, including stuff we don't need 
# for some reason it also uses camelCase, which MAFtools does not like
# convert column names
maf_conv_cols <- c("Hugo_Symbol", "Entrez_Gene_Id", "Sample_Key", "Patient_Key", 
                   "Molecular_Profile_ID", "Tumor_Sample_Barcode", "Patient_Id", 
                   "Study_Id", "Centre", "Mutation_Status", "Validation_Status", 
                   "t_alt_count", "t_ref_count", "Start_Position", "End_Position", 
                   "Reference_Allele", "Protein_Change", "Variant_Classification", 
                   "NCBI_Build", "Variant_Type", "Keyword", "Chromosome", 
                   "Tumor_Seq_Allele2", "RefSeq", "Protein_position", "Protein_end_position")

# set new column names and create a new VAF (variant allele frequence) column
colnames(tcga_blca_muts) <- maf_conv_cols
tcga_blca_muts <- mutate(tcga_blca_muts, VAF = t_alt_count/(t_alt_count + t_ref_count))

# now subset the cbioportalR-downloaded date for the columns you actually want in the right order
maf_cols <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", 
              "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification", 
              "Variant_Type", "Tumor_Sample_Barcode", "Protein_Change", "VAF", 
              "NCBI_Build", "Patient_Id")
tcga_blca_muts_req_maf_cols <- tcga_blca_muts[, maf_cols]

# convert to MAF
tcga_blca_maf <- read.maf(tcga_blca_muts_req_maf_cols)

# we get some warnings here about genes which are commonly mutated because they are massive
# now in MAF, we can extract summary information on the whole cohort, by patient and by gene
tcga_blca_maf
getSampleSummary(tcga_blca_maf)

gene_tots <- getGeneSummary(tcga_blca_maf)
gene_tots[grep("FGFR3", gene_tots$Hugo_Symbol),]

# we can also make a nice summary plot
plotmafSummary(maf = tcga_blca_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

# the problem is all of those genes which are not likely to be important in the cancer, but they are mutated lots

### 3) Filter the MAF file using the oncoKB gene list
### gene sets are updated regularly and you can download from here - https://www.oncokb.org/cancer-genes
### for this workshop I have extracted some pertinent columns before loading into R
### you can access my file here - https://asmasonomics.github.io/assets/coursefiles/2024-03-21_cBioPortal/2024-02-08_oncoKB_list.tsv

# load in the oncoKB file
oncoKB <- read.table("2024-02-08_oncoKB_list.tsv", header = TRUE)
head(oncoKB)

# extract just the genes
oncoKB_genes <- oncoKB$Hugo_Symbol

# filter cbioportalR download using oncoKB gene list
tcga_blca_muts_req_maf_cols_oncoKB <- tcga_blca_muts_req_maf_cols[tcga_blca_muts_req_maf_cols$Hugo_Symbol %in% oncoKB_genes,]

# build new maf object
tcga_blca_maf_oncoKB <- read.maf(tcga_blca_muts_req_maf_cols_oncoKB)

# summary and plot
tcga_blca_maf_oncoKB
plotmafSummary(maf = tcga_blca_maf_oncoKB, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

## time to make some nice plots using the MAFtools functionalities
# Draw an oncoplot/oncoprint
# by top number of genes
oncoplot(maf = tcga_blca_maf_oncoKB, top = 10)

# by genes mutated above a given cohort frequency
oncoplot(maf = tcga_blca_maf_oncoKB, minMut = 0.1)

# by a specific gene list
oncoplot(maf = tcga_blca_maf_oncoKB, genes = c("KDM6A", "ARID1A", "KMT2D"))

# be careful here, as oncoplot() only shows altered samples by default
# so it looks like the full cohort is mutated at first glance
oncoplot(maf = tcga_blca_maf_oncoKB, genes = c("KDM6A", "ARID1A", "KMT2D"), removeNonMutated = FALSE)

# lollipop plots
lollipopPlot(maf = tcga_blca_maf_oncoKB, gene = "FGFR3", AACol = "Protein_Change", labelPos = 249)

# assess mutual exclusivity
mut_excl <- somaticInteractions(maf = tcga_blca_maf_oncoKB, top = 20, pvalue = c(0.05, 0.1))


### DEMO COMPLETE
# Thanks for getting this far - I hope it was informative!
# For the purposes of the demo we only used a single dataset to match the website demo earlier in the course.
# The real power of the API is that you can combine multiple sequencing modalities from the same project
# and combine multiple projects. Check out the cbioportalR documentation for tips on this (link at the top).

### Andrew Mason
### The University of York
### https://asmasonomics.github.io
### 21/03/2024
### Elixir-UK