---
layout: page
permalink: /courses/BABS4_Biochem_DataWorkshop3_March2024
---

![BABS4 banner](/assets/images/babs4_banner.jpg){:class="img-responsive"}

<span style="font-size:1.6em;">**BABS4 - Data Workshop 3**</span><br/>

<p align="justify">
Welcome to Data Workshop 3! This is the first of two RNAseq workshops as part of the BABS4 (66I) "Gene expression and biochemical interactions strand". If you're from the future and are completing Workshop 4,  <a href="https://asmasonomics.github.io/courses/BABS4_Biochem_DataWorkshop4_March2024">please follow this link to the correct material</a>.<br/>
The material below will cover many of the R commands needed to fully analyse these data. If you're feeling a bit rusty, <a href="https://3mmarand.github.io/R4BABS/r4babs4/week-1/workshop.html">please consult Emma's material from the BABS4 core data workshop in week 1</a>.<br/>
</p>

### Introduction
<p align="justify">
In this workshop you will work on publicly available RNAseq data from wildtype <i>Haemophilus influenzae</i>. Remember, this bacterium is naturally competent. During stress, gene networks which regulate competence should be activated to help with survival in difficult conditions.<br/>
In their <a href="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0217255">2020 paper, Black <i>et al.,</i></a> aimed to measure this competence response by switching happily proliferating <i>Haemophilus influenzae</i> in BHI medium, to a "starvation" medium called MIV. They took RNA at multiple time points and had cultures in triplicate to assess the deviation in response and enable statistical testing. We will focus initially on the control (<i>t</i>=0 minutes in MIV; kw20-MIV0) compared with <i>t</i>=30 minutes (kw20-MIV2).<br/>
</p>
![Hi experimental setup](/assets/images/2024-03-11_66I-DW3_Hi_exp_setup.png){:class="img-responsive"} <br/>
<p align="justify">
Hopefully this all sounds quite familiar! If not, after this morning, make time to watch the videos we've generated to support these two workshops. These are embedded below and on the VLE.<br/><br/>
You can also download the introductory slides <a href="/assets/coursefiles/2024-03_66I/66I-DW3_introductory_slides.pdf" download>here as a PDF</a>.<br/>
</p>
[Introduction to transcriptomics](#introduction-to-transcriptomics)<br/>
[How has the data in this workshop been processed so far?](#data-processing-before-this-workshop)
<br/>

### The workshop
#### Set up your RStudio project
<p align="justify">
Start RStudio from the Start menu.<br/><br/>
Make a new RStudio project. Put it in a sensible place with a sensible name. Remember that if you are using a university machine, select a location via the M/H drives, not just using the "Documents" shortcut. This <b>will</b> create problems for you!<br/><br/>
Use the Files pane to make subdirectories for your <code>raw_data</code>, <code>proc_data</code> and <code>plots</code>. These names are suggestions only. Make a new script file, perhaps <code>babs4_rnaseq_workshop.R</code>, to complete your work.<br/><br/>
Remember, you will eventually submit a <b>single</b> RStudio project which should represent the entire project over both workshops (3&4).<br/><br/>
</p>

#### Access the data
<p align="justify">
During workshop 3 you will need 3 <i>Haemophilus influenzae</i> datafiles. Download and add to your <code>raw_data</code> subdirectory. From here on, <i>Haemophilus influenzae</i> will be <b>abbreviated to Hi</b>.<br/><br/>

<a href="/assets/coursefiles/2024-03_66I/Hi_PRJNA293882_counts.tsv" download>Hi_PRJNA293882_counts.tsv</a>. These are the RNAseq data from Hi. The data includes 11 different conditions, each with three replicates (explained below). This is a tab separated file with sample names in the first row and feature IDs in the first column. The data are counts of the number of reads.<br/>
kw20 is the wildtype strain. OD<sub>600</sub> is a measure of cell suspension density. Sxy is a transcription factor (encoded by <i>tfoX</i>) which regulates competence. The Sxy- strain is a null mutant strain where Sxy is non-functional. <br/><br/>
</p>

| Condition | Strain | Media | OD<sub>600</sub> | time (mins) | replicates |
| --- | --- | --- | --- | --- | --- |
| kw20-BHI1 | kw20 | BHI | 0.02 | NA | -F,-G,-H |
| kw20-BHI2 | kw20 | BHI | 0.60 | NA | -F,-G,-H |
| kw20-BHI3 | kw20 | BHI | 1.00 | NA | -F,-G,-H |
| kw20-MIV0 | kw20 | MIV | 0.25 | 0 | -A,-B,-C |
| kw20-MIV1 | kw20 | MIV | 0.25 | 10 | -A,-B,-C |
| kw20-MIV2 | kw20 | MIV | 0.25 | 30 | -A,-B,-C |
| kw20-MIV3 | kw20 | MIV | 0.25 | 100 | -A,-B,-C |
| sxyx-MIV0 | Sxy- | MIV | 0.25 | 0 | -A,-B,-D |
| sxyx-MIV1 | Sxy- | MIV | 0.25 | 10 | -A,-B,-D |
| sxyx-MIV2 | Sxy- | MIV | 0.25 | 30 | -A,-B,-D |
| sxyx-MIV3 | Sxy- | MIV | 0.25 | 100 | -A,-B,-D |

<p align="justify">
<a href="/assets/coursefiles/2024-03_66I/Hi_feature_names.tsv" download>Hi_feature_names.tsv</a>. This is a tab separated file with the Hi feature ID in column 1, the symbol for that feature in column 2 (if there is one, else "." is used as a placeholder), and a description in column 3 (vague at times!). There is a header. <br/><br/>
<a href="/assets/coursefiles/2024-03_66I/Hi_feature_locations.bed" download>Hi_feature_locations.bed</a>. This is a tab separated file with a particular format called BED format. This means the columns have standard data types, so there is no header. Column 1 is the sequence a feature is located. Column 2 is the start position and column 3 the end position. Column 4 has the Hi feature ID. Column 5 has the biotype of the feature. Column 6 has the strand.<br/><br/> 
</p>

#### Load your libraries
<p align="justify">
We need the following packages for this workshop:<br/>
<ul>
	<li><code>tidyverse</code> <a href="https://doi.org/10.21105/joss.01686">Wickham <i>et al.,</i> 2019</a></li>
	<li><code>DESeq2</code> <a href="https://doi.org/doi:10.18129/B9.bioc.DESeq2">Love <i>et al.,</i> 2014</a></li>
	<li><code>EnhancedVolcano</code> <a href="https://github.com/kevinblighe/EnhancedVolcano">Blighe <i>et al.,</i> 2018</a></li>
	<li><code>ggplot2</code> <a href="https://ggplot2-book.org/">Wickham <i>et al.,</i> 2016</a></li>
	<li><code>dplyr</code> <a href="https://doi.org/10.21105/joss.01686">Wickham <i>et al.,</i> 2019</a> (this is part of tidyverse, but I always find loading it explicably helps with redundancy issues later on)</li>
</ul>
<br/>
Load these libraries.
</p>

```R

library(tidyverse)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(EnhancedVolcano)

# if EnhancedVolcano is not found as a package, do the following:
install.packages("remotes")
library("remotes")
# press enter to skip updates when prompted
remotes::install_github("kevinblighe/EnhancedVolcano")
library(EnhancedVolcano)


```
<br/>

#### Load datasets
```R

# counts data
counts <- read.table("raw_data/Hi_PRJNA293882_counts.tsv", row.names = 1, header = TRUE)

# feature IDs and symbols where possible
featname <- read.table("raw_data/Hi_feature_names.tsv", row.names = 1, header = TRUE)

# feature locations, setting the column names for BED format
featlocs <- read.table("raw_data/Hi_feature_locations.bed", header = FALSE, col.names = c("chr","start","end","feat_ID","biotype","strand"))

```
<br/>

#### Explore the datasets
<p align="justify">
Take some time to explore the data using <code>view()</code>, <code>summary()</code> and <code>head()</code>. Does the data look right? Do the column names make sense? Do you understand what each file is showing you?<br/><br/>
Let's start to see what the data looks like, and how we can bring the datasets together.<br/>
</p>

```R

# explore the counts data by making a histogram of the read counts in the first column
counts |> ggplot(aes(x = kw20.BHI1.F)) + geom_histogram()

```
![kw20.BHI1.F histogram](/assets/coursefiles/2024-03_66I/plots/03_explore_001.png){:class="img-responsive"}
<p align="justify">
<br/>
<details>
   <summary>What is this histogram telling you about the distribution of read counts in this sample?</summary>
   Most features have a "low" read count. But a few features have very high read counts.<br/>
   So how can we make the plot easier to interpret?
</details><br/>
</p>

```R

# try a log10 transform of the same column
counts |> ggplot(aes(x = log10(kw20.BHI1.F + 1))) + geom_histogram()

# we do +1 inside the transform because log10(1) is 0 
# this means the graph is not impacted by 0s and small decimals (missing and negatives respectively)
# less of a problem here because read counts are integers, but it may be relevant later on.

```
![kw20.BHI1.F log10 histogram](/assets/coursefiles/2024-03_66I/plots/03_explore_002.png){:class="img-responsive"}

<p align="justify">
<br/>
Much clearer! The count data has a fairly normal distribution with the means just above a thousand counts (10<sup>3</sup> on the x axis).<br/><br/>
Now let's explore the other two datasets we loaded in.<br/>
</p>

```R

# take a look at the top of the feature ID file with symbols and gene descriptions
head(featname)

# OK, so we get some symbols - maybe we can look for the muA, muB and gam genes we have been working on?
# use grep (pattern matcher) to check the symbol column for the three genes
featname[grep("muA|muB|gam", featname$symbol),]

# Success!
# those IDs suggest the features are quite close together. Let's use the location data to check.
# use the same pattern match as before, but extract the IDs and use them as a search term in the location data
mu_feats <- rownames(featname[grep("muA|muB|gam", featname$symbol),])
mu_feats_grep <- paste(mu_feats, collapse = "|")
featlocs[grep(mu_feats_grep, featlocs$feat_ID),]

# Yes - very close together
# in the BLAST/PHASTER workshop we were trying to annotate the prophage region annotated by PHASTER - we can use the location data now to speed this up
# use the coordinates from PHASTER (data workshop 2) to extract all IDs in the region
hi_prophage_region <- featlocs |> filter(between(start, 1558774, 1597183))

# think how you can use the code above to take these IDs from the prophage region to extract the gene symbols
# this will help your annotation of the prophage region from last session

# so now we have counts and feature IDs for genes we care about - what does gam look like?
# filter the counts data for the gam feature ID, transpose to turn the row into a column, and then convert back to a dataframe
gam_counts <- counts |> filter(row.names(counts) %in% c("gene-HI_1483")) |> 
  t() |> 
  as.data.frame()

# set the column name
colnames(gam_counts) <- "gamcounts"

# plot
ggplot(gam_counts, aes(x = rownames(gam_counts), y = gamcounts)) + 
  geom_bar(stat = "identity") + 
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = "samples", y = "Hi gam total reads")

```
![gam bar chart](/assets/coursefiles/2024-03_66I/plots/03_explore_003.png){:class="img-responsive"}

<p align="justify">
<br/>
Great - we've already started to explore one of our genes of interest. But, we can't just compare raw counts on their own.<br/>
<details>
   <summary>Why not? (Hint: are there any patterns in the bar chart related to the samples?)</summary>
   There are patterns of high/low/medium <i>etc</i> across the replicates. We don't know whether this is biological or technical at this point.<br/>
   But, before we can compare between genes, we need to know if each sample was sequenced to the same total depth (<i>i.e.</i> total number of reads). Otherwise, samples with more reads may just have higher counts - there may be no biological difference. Let's check.<br/>
</details><br/>
</p>

```R

# create a dataframe with sums for each column
count_tots <- data.frame(colSums(counts))
colnames(count_tots) <- "sum"

# let's check the ratio between the highest and lowest total read counts
max_reads <- max(count_tots$sum)
min_reads <- min(count_tots$sum)
max_reads / min_reads

# the highest has almost 4x as many reads as the lowest - no wonder there were differences in our last bar chart

# plot all the read totals to look for variance
ggplot(count_tots, aes(x = rownames(count_tots), y = sum)) + 
  geom_bar(stat = "identity") + 
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = "samples", y = "total reads")

```
![read count totals bar chart](/assets/coursefiles/2024-03_66I/plots/03_explore_004.png){:class="img-responsive"}

<p align="justify">
<br/>
When we run differential expression analysis (DEA) later, DESeq2 accounts for the differences in read counts automatically. But, when we want to plot individual genes and make nice graphs we need to normalise the data. We are going to convert our read count data to <b>TPM</b> - transcripts per million.<br/><br/>
I briefly covered TPM in the Introduction to Transcriptomics video, but I also have a <a href="https://www.youtube.com/watch?v=3Pe9xcGF_Wo&t=5s">specific video on RNAseq units</a> if you want to know a bit more.<br/>
</p>
<br/>

#### Normalise the count data
<p align="justify">
Normalising count data gives us a fairer way to compare between samples and even between genes in the samples.<br/><br/>
<details>
   <summary>What three things influence the read count for any feature?</summary>
   <ol>
     <li><b>Expression level</b>. The biologically interesting one! Features which have higher expression should generate more sequencing reads. </li>
	 <li><b>Feature length</b>. Longer features expressed at similar levels to shorter ones would generate more reads. </li>
	 <li><b>Total sequencing depth</b>. As we saw above, physically generating more reads may give artificially higher counts. </li>
   </ol>
</details><br/>
The TPM of any feature can be calculated by:
<ol>
<li>dividing the feature counts by the feature length (<i>i.e.</i> counts per base)</li>
<li>dividing counts/base by the total counts for the sample</li>
<li>multiplying by 1 million to get numbers we can work with (otherwise everything is very small)</li>
</ol>
So let's get on with it. 
</p>

```R

# let's first subset our data to only our core comparison: the impact of starvation media on Hi
core_counts <- counts |> select(starts_with("kw20.MIV"))

# we have the position information on all features in featlocs, so we can create a new column for the feature lengths
featlens <- data.frame(mutate(featlocs, feat_len = (featlocs$end+1)-featlocs$start))
rownames(featlens) <- featlens$feat_ID

# we now need to combine our datasets so we can divide the feature counts by feature lengths for each sample
# to join two datasets you need a key to link associated data together
# that's why we've made the rownames the feature IDs (feat_ID) in both dataframes (hence column 0 in the merge)
withlens <- merge(core_counts, featlens, by=0)
rownames(withlens) <- withlens$Row.names

# create a new dataframe with the sample counts divided by the lengths
counts_per_base <- subset(withlens, select = colnames(core_counts)) / withlens$feat_len

# now we use apply() to divide each value by the sum of its column, then multiply by a million to make it TPM
# the 2 indicates the function is applied on columns, rather than rows (where it would be 1)
tpms <- data.frame(apply(counts_per_base, 2, function(x){(x/sum(x))*1000000}))

# check each column now totals 1 million 
colSums(tpms)

# round the data to 2 decimal places to make it more human readable (does each column still total 1 million?)
tpms <- round(tpms, 2)

# now, let's make a quick assessment as to how consistent TPMs and counts are with each other, using gam
# extract just the gam TPMs, as we did with counts above
gam_tpms <- tpms |> filter(row.names(tpms) %in% c("gene-HI_1483")) |> t() |> as.data.frame()
colnames(gam_tpms) <- "gamtpms"

# merge the TPM and count gam subset dataframes and plot as a scatter
gamcor <- merge(gam_tpms, gam_counts, by=0)
rownames(gamcor) <- gamcor$Row.names

ggplot(gamcor, aes(x=gamtpms, y=gamcounts)) + 
  geom_point() + 
  geom_text_repel(label=rownames(gamcor)) + 
  labs(x = "gam TPMs", y = "gam read counts") +
  ylim(0, max(gamcor$gamcounts))
  
# calculate the correlation
cor(gamcor$gamcounts, gamcor$gamtpms, method = c("pearson"))

```
![gam count vs tpm](/assets/coursefiles/2024-03_66I/plots/03_normalise_001.png){:class="img-responsive"}

<p align="justify">
<br/>
The correlation is pretty good for this gene (and it's usually quite good), though both MIV0.A and MIV3.C have quite different TPMs to what their read counts might indicate.<br/>
</p>
<br/>

#### Principal component analysis
<p align="justify">
Before we do our stat testing for differences between conditions, it is always a good idea to perform principal component analysis (PCA).<br/><br/>
Our data is "high dimensional data", because we have lots more rows than we could plot simultaneously to describe the samples. We can just about interpret plots on 3 axes, but no further. PCA is a method for "dimension reduction". It works by looking at correlated variance across the samples in each row of data. In biology this is quite easy to rationalise. Genes work in pathways, or may be regulated by the same transcription factor, so you can imagine that these genes would go up and down in expression together (<i>i.e.</i> they are not each truly independent measures of the samples). This means you can reduce your data to a smaller number of "principal components" where correlated measurements are collapsed together. You can then plot the most informative principal components (<i>i.e.</i> those accounting for the most correlated variance in the data) to get a feel for how your samples group together (or not).<br/><br/>
The important thing here is that differential expression (like with any stat test) is only able to detect significant changes if your data are not too noisy (<i>i.e.</i> high in variance). So if your samples do not group together nicely by PCA it may be informative for removing samples (outliers) or explain why some comparisons will not give significant differences.<br/><br/>
Hopefully this becomes clearer with a plot.<br/>
</p>

```R

# perform PCA on the tpm dataframe
res_pca <- princomp(tpms)

# extract the variance accounted for by each principal component (squared standard deviation as proportion of total)
pca_var <- round((data.frame(res_pca$sdev^2/sum(res_pca$sdev^2)))*100,2)
colnames(pca_var) <- "var_perc"
pca_var

# most of the variance is accounted for by PCs 1 and 2. You can now plot a scatter of these samples to see if they group based on these components
# extract PCs 1 and 2
pca_comps <- data.frame(res_pca$loadings[, c("Comp.1", "Comp.2")])

# create axis labels which include the % variance accounted for by each component
pca_x <- paste(c("PC1 ("), pca_var["Comp.1",], c("%)"))
pca_y <- paste(c("PC2 ("), pca_var["Comp.2",], c("%)"))

# scatter
ggplot(pca_comps, aes(x=Comp.1, y=Comp.2)) + 
  geom_point() + 
  geom_text_repel(label=rownames(data.frame(res_pca$loadings[, 1:2]))) +
  labs(x = pca_x, y = pca_y)

# you could plot PC1 vs PC3 or PC2 vs PC3 to further interogate the data and the relationships

```
![PCA](/assets/coursefiles/2024-03_66I/plots/03_normalise_002.png){:class="img-responsive"}
<p align="justify">
<br/>
<details>
   <summary>What does this plot tell you about the samples?</summary>
   <ol>
     <li>The control samples are tightly grouped and well away from the rest. This suggests consistency in this group (good for stats), and that exposure to the starvation media <b>does</b> appear to have a measurable impact on the Hi transcriptome - good news!</li>
	 <li>The t=10mins (MIV1) and t=30mins (MIV2) conditions are quite similar by this plot.</li>
	 <li>The t=100mins (MIV3) samples are very variable. This may reflect differential responses of the cultures (which <i>should</i> be genetically identical...) but also noise in the data.</li>
   </ol>
</details>
<br/>
Genes with higher expression also have higher variance, so a log transformation will reduce that weighting. Try rerunning the PCA code with: <code>res_pca_log10 <- princomp(log10(tpms+1))</code><br/><br/>
<details>
   <summary>Does the log transformation give a different result, or tell you anything new about the data?</summary>
   The log transformation strongly suggests that the "C" replicates across conditions were quite different to the rest. It's possible these were processed all at the same time, different to those of A and B, maybe even by a different person. We may be able to control for this technical artefact during differential expression.
</details>
<br/>
Your PCA plots would be good ones to save. Make sure you give them sensible names and put them in your <code>plots</code> directory to keep things organised.<br/>
You can use the Export function on the Plots pane of RStudio, or code it. See below the PCA image.<br/>
<br/>
</p>

![log10 PCA](/assets/coursefiles/2024-03_66I/plots/03_normalise_003.png){:class="img-responsive"}
<br/>

```R 

ggplot(pca_comps, aes(x=Comp.1, y=Comp.2)) + 
  geom_point() + 
  geom_text_repel(label=rownames(data.frame(res_pca$loadings[, 1:2]))) +
  labs(x = pca_x, y = pca_y)
ggsave("plots/tpm_pca.pdf")

```

<br/><br/>

#### Differential expression analysis
<p align="justify">
We're now going to focus in on the main comparison for this workshop: the impact of the changed media after 30 mins (<i>i.e.</i> comparing <i>t</i>=30 against <i>t</i>=0; MIV2 vs MIV0).
<br/>
</p>

```R

# reduce counts matrix to just the groups we want
comp_counts <- counts[, grep("kw20.MIV0|kw20.MIV2", colnames(counts))]

# DESeq2 needs a dataframe detailing the experimental setup
sample_info <- data.frame(colnames(comp_counts))
colnames(sample_info) <- c("sample")
rownames(sample_info) <- sample_info$sample
sample_info

# currently sample_info looks very sparse, but we actually have all the information we need because of our consistent sample naming
# we can split the sample name and this gives us the genotype (kw20), condition (MIV0, MIV2) and replicates (A,B,C) 
sample_info <- sample_info |> separate(sample, c("genotype", "condition", "replicate"))

# crucial check needed now - are the columns in the counts data all found in the rownames of the sample_info (and in the same order)
all(rownames(sample_info) %in% colnames(comp_counts))
all(rownames(sample_info) == colnames(comp_counts))

# these are both TRUE, which is good. If these are FALSE, you need to reorder your columns/rows to make them match

# now we can run the differential expression
# the design parameter tells DESeq2 what comparison you want to make
dds <- DESeqDataSetFromMatrix(countData = comp_counts, colData = sample_info, design = ~ condition)

# explicitly set the control condition so fold changes are in the direction you want/makes biological sense
dds$condition <- relevel(dds$condition, ref = "MIV0")

# run DESeq2
dds <- DESeq(dds)
dds_results <- results(dds)

# get some quick summaries of 'significant' changes and take a look at the results
summary(results(dds, alpha=0.05, lfcThreshold = 1))
dds_results

```
<p align="justify">
<br/>
Good news! There are lots of significant changes. We use as log<sub>2</sub> fold change as it makes positive and negative changes symmetrical - see the table below for an explainer. We also need to have an 'adjusted' p value, not just the normal p value. The adjustment accounts for multiple testing of the data, where chance differences could come through as false positives. Effectively the adjustment makes our threshold for significance more stringent.<br/><br/>
When we compare values in two groups, we divide our experimental group by our control group to generate fold change values. So if something goes up in the experimental group, fold change will be greater than 1. Doubling expression or halving expression should be equivalent in our heads (they are <i>reciprocal</i> values), but they differ in their arithmetic distance from 1 (see the table). This means increasing expression could have a fold change value from >1 all the way to infinity. But a decrease can only go from <1 to approaching 0 - see the left histogram below.<br/>
A log<sub>2</sub> transformation of the fold change brings the equivalence back, making it symmetrical around 0 (see right hand histogram). This is why volcano plots use a log<sub>2</sub> x axis.<br/>
</p>

| gene | test<sub>avg TPM</sub> | control<sub>avg TPM</sub> | test/control | log<sub>2</sub>(test/control) |
| --- | --- | --- | --- | --- | 
| geneA | 100 | 50 | 2 | 1 | 
| geneB | 50 | 100 | 0.5 | -1 | 
| geneC | 20 | 20 | 1 | 0 |

![FC log2FC comparison](/assets/coursefiles/2024-03_66I/plots/03_dea_001.png){:class="img-responsive"}

<p align="justify">
That's the maths explainer over. Let's get back to interrogating the results of our DEA.<br/>
</p>

```R 

# currently our results all have the feature ID, so let's add in the gene symbol and TPMs
# first merge with the featnames to get symbols
dds_results <- merge(as.data.frame(dds_results), featname, by=0)
rownames(dds_results) <- dds_results$Row.names
dds_results <- dds_results[,-1]

# then subset the TPMs, calculate log2FC and merge together
comp_tpms <- tpms[, grep("kw20.MIV0|kw20.MIV2", colnames(tpms))]
comp_tpms <- mutate(comp_tpms, MIV0_avg = rowMeans(select(comp_tpms, contains("MIV0"))), 
                    MIV2_avg = rowMeans(select(comp_tpms, contains("MIV2"))))
comp_tpms <- mutate(comp_tpms, log2FC = log2((MIV2_avg+1)/(MIV0_avg+1)))
comp_red_tpms <- comp_tpms |> select(-starts_with("kw20.MIV"))

dds_tpm <- merge(dds_results, comp_red_tpms, by=0)
rownames(dds_tpm) <- dds_tpm$Row.names
dds_tpm <- dds_tpm[,-1]

# we now have log2FC values calculated from the DESeq2-normalised counts and our TPMs
# are they well correlated?
cor(dds_tpm$log2FoldChange, dds_tpm$log2FC, method = c("pearson"))

# now you have some context (particularly gene names), order by most significant and take a look at the top 20
head(dds_tpm[order(dds_tpm$padj),], 20)

```

<p align="justify">
<br/>
<details>
   <summary>Are there genes here which make sense?</summary>
   The <i>com</i> genes are <b>competence</b> genes, so that makes a lot of sense. We also have DNA (<i>dprA</i>) or recombination (<i>rec2</i>) regulators.<br>
   There are also lots of genes without symbols. This is fine! It may represent new biology to explore. Some of the descriptions for these unnamed ones make sense though - protein transport <i>etc</i>.
   <br/>
</details>
<br/>
Looking at a table like this is tricky, so let's make a volcano plot!<br/>
</p>

```R 

# first extract a list of genes (with symbols) to annotate the biggest changes on our volcano
# subset the dataframe where the symbols column does not equal "."
# we need the square brackets here otherwise "." takes on its special meaning of "match any character"
withsymbols <- dds_tpm[- grep("[.]", dds_tpm$symbol),]

# take the top 60 based on absolute log2FC (so this includes big negatives if there are any)
siggenes <- head(withsymbols |> arrange(desc(abs(log2FC))), 60)$symbol

# make the volcano plot - you can personalise these forever, check out - https://github.com/kevinblighe/EnhancedVolcano
# the cutoffs define the dotted line thresholds and default colours
# the second line is all about the labels used 
# the third line are some personal preferences on the plot layout - play with these!
EnhancedVolcano(dds_tpm, x = "log2FC", y = "padj", FCcutoff = 1, pCutoff = 0.05,
                lab = dds_tpm$symbol, selectLab = siggenes, labSize = 4.5, max.overlaps = 1000, drawConnectors = TRUE,
                legendPosition = 0, gridlines.major = FALSE, gridlines.minor = FALSE)

# there are a lot of big, significant changes in this dataset - I would perhaps suggest making the FCcutoff harsher. The plot below is for 2.5

```
![MIV0 MIV2 volcano](/assets/coursefiles/2024-03_66I/plots/03_dea_002.png){:class="img-responsive"}

<p align="justify">
<br/>
The volcano plot highlights the large number of very significant changes in the dataset after culture in the starvation media, MIV. We've already mentioned some of the genes such as <i>dprA</i> and <i>comA</i>. Next session we will work more on exactly what the biological impact is, but for now investigate a few more of the groups. This can be from the description column in the results, or by some tactical research on google.<br/>
<details>
   <summary>What do the <i>pur</i> and <i>trp</i> genes do?</summary>
   Genes such as <i>purC</i> are involved in making nucleotides. Similarly, genes like <i>trpC</i> make the amino acid tryptophan. MIV lacks both of these metabolites.<br/><br/>
</details>
<details>
   <summary>Don't forget about the downregulated genes - any common themes?</summary>
   The <i>rpL</i> and <i>rpS</i> genes encode the large and small ribosomal subunits. <i>fis</i> activates rRNA transcription. It looks like there is a general downregulation of translation machinery.<br/><br/>
</details>
<details>  
   <summary>Do you recognise <i>tfoX</i>?</summary>
   <i>tfoX</i> encodes the Sxy protein which is involved in competence response (directly via CRP). It is upregulated here and we have another Hi dataset which has a Sxy null mutant. Can you make any hypotheses about what may happen in Sxy mutants when they are grown in MIV?<br/><br/>
</details>
<details>   
   <summary>What about <i>muA</i>, <i>muB</i> and <i>gam</i>?</summary>
   They don't come up in this most significant list. You could change the <code>siggenes</code> variable to include only these genes and see where they are on the volcano. Or use your <code>grep</code> commands from earlier in the workshop to extract these genes from the results file.<br/>
   Modest but statistically significant upregulation is seen across all genes. So there is constituitive expression, enhanced in MIV.<br/><br/>
</details>
<br/>
</p>

#### Finishing up for today
<p align="justify">
Fantastic work! Lots of R coding and you've taken a big dataset, normalised it, extracted the columns you're biologically interested in, and performed differential expression analysis. Next session you'll be focused on what this actually means, relating this to genome location and how these genes may be regulated through regulons. Before that you need to make sure you've got everything saved so you don't need to do it again.<br/><br/>
Make sure you have saved the plots you want to save - particularly the PCA and volcano plots - and also the relevant datasets which you've modified.
<br/>
</p>

```R 

write.table(dds_tpm, file="proc_data/MIV0vsMIV2_DEA_full_results.tsv", sep="\t", row.names=TRUE)
write.table(tpms, file="proc_data/full_dataset_TPMs.tsv", sep="\t", row.names=TRUE)
write.table(comp_tpms, file="proc_data/MIV0vsMIV2_TPMs.tsv", sep="\t", row.names=TRUE)

# and anything else you've done which is relevant to you

```
<br/>
<p align="justify">
<b>GREAT WORK!</b> See you next week!
<br/>
</p>
<br/><br/><br/>

#### Consolidation and extending the workshop analysis (if you're interested)
<p align="justify">
If you want to check your understanding, the best bet is to look back through the workshop section headers and make sure you understand what each step was for. It's easy to put your head down and ignore the biology - so always bring it back to that. Remember, the question was how <i>Haemophilus influenzae</i>, a naturally competent bacterium, would respond when placed in starvation conditions and what this says about the regulation of the competence response.<br/><br/>
You will need to weave the RNAseq data analysis into your project report. So, how does this data support your wet lab work, and the BLAST/PHASTER practical in data workshop 2? How does it fit together? Can you check the wet lab results in the RNAseq data?<br/>
There is also the opportunity to use the other datasets available to you to further explore these questions, and to explore the data more too.<br/><br/>

One outstanding question is why the C replicates were so different from the rest (look back at your PCA plots). I mentioned that we might be able to control for this difference. Let's take a look.
<br/>
</p>

```R 

# revisit the DEA section - we're going to run another type of DEA
# we can include the replicate information (which is in sample_info) to control for variance based on the replicate
# when we instruct the DEA design, put variable you want to control before the experimental variable you want to measure
dds2 <- DESeqDataSetFromMatrix(countData = comp_counts, colData = sample_info, design = ~ replicate + condition)

# you can then follow the rest of the analysis we did, include summaries and making the volcano plot

```
<p align="justify">
Often we control for replicate if there is a genuine reason to do so, such as replicates coming from different labs, isolates or donor backgrounds. In this case, "C" is not really a separate grouping and we don't know for certain it was treated differently - this may simply reflect genuine biological variance (rather than a technical annoyance). We do get a few more significant hits, but we already had more significant hits than we could feasibly validate in the lab, so it's unlikely we've gained anything properly new.<br/><br/>
However, we can flip the DESeq2 design to work out what it was about the replicates that caused the variance on the PCA.<br/>
If you try this, you'll see that lots of ribosomal RNA comes up as significant. Remember, our dataset was derived after rRNA depletion, so any rRNA leftover is simply a measure of how (in)efficient the library prep was. It looks like the C samples had a worse prep, and therefore have more rRNA leftover, hence the variance on the plot. Again, we have lots of significant results despite this, but a solution could be to remove any rRNA genes from the count matrix, removing this problem. There is a good reason to do this. The danger is <b>overdoing it</b>. This is where good data management and note taking comes in.<br/><br/>

Another good question is whether the response we see at <i>t</i>=30 happens earlier (MIV1; <i>t</i>=10) and/or is maintained (MIV3; <i>t</i>=100). Things you now have the skills to answer. <br/><br/>
Could you make sensible hypotheses about the other datasets you have access to, that may also help your understanding of the Hi competence response under stress? I look forward to seeing what you might come up with...
</p>

<br/><br/><br/><br/><br/><br/>
[Return to tutorial introduction.](#introduction)<br/>
<br/><br/>

### Video explainers
#### Introduction to transcriptomics
<p align="justify">
This video was included in the week 4 teaching materials on the VLE to view before the workshop - hopefully this is just here as a reminder. Here we introduce you to some of the concepts of transcriptomics - where the data comes from, what it looks like, and what you'll be doing with it in the data workshops.<br/>
</p>
<iframe src="https://york.cloud.panopto.eu/Panopto/Pages/Embed.aspx?id=ff92dbc6-c7a3-47c1-96ac-b1260104d15f&autoplay=false&offerviewer=true&showtitle=true&showbrand=true&captions=true&interactivity=all" height="405" width="100%" style="border: 1px solid #464646;" allowfullscreen allow="autoplay" aria-label="Panopto Embedded Video Player"></iframe>
<br/>
[Return to tutorial introduction.](#introduction)<br/>
[Begin the workshop.](#the-workshop)<br/><br/>

#### Data processing before this workshop
<p align="justify">
This video was included in the week 5 teaching materials on the VLE to support the workshop material. We introduce you to the data you will analyse in these two data workshops (3 and 4), and tell you what has happened to go from raw sequencing data (files for each sample total c.40 million lines) to the counts matrix you will work on (1745 rows and 6 columns, in the first instance).<br/>
</p>
<iframe src="https://york.cloud.panopto.eu/Panopto/Pages/Embed.aspx?id=117ba13f-1f0c-4c2a-8f48-b12a00929d38&autoplay=false&offerviewer=true&showtitle=true&showbrand=true&captions=false&interactivity=all" height="405" width="100%" style="border: 1px solid #464646;" allowfullscreen allow="autoplay" aria-label="Panopto Embedded Video Player"></iframe>
<br/>
[Return to tutorial introduction.](#introduction)<br/>
[Begin the workshop.](#the-workshop)<br/><br/>
