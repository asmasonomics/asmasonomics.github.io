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
Hopefully this all sounds quite familiar! If not, after this morning, make time to watch the videos we've generated to support these two workshops. These are embedded below and on the VLE.<br/>
</p>
[Introduction to transcriptomics](#introduction-to-transcriptomics)<br/>
[How has the data in this workshop been processed so far?](#data-processing-before-this-workshop)
<br/>

### The workshop
#### Set up your RStudio project
<p align="justify">
Start RStudio from the Start menu.<br/><br/>
Make a new RStudio project. Put it in a sensible place with a sensible name. Remember that if you are using a university machine, select a location via the M/H drives, not just using the "Documents" shortcut. This <b>will</b> create problems for you!<br/><br/>
Use the Files pane to make subdirectories for your <code>raw_data</code>, <code>proc_data</code> and <code>plots</code>. These names are suggestions only. Make a new script file, perhaps <code>data_workshop_3.R</code>, to complete your work.<br/><br/>
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
library(EnhancedVolcano)
library(ggplot2)
library(dplyr)

```
<br/>

#### Load datasets
```R

# counts data
counts <- read.table("raw_data/Hi_PRJNA293882_counts.tsv", row.names = 1, header = TRUE)

# feature IDs and symbols where possible
featname <- read_tsv("raw_data/Hi_feature_names.tsv")

# feature locations, setting the column names for BED format
featlocs <- read_tsv("raw_data/Hi_feature_locations.bed", col_names = c("chr","start","end","feat_ID","biotype","strand"))

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
mu_feats <- featname[grep("muA|muB|gam", featname$symbol),]$feat_id
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



















<p align="justify">
<br/>
</p>


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
