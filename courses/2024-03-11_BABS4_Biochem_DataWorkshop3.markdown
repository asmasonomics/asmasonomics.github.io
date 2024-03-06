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
	<li><code>tidyverse</code> <a href="https://doi.org/10.21105/joss.01686">Wickham <i>et al.,</i> 2019</a>
	<li><code>DESeq2</code> <a href="https://doi.org/doi:10.18129/B9.bioc.DESeq2">Love <i>et al.,</i> 2014</a>
	<li><code>EnhancedVolcano</code> <a href="https://github.com/kevinblighe/EnhancedVolcano">Blighe <i>et al.,</i> 2018</a>
	<li><code>ggplot2</code> <a href="https://ggplot2-book.org/">Wickham <i>et al.,</i> 2016</a>
	<li><code>dplyr</code> <a href="https://doi.org/10.21105/joss.01686">Wickham <i>et al.,</i> 2019</a> (this is part of tidyverse, but I always find loading it explicably helps with redundancy issues later on)
</ul><br/>
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
   Most features have a "low" read count, less than ~10000. But a few features have very high read counts.<br/>
   So how can we make the plot easier to interpret?
</details><br/>
</p>

```R

# try a log10 transform of the same column
counts |> ggplot(aes(x = log10(kw20.BHI1.F + 1))) + geom_histogram()

# we do +1 inside the transform because log10(1) is 0. This means the graph is not disrupted by large minus numbers (representing 0.0001 values etc) and missing 0 values.
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
colnames(gam_counts) <- "gam"

# plot
ggplot(gam_counts, aes(x = rownames(gam_counts), y = gam)) + 
  geom_bar(stat = "identity") + 
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = "samples", y = "Hi gam total reads")


```






















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

























<p align="justify">
Welcome to Workshop 4! This workshop will build on some of the skills you developed in previous courses using gene expression data. In previous courses you have worked with <b>count</b> data to look at gene expression, and to explore genes which are differentially expressed between conditions (such as in health and disease). This data has already been processed from raw reads: performing quality control, aligning to a reference genome or transcriptome, and then summarising to gene level. In this workshop you will do the whole process, from raw reads to differential expression analysis.<br/>
If you choose to base your final report on <i>this</i> workshop, you will need to use the hypoxia exposed urothelium dataset and  <b>not the BKPyV dataset used in this workshop </b>. The bioinformatic approach will be very similar, but you will need to bring in relevant biology. More details on this dataset are at the end of the workshop material.<br/><br/>
As ever, the workshop is aimed towards those with Windows machines, including the managed university machines. If you are using your own machine, you have a little more freedom. If you have a Mac, you should be aware of the differences by now (Terminal rather than PowerShell, <i>etc.</i>).<br/>
</p>

### Before you get stuck in
<p align="justify">
As part of this workshop you will use the teaching0.york.ac.uk server, like last time. You will use linux modules and command line R. This should feel nice and familiar! This session will generate plots at the command line which you can then download. If you do RNAseq for your final report, you may want to take your final data off the server and explore with RStudio - this is absolutely fine! Always remember there are lots of different ways to compete the same task when you're coding.<br/>
As soon as you can, get logged on to the server and get the following R libraries installed. You will only need to do this once.<br/>
</p>

Open Windows PowerShell and log on to the server as before:

```sh
ssh USERID@teaching0.york.ac.uk
```

First, create a new shortcut to help with very long path names:

```sh
ln -s /shared/biology/bioldata1/bl-00087h/ genomics

# this means when you log in you can now do this to get to the genomics space
cd ~/genomics/

# you can also use ~/genomics in paths, rather than /shared/biology/bioldatat1/bl-00087h/
# I will use this throughout the workshop material
```

Change directory into your user space within the ~/genomics/students/ area, and make a new directory for this workshop, for example: workshop4. Change into that new directory. We are now ready to set up your R space.

```sh
# create directories for your newly installed R libraries to go (this helps with version control and installing on a managed machine)
mkdir ~/Rlibs ~/Rlibs/R_4.1.2

# type R to launch R from the linux terminal
R 
```

Within R, you can run the following lines one by one, or copy and paste in one go. If any package installations ask for updates, you can skip. Whilst these are installing, read the introductory material for the workshop below.

```R

install.packages("ggrepel", lib="~/Rlibs/R_4.1.2")
install.packages("BiocManager", lib="~/Rlibs/R_4.1.2")
library("BiocManager", lib.loc="~/Rlibs/R_4.1.2")

BiocManager::install("tximport", lib="~/Rlibs/R_4.1.2")
BiocManager::install("rhdf5", lib="~/Rlibs/R_4.1.2")

install.packages("remotes", lib="~/Rlibs/R_4.1.2")
library("remotes", lib.loc="~/Rlibs/R_4.1.2")
remotes::install_github("pachterlab/sleuth", lib="~/Rlibs/R_4.1.2")
remotes::install_github("kevinblighe/EnhancedVolcano", lib="~/Rlibs/R_4.1.2")
remotes::install_github("ctlab/fgsea", lib="~/Rlibs/R_4.1.2")

q()

```

The q() command exits R and brings you back ready for the workshop. <br/><br/>


### Introduction to the material and research question in this workshop
<p align="justify">
Bladder cancer is the 10th most common global cancer and is one of the most expensive to treat as most cancers recur and when the disease progresses, standard practice is to remove the bladder. Despite such radical surgery, 5-year survival is <50% in muscle-invasive disease.<br/>
There is therefore an unmet need to understand how bladder cancers start, and what can be done to prevent them. The long-standing risk factor for bladder cancer is smoking, but mutational signatures (see lecture 7 if you're unsure what these are) showed that the classic mutational profile of smoking seen in lung cancer is <b>not</b> present in bladder cancer. Instead, there are signatures of APOBEC mutagenesis - a family of enzymes which defend against viruses. These signatures are very prevalent in HPV-driven cervical cancers. However, there is no obvious viral cause for bladder cancer, and viral genomes are not found within bladder cancer genomes (as they are with HPV positive cervical cancer).<br/>
</p>
<br/>
![Bladder cancer mutational signature showing APOBEC mutations, but no smoking signature - SBS4](/assets/images/BLCA_SBS_mutational_signature.png){:class="img-responsive"}
<br/><span style="font-size:0.8em;">*Deconvolution of mutational signatures. Top left plot shows the original proportion of single base substitutions in the sample. The right hand side shows the deconvolution into separate SBS derivations and their relevant proportions, and as proof the bottom left plot shows the reconstruction of the signature using the deconvoluted plots.*</span><br/>

<p align="justify">
Based on epidemiolgical data, and high incidence of bladder cancer in kidney transplant patients, researchers at York hypothesised that BK Polyomavirus (BKPyV) may be the cause. The RNAseq data in this workshop was our first effort to explore this association and formed part of a <a href="https://doi.org/10.1038/s41388-022-02235-8">publication in <i>Oncogene</i></a> in 2022.<br/><br/>
</p>

### The workshop
<p align="justify">
Hopefully your R libraries have finished installing, if not it hopefully won't be too much longer.<br/>
All the workshop material is in: <br/>
</p>

```sh
~/genomics/rnaseq_data/01_workshop4_BKPyV-infection/
```

<p align="justify">
As before, please don't copy this data. Either make symbolic links in your Workshop4 folder in your own student area, or just use the paths to the original data. I will show the latter. And don't worry, you don't have permissions to delete the data of other people.<br/><br/>
For this workshop you will work with paired end RNA sequencing data. This means each sample has both a <b>read1.fq.gz</b> and a <b>read2.fq.gz</b> file. Data were derived from cell cultures of biomimetic human urothelium (the epithelial lining of the bladder) which were either infected with BKPyV or not. Cells originated from three different people. Cells were expanded in the lab, split into two dishes where one was infected and the other wasn't. This experimental design allows us to control for the different anti-viral response seen among different people.<br/><br/>
</p>

#### Workshop Aims
1. Run FastQC on a subset of the data, learning how to interpret the output files
2. Align a subset of the data to the human transcriptome, understanding where different gene expression metrics come from
3. Perform diffential expression analysis and gene set enrichment analysis on the entire dataset, exploring the actual biology behind the data

<br/><br/>

### 1 Quality control of raw RNA sequencing data with FastQC
<p align="justify">
Raw sequencing data is in FASTQ format. For each sequencing read (commonly 75, 100 or 150 base pairs in length), there is a name, the As, Ts, Cs and Gs of sequence, and empty line (kept as a "just in case" by the format inventors) and a quality score. When we QC our raw data, it is these quality scores that we are largely assessing - how much do we trust the data we got from the sequencing machine?<br/>
We won't spend more time on the format here, but I have produced a video with <a href="https://elixiruknode.org/">Elixir-UK</a> which I have included below.<br/>
</p>
{% include youtube.html id="tO2H3zuBouw" %} <br/>

Now to run FastQC yourself. Get back to the directory you created for this workshop, and run the following:

```sh
# make a subdirectory for this part of the workshop
mkdir 1_fastqc_test
cd 1_fastqc_test
 
# load the required module
module load bio/FastQC/0.11.9-Java-11

# run on one reduced file (1 million reads only) from the infected and uninfected datasets
fastqc -o ./ ~/genomics/rnaseq_data/01_workshop4_BKPyV-infection/01_test_files_for_fastqc/BKPyV-fastqctest_read1.fq.gz
fastqc -o ./ ~/genomics/rnaseq_data/01_workshop4_BKPyV-infection/01_test_files_for_fastqc/uninfected-fastqctest_read1.fq.gz
```

<p align="justify">
These should finish pretty quickly, creating a zip file (which we don't care about), and an html file. You need to download both of these html files to your local machine, either using sftp at the command line or the WinSCP or FileZilla GUIs (check the Linux II video if you're unsure).<br/>
To use sftp at the command line, open a new PowerShell and do the following (make sure to replace USERid with your user ID):<br/>
</p>

```sh
# or similar path to wherever you want it
chdir C:\Users\USERid\Downloads

sftp USERid@teaching0.york.ac.uk
cd /shared/biology/bioldatat1/bl-00087h/students/USERid/workshop4/1_fastqc_test/
get *html

# this gets you out of sftp
exit

# this closes the PowerShell 
exit
```

<p align="justify">
Open these html files in a browser of your choice, and look at the uninfected first. Immediately you can see that the summary is all ticks (i.e. good!) except one. Take a look through all the plots, particularly the "fail" - what is this graph telling you?<br/><br/>
Now look at the BKPyV infected report. First look at the overrepresented sequence and copy the first one and <a href="https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch">check what it is using BLASTn</a>. What is it?<br/>
Now look at the GC content graphs of both. These are both RNAseq datasets from human, so why do they look different?<br/><br/>
In summary, these data were really high quality (as are all the data we have provided to you). If quality was lacking, we could trim parts of reads, or remove entire reads. Again, I've made some videos with <a href="https://elixiruknode.org/">Elixir-UK</a> on these topics.<br/>
</p><br/>
{% include youtube.html id="0nFwZC6VZyQ" %} <br/>
{% include youtube.html id="wXKxVhOSVa0" %} <br/>
{% include youtube.html id="megMSTmQN7g" %} <br/><br/>

### 2 (Pseudo)alignment of reads to the human transcriptome with kallisto
<p align="justify">
In many cases we align/map our RNAseq reads to the reference genome, but this can be slow and error prone. With a genome as well annotated as human, we can map to the reference transcriptome - a fasta file with every known and predicted transcript from the human genome. This is much much faster and much more accurate.<br/>
In this workshop we will use Gencode v44 protein-coding genes. <a href="https://www.gencodegenes.org/">Gencode</a> is a genome annotation consortium used to inform the <a href="https://www.ensembl.org/index.html">Ensembl genome browser</a>. These are all major tools used by the community. I will cover these concepts in videos in weeks 7, 8 and 9.<br/><br/>
</p>

```sh
# create a directory for your work
cd ../
mkdir 2_kallisto_subset
cd 2_kallisto_subset/

# load the required module
module load bio/kallisto/0.48.0-gompi-2022a

# run kallisto using the index and relevant read files
kallisto quant --index ~/genomics/rnaseq_data/gencode.v44.pc_transcripts-kallisto --output-dir=BKPyVinfected-01 ~/genomics/rnaseq_data/01_workshop4_BKPyV-infection/02_reduced_files_for_workshop_mapping/BKPyVinfected-01-2M_read1.fq.gz ~/genomics/rnaseq_data/01_workshop4_BKPyV-infection/02_reduced_files_for_workshop_mapping/BKPyVinfected-01-2M_read2.fq.gz
kallisto quant --index ~/genomics/rnaseq_data/gencode.v44.pc_transcripts-kallisto --output-dir=uninfected-01 ~/genomics/rnaseq_data/01_workshop4_BKPyV-infection/02_reduced_files_for_workshop_mapping/uninfected-01-2M_read1.fq.gz ~/genomics/rnaseq_data/01_workshop4_BKPyV-infection/02_reduced_files_for_workshop_mapping/uninfected-01-2M_read2.fq.gz
```

<p align="justify">
Each kallisto quant line should take ~2 minutes to run. After you have both - look at the rates (%) of pseudoalignment. How do they differ, and can you think why?<br/><br/>
</p>

### 3 Convert kallisto output to gene-level expression matrix
<p align="justify">
As we have mapped to the transcriptome, we now want to bring our data to the gene-level. This is typical for differential expression analysis, as more is known about gene function rather than that of individual transcripts. This does lose information, particularly if gene expression is regulated using antisense transcripts or function is changed by alternative splicing. That's for another time.<br/><br/>
We now want our data to be in TPMs - transcripts per million. This is a really good metric for RNAseq data as it gives a proportional relationship between transcripts in a population of cells, but because the denominator is big (a million) it is quite robust to changes unless they're real. If you're not sure about counts vs TPMs vs FPKMs and why we use different ones (and shouldn't use the latter anymore) - you can check out another <a href="https://elixiruknode.org/">Elixir-UK</a> video below.<br/>
</p>
{% include youtube.html id="3Pe9xcGF_Wo" %} <br/>

<p align="justify">
Now we're going to combine all the kallisto output, sum to gene-level and start to look at changes. To do this (and to save time) our coding needs to get more complex. Here we will use a <b>for</b> loop.<br/>
The logic of a for loop is to say, for each item in a list, do the same set of commands, until you run out of items in your list. The list could be numbers in a range, lines from a file, or files in a directory (and less commonly the individual letters in a string). This is a really important skill as it saves you typing the same line multiple times (with the chance for typos) for different samples (like the kallisto quant lines from above).<br/>
</p>

The standard format for a for loop is:
```sh
for item in list
  do
    command1
    command2
done
```
<p align="justify">
The <b>for</b> is a keyword which tells the command line that a loop statement is coming. The <b>do</b> keyword indicates the start of the loop, and the <b>done</b> shows the end. This means you can hit enter after each line and the terminal waits to execute (like with a backslash). We will use this now:<br/>
</p>

```sh
for dir in *
  do
    cp ${dir}/abundance.tsv ${dir}_abundance.tsv
done
```

<p align="justify">
Now it's time to go into R, using the libraries you installed before. We will use an R package called tximport to combine the different datasets together and summarise to the gene level.<br/>
</p>

```sh
# check what an abundance.tsv file looks like
head BKPyVinfected-01_abundance.tsv

# take the first transcript ID, missing the version (.2) at the end and check for its gene name
grep "ENST00000641515" ~/genomics/rnaseq_data/gencode.v44.pc_transcripts.t2g

# now let's get into R
R
```
```R
library("tximport", lib.loc="~/Rlibs/R_4.1.2")

# load the transcript to gene (t2g) file
t2g <- read.csv("~/genomics/rnaseq_data/gencode.v44.pc_transcripts.t2g", sep="\t", header=TRUE)

# create a list of abundance files
files <- list.files(".","tsv$")
names(files) <- gsub("_abundance.tsv", "", files)

# use tximport to combine the data and write gene-level TPMs to file
txi <- tximport(files, type = "kallisto", tx2gene = t2g, ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)
txiDF <- as.data.frame(txi$abundance)
txiDF <- cbind(genes = rownames(txiDF), txiDF)
rownames(txiDF) <- NULL
write.table(txiDF,file="allTPMs.tsv",sep="\t", quote=FALSE, row.names=FALSE)

#leave R 
q()
```

<p align="justify">
Now we have gene-level TPM values for our samples - nice work! This is approximately where you would have started when you have done gene expression work before, but now you have experience of data QC, alignment and quantification.<br/><br/>
Think back to our hypothesis for these data. We are trying to work out if BKPyV can infect urothelium (yep, tick), what the impact is on the transcriptome (in progress), and whether this is consistent with BKPyV causing bladder cancer (can we answer this with this experiment?).<br/>
Often with RNAseq data you have a few key indicator genes in mind which could let you know if your experiment has worked. Before doing differential expression it is always nice to check these and get a feel for the data. This is the first time we can check the data is biologically meaningful, not just that the data we got off the machine is technically fine (our initial QC). The data could be of good technical quality, but only have sampled <i>GAPDH</i> millions of times, for example.<br/>
We will now manually check a few genes: <i>APOBEC3A</i> and <i>APOBEC3B</i> (viral response genes thought to cause the mutational signatures we see in bladder cancer), <i>MKI67</i> (marker of active proliferation, when bladder urothelium is typically quiescent; out of cell cycle) and <i>KRT13</i> (a marker of urothelial differentiation).<br/> 
</p>

```sh 
head -n 1 allTPMs.tsv; egrep "APOBEC3A|KRT13|MKI67" allTPMs.tsv
```

<p align="justify">
What does the (limited) data suggest to you?<br/>
</p><br/>

### 4 Differential Expression Analysis (DEA) with Sleuth
<p align="justify">
Taking a quick look at TPMs is one thing, but we want to be unbiased and use the full statistical power of the dataset. For that we need DEA (done in R) using the full kallisto folders similar to those you made in part 2 (not the abundance.tsv files) and an input file showing our experimental design, which we will make now. This is another for loop, with each line doing text manipulation of the experiment metadata stored in our read names - this is why consistent and informative naming is so useful.<br/><br/>
A quick confession. You will now use some kallisto output which I created from the full dataset (not just the first 2 million reads as you did). This is because it would take more like 20-25 minutes per sample to run (rather than 2 minutes). If you do the RNAseq data for your report, you should run kallisto in full for all samples (I will give the full command at the end of the workshop material below). Part of this is we use a statistical method called bootstrapping to give us confidence in our expression values. Bootstrapping is where we re-run the data alignment multiple times to get consistent answers (like you would when building phylogenetic trees). In the full run you will do 20 bootstraps compared to the 0 you did in part 2. This explains why you <i>might</i> see (very) slight differences between your TPMs from part 3 and those of your classmates around you.<br/>
</p>

```sh
cd ../
mkdir 3_DEA
cd 3_DEA

# create the experimental design file for sleuth
# first write the headers
echo -e "sample\texpress\tbiorep\tpath" > infectionDEA.info

# now loop through the kallisto output directory filenames to get the information needed
for sampledir in ~/genomics/rnaseq_data/01_workshop4_BKPyV-infection/03_full_kallisto_output_for_workshop_DEA/*/abundance.h5
  do
    sample=`echo $sampledir | rev | cut -d'/' -f2 | rev`
    donor=`echo $sample | cut -d'-' -f2`
    exp=`echo $sample | awk -v sample=$sample '{if (sample~/^BKPyV/) {print "0"} else {print "1"}}'`
    echo -e $sample'\t'$exp'\t'$donor'\t'$sampledir >> infectionDEA.info
done
```

Now we can use this and head back into R.
```sh 
R
```
```R
library(tidyverse)
library(readr)
library("sleuth", lib.loc="~/Rlibs/R_4.1.2")

# create an object for the experimental design and t2g
s2c <- read.table("infectionDEA.info", header=TRUE, stringsAsFactors=FALSE)
t2g <- dplyr::select(read.table("~/genomics/rnaseq_data/gencode.v44.pc_transcripts.t2g", header=TRUE, stringsAsFactors=FALSE), target_id = ensembl_transcript_id, ext_gene = external_gene_name)

# build a sleuth object, aggregating to gene level
so <- sleuth_prep(s2c, ~express + biorep, extra_bootstrap_summary=TRUE, num_cores=1, target_mapping=t2g, transformation_function = function(x) log2(x+1.5), gene_mode=TRUE, aggregation_column = 'ext_gene')

# fit the statistical models
# start with the 'full' model accounting for variation due to experimental condition AND donor background
# then account for variance just due to difference between donors
# perform a likelihood ratio test (LRT) to identify significant changes just associated with the experimental condition
so <- sleuth_fit(so, ~express + biorep, 'full')
so <- sleuth_fit(so, ~biorep, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table <- sleuth_results(so, 'reduced:full', test_type='lrt', show_all = FALSE)

# write results to file
write.table(distinct(as.data.frame(sleuth_table)[,c("target_id","pval","qval")]),file="infectionDEA_results.tsv",sep="\t",row.names=FALSE,col.names=TRUE)

q()
```

<p align="justify">
As we are running statistical tests for so many genes, we need to have a p value correction. Sleuth uses Benjamini-Hochberg, a ranking system to reduce false positives. We then use these 'q' values, again with a stringent threshold such as <0.05. Use head, less, grep or other tools you know to look at the significant genes - do any come up that you recognise? Do any of our "quick check" genes come up? If not, why do you think this is?<br/><br/>
The other thing to consider here is that q<0.05 genes are <i>statistically</i> different, but that doesn't necessarily equal biologically significant - the fold change in expression could be very small, but appear very significant just due to the power of the test (i.e. the ability of the test to see changes). This is particularly common when you have replicates from cell lines (as these are more like technical repeats than biological ones).<br/>
Now we're going to combine our TPMs and our sleuth results and make a volcano plot.<br/>
</p>

```sh 
# create symbolic link for TPMs
ln -s ~/genomics/rnaseq_data/01_workshop4_BKPyV-infection/03_full_kallisto_output_for_workshop_DEA/allTPMs.tsv allTPMs.tsv

# now get back into R
R
```

```R
library(tidyverse)
library(dplyr)

# load TPMs and look at it
tpms <- read.table("allTPMs.tsv", header=TRUE)
head(tpms)

# add columns for the average condition TPMs for each gene
tpms <- mutate(tpms, BKPyVinfected_avg = rowMeans(select(tpms, starts_with("BKPyVinfected"))))
tpms <- mutate(tpms, uninfected_avg = rowMeans(select(tpms, starts_with("uninfected"))))

# calculate the log2(fold change) so that positive values are genes with higher expression in BKPyV
# we do a "+1" so that genes with very low expression don't get big fold changes
tpms <- mutate(tpms, log2FC = log2((BKPyVinfected_avg+1)/(uninfected_avg+1)))

# round all value columns to 2 decimal places to make it more manageable(!) and take a look
tpms <- mutate(tpms, across(2:ncol(tpms), round, 2))
head(tpms)

# now load the differential expression analysis results and take a look
dea <- read.table("infectionDEA_results.tsv", header=TRUE)
head(dea)

# now we combine our TPMs and our DEA results using a left join on the gene names
# in our TPMs dataframe our gene names are in a column called "genes"
# in our DEA dataframe gene names are in a column called "target_id"
res <- tpms %>% left_join(dea, by=c("genes"="target_id"))
head(res)

# if sleuth cannot do a stat comparison, it doesn't produce a test value so the p and q values are NA after the join
# let's replace NA values with 1 (i.e. not significant)
res[is.na(res)] <- 1

# we can create a combined metric of fold change and significance called a pi value
res <- mutate(res, pi = (log2FC * (-1 * log10(qval))))

# now we have all our data in once table, let's save it
write.table(res, file="TPMs_and_DEA_results_file.tsv", sep="\t", row.names=FALSE, col.names=TRUE)

# don't close R!
```

<p align="justify">
OK! That was a great block of R coding, and you're almost (almost) at the plotting stage.<br/><br/>
As a recap, you took the full kallisto output which allowed you to get TPMs for each gene - a metric for gene expression. You then used the kallisto output to run sleuth - allowing you to use proper statistics to see which genes have significantly changing expression due to infection with BKPyV. Then you took the output from these two steps and merged them using an inner join. This means you have expression values and significance values in the same place.<br/><br/>
Now you can plot!
<br/>
</p>

```R 
# carrying on with our R session - if you accidentally closed it, just reload TPMs_and_DEA_results_file.tsv using read.table()
library("ggrepel", lib.loc="~/Rlibs/R_4.1.2")
library("EnhancedVolcano", lib.loc="~/Rlibs/R_4.1.2")

# check how many significantly different genes we have, using log2FC threshold of >=0.58 (50% increase) or <=-0.58 (50% decrease) and qvalue of <0.05
sum(res$log2FC>=0.58 & res$qval<0.05)
sum(res$log2FC<=-0.58 & res$qval<0.05)

# create list of most significantly different genes for labelling the volcano - change the 100 to get more/fewer genes
siggenes <- (res %>% arrange(desc(abs(pi))) %>% slice(1:100))$genes

# let's create our volcano plot
EnhancedVolcano(res, x = "log2FC", y = "qval", FCcutoff = 0.58, pCutoff = 0.05,
	lab = res$genes, selectLab = siggenes, labSize = 2.0, max.overlaps = 1000, drawConnectors = TRUE,
	ylim = c(0, (max(-1*log10(res$qval)) + 0.5)), xlim = c((max(abs(res$log2FC))*-1), max(abs(res$log2FC))+1),
	legendPosition = 0, gridlines.major = FALSE, gridlines.minor = FALSE)

# save it
# you can't view it on the server, but use sftp or WinSCP/FileZilla to download and view
ggsave("volcano.pdf")

# don't close R!
```

<p align="justify">
You now have a volcano plot - the hallmark of an RNAseq experiment! In a volcano you plot log2 of the fold change (x axis) against the -log10 of the significance value (y axis). These two log transformations are really important for actually understanding the data.<br/><br/>
We use log2FC so that up and down changes are symmetrical around zero. If your expression goes from 10 to 20, this is a fold change of 2, but if you go from 10 to 5 this is a fold change of 0.5. Whilst the fold change should be symmetrical, 0.5 is closer to 1 than 1 is to 2. This makes your plots wonky, and makes the up changes seem more important. If you do a log2 transform, a fold change of 1 (i.e. no change) gives a log2FC of 0. log2 of 0.5 is equal to -1 and log2 of 2 is equal to 1. Now we have a symmetrical plot.<br/><br/>
We want to have small p/q values to give us confidence differences between conditions are not just due to chance. But, we want these important points to be highlighted at the top of our graph. If we plotted fold change against p or q values, all the important stuff would be on or very near the x axis (squished at the bottom) and all the non-significant stuff would be at the top. So, we do -log10 of the stat test value, and this puts the rubbish at the bottom and the (hopefully) interesting genes at the top.<br/><br/>
One important note. We did include a +1 transformation when calculating fold change to get rid of seemingly huge changes between decimals (the fold change between 0.000005 and 0.0001 is huge!) but some of your most significant changes may come from very low expression. Is this biologically meaningful? Think about times when genes with low expression can still be incredibly important for cell function.<br/><br/>
Now, our comparison doesn't have that many massively significant changes. One step is to go through the genes which have changed and try and see patterns - do they make sense? A more powerful and unbiased technique is to use gene set enrichment analysis (GSEA). Instead of just focusing on significant changes, GSEA looks for patterns in the dataset as a whole. We'll do this now.
<br/>
</p>

```R 
library("fgsea", lib.loc="~/Rlibs/R_4.1.2")

# use our pi values to rank the genes biologically from most up to most down
prerank <- res[c("genes", "pi")]
prerank <- setNames(prerank$pi, prerank$genes)
str(prerank)

# run GSEA using a list of gene sets curated by MSigDB (Broad Institute)
genesets = gmtPathways("~/genomics/rnaseq_data/h.all.v2023.2.Hs.symbols.gmt")
fgseaRes <- fgsea(pathways = genesets, stats = prerank, minSize=15, maxSize=500)

# check the top most significant hits
head(fgseaRes[order(pval), ], 10)

# plot an enrichment plot for the top hit, and explore some of the others
# google the terms, or check the MSigDB website - do the results make sense?
plotEnrichment(genesets[["HALLMARK_E2F_TARGETS"]], prerank) + labs(title="E2F targets")
ggsave("E2Ftargets.pdf")

# save your results
write.table(res, file="GSEA_results_file.tsv", sep="\t", row.names=FALSE, col.names=TRUE)

# when you're happy, the workshop is over - you can quit!
q()
```
<br/>
<p align="justify">
Excellent work! You have now taken a dataset from raw reads, to gene expression values, to differential expression and then to investigating the biology to answer  pertinent research question. This is genomics and bioinformatics in action!
<br/><br/>
</p>

### Concluding remarks
<p align="justify">
In terms of the biology of this dataset, you have seen that infection with BKPyV causes the urothelium to alter its transcriptome (significant gene changes). These changes are to do with the antiviral response (interferon gamma and alpha responses in GSEA) and changes to the cell cycle. Normally urothelium is arrested in G0, with only about 1% of the cells in cycle at any time. The virus causes the cells to re-enter but not complete the cell cycle, instead sitting at the G2M checkpoint where BKPyV is best placed to replicate.<br/>
The APOBEC response is there, but this is very donor-dependent (resting APOBEC3A levels vary a lot between donors) and therefore does not come up with a significant q value. But, we did see that the urothelium does actually become infected and that it can induce changes to the cell cycle and DNA replication. APOBEC responds to the presence of viral DNA/RNA. We have continued with the work to see if BKPyV infections led to the damage in DNA you would associate with the APOBEC damage in tumours - very much work still in development.<br/><br/>
If you're interested in reading more, check out our <a href="https://doi.org/10.1038/s41388-022-02235-8">publication in <i>Oncogene</i></a>. Also, look out for MBiol projects in this area next year.
<br/><br/>
</p>

### What to do if you want to do RNAseq for your report
<p align="justify">
Hopefully you've found this session useful, interesting and inspirational (plus the other learning objectives...). If you would like to do RNAseq analysis for your Genomics3 assessment, you <b>should not</b> use the data from above. Instead I have provided with a similar dataset of human urothelial cell cultures grown in either normoxic (20% O2) or hypoxic (1% O2) conditions:<br/>
</p>

```sh
~/genomics/rnaseq_data/02_hypoxia_assessment_dataset/
```

<p align="justify">
This is a 4vs4 dataset where urothelial cells from 4 different people have been used for cell culture. Cells were re-differentiated to form a biomimetic tissue and then split, with one half continuing to be grown in standard normoxic conditions (20% O2) and the others grown in hypoxia (1% O2). We are able to measure the ability of biomimetic urothelium to form a tight barrier (to urine, if it was in the body), and the barrier ability was highly compromised when the cells were taken into hypoxia.<br/>
In healthy people, the urothelium is an epithelial layer some 3-6 cells thick (depending on how full the bladder is) which sits on top of a basement membrane associated with a capillary bed (i.e. it has good oxygen supply). For this report consider what reductions in oxygen would mean for a tissue, its physiology and the resultant transcriptome. Consider how this could impact bladder health, or in which bladder diseases hypoxia may play a role. Are there major regulators of hypoxia which do not change at the transcript level?<br/><br/>
For your analysis, run through the same process as in this workshop: fastQC for looking at the reads, kallisto for the alignment, tximport in R for getting the gene-level TPM files, sleuth in R for running differential expression analysis, and then interpretation by plotting individual genes, making a volcano plot and using GSEA to report a little on the biology of the dataset. Look at the very bottom of this page for a couple of extra commands to help you.<br/>
There is more information on the assessment criteria on the VLE. Essentially we are looking for a decent introduction to the topic and relevance of the dataset, an explanation of your methods, and then your attempt to interpret the results and suggest how they address the purpose of the study.<br/>
Good luck, remember to support each other using the discussion boards, and do ask us for help and guidance if you need it.<br/><br/>
</p>

#### Some extra commands to help
<p align="justify">
Remember, you should run fastQC on each read file (all 16). You don't need to put any of the graphs into your report, but it is always good practice to check the quality of your data, and putting information on average read number is often good practice.<br/>
You then need to run kallisto on each sample (all 8) using both read 1 and read 2. You will also include bootstrapping here (which you didn't do in the workshop) to make your DEA more robust. The command for that is below. Then your Sleuth, plotting and GSEA commands should be very similar, just changing the read names and any relevant paths to data.<br/>
fastQC on these full files will take 10-20 minutes per file. kallisto will take 40-60 minutes per sample. Use loops, and set the commands running in the background, or using screen. Remember, you can set jobs running, and then go off and do something else while they run. Setting something running at 5pm often means it is ready for you in the morning.<br/>
</p>

```sh 
# a fastqc loop
# this would loop through any file ending in gz in the current directory and run fastqc one after the other
for readfile in *gz
  do
    fastqc -o ./ $readfile
done

# this fastqc loop would submit each fastqc job into the background 
for readfile in *gz
  do
    fastqc -o ./ $readfile &
done

# this is the command for kallisto with the bootstraps
# you could run each of the 8 samples with an '&' at the end to run in the background
# you could put this in a loop too, if you're careful about how you define each file 
# remember, sometimes it is quicker (if you're doing something once) to type it out vs spending an hour getting a loop right...
kallisto quant --index ~/genomics/rnaseq_data/gencode.v44.pc_transcripts-kallisto --output-dir=SAMPLENAME --bootstrap-samples=20 ./SAMPLENAME_read1.fq.gz ./SAMPLENAME_read2.fq.gz 
```

