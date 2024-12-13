---
layout: page
permalink: /courses/Genomics3_Workshop4_RNAseq_Nov2024
---

![Genomics3 banner](/assets/images/genomics_banner.jpeg){:class="img-responsive"}

<span style="font-size:1.6em;">**Genomics 3 - Workshop 4: RNAseq**</span><br/>

<p align="justify">
Welcome to Workshop 4! This workshop will build on some of the skills you developed in previous courses using gene expression data. In previous courses you have worked with <b>count</b> data to look at gene expression, and to explore genes which are differentially expressed between conditions (such as in health and disease). This data has already been processed from raw reads: performing quality control, aligning to a reference genome or transcriptome, and then summarising to gene level. In this workshop you will do the whole process, from raw reads to differential expression analysis.<br/>
If you choose to base your final report on <i>this</i> workshop, you will need to use the hypoxia exposed urothelium dataset and  <b>not the BKPyV dataset used in this workshop </b>. The bioinformatic approach will be very similar, but you will need to bring in relevant biology. More details on this dataset are at the end of the workshop material.<br/><br/>
As ever, the workshop is aimed towards the dual-boot Linux machines in G/N/169. You will work at the Linux command line, including working in <code>R</code> from the command line.<br/><br/>
</p>

### Introduction to the material and research question in this workshop
<p align="justify">
Bladder cancer is the 10th most common global cancer and is one of the most expensive to treat as most cancers recur and when the disease progresses, standard practice is to remove the bladder. Despite such radical surgery, 5-year survival is <50% in muscle-invasive disease.<br/>
There is therefore an unmet need to understand how bladder cancers start, and what can be done to prevent them. The long-standing risk factor for bladder cancer is smoking, but mutational signatures (see lecture 7 if you're unsure what these are) showed that the classic mutational profile of smoking seen in lung cancer is <b>not</b> present in bladder cancer. Instead, there are signatures of APOBEC mutagenesis - a family of enzymes which defend against viruses. These signatures are very prevalent in HPV-driven cervical cancers. However, there is no obvious viral cause for bladder cancer, and viral genomes are not found within bladder cancer genomes (as they are with HPV positive cervical cancer).<br/>
</p>
<br/>
![Bladder cancer mutational signature showing APOBEC mutations, but no smoking signature - SBS4](/assets/images/BLCA_SBS_mutational_signature.png){:class="img-responsive"}
<br/><span style="font-size:0.8em;">*Deconvolution of mutational signatures. Top left plot shows the original proportion of single base substitutions in the sample. The right hand side shows the deconvolution into separate SBS derivations and their relevant proportions, and as proof the bottom left plot shows the reconstruction of the signature using the deconvoluted plots.*</span><br/>

<p align="justify">
Based on epidemiolgical data, and high incidence of bladder cancer in kidney transplant patients, researchers at York hypothesised that BK Polyomavirus (BKPyV) may be the cause. The RNAseq data in this workshop was our first effort to explore this association and formed part of a <a href="https://doi.org/10.1038/s41388-022-02235-8">publication in <i>Oncogene</i></a> in 2022. This data is paired end RNA sequencing data. This means each sample has both a <b>read1.fq.gz</b> and a <b>read2.fq.gz</b> file. Data were derived from cell cultures of biomimetic human urothelium (the epithelial lining of the bladder) which were either infected with BKPyV or not. Cells originated from three different people. Cells were expanded in the lab, split into two dishes where one was infected and the other wasn't. This experimental design allows us to control for the different anti-viral response seen among different people.<br/><br/>
</p>

### The workshop
#### Workshop Aims
1. Run FastQC on a subset of the data, learning how to interpret the output files
2. Align a subset of the data to the human transcriptome, understanding where different gene expression metrics come from
3. Perform diffential expression analysis and gene set enrichment analysis on the entire dataset, exploring the actual biology behind the data

<br/><br/>

#### 0 Set up your directories
<p align="justify">
<b>Remember, you should be working on a Linux machine in G/N/169 for this workshop.</b><br/>
Assuming that you are, open up a terminal window, change into your student directory, and set up the directories for this workshop:<br/>
</p>

```sh
# change into your student directory
# this assumes you have made a directory which matches your username exactly - change if needed
cd /shared/biology/bioldata1/bl-00087h/students/${USER}/

# create directories for this workshop, and set their permissions
mkdir workshop4 workshop4/1_fastqc_test workshop4/2_kallisto_subset workshop4/3_DEA
chmod -R 775 workshop4

# enter the workshop 4 directory
cd workshop4/

# create a symbolic link for the rnaseq data to make your paths more manageable
ln -s /shared/biology/bioldata1/bl-00087h/data/rnaseq_data rnaseq_data

```
<br/>
<p align="justify">
<details>
   <summary>What to do if you want to use teaching0 rather than work in G/N/169.</summary>
   <br/>
   The <b>best</b> option is to work in G/N/169. But if this isn't possible you <i>can</i> complete the work using teaching0. I have put some things in place to help, but I can't support this workshop over multiple platforms. I say again - use G/N/169.<br/><br/>
   1) You will need to use the following FastQC and kallisto versions:
   <pre><code class="language-bash">
   module load bio/FastQC/0.11.9-Java-11
   module load bio/kallisto/0.48.0-gompi-2022a
   </code></pre>
   <br/>
   2) You will need to use this index for kallisto:<br/>
   <code>../rnaseq_data/kallist0.48.0_gencode.v44.pc_transcripts.processed-kallisto</code>
   <br/><br/>
   3) When you use <code>R</code> it will be version 4.1.2 (not v4.3.3). You will need to change the <code>.libPaths()</code> commands across the workshop:
   <pre><code class="language-r">
   .libPaths("../rnaseq_data/R_4.1.2")
   </code></pre>
   <br/>
   Everything else should work, but you will need to download any figures you make (or the FastQC html files) using WinSCP/FileZilla/<code>ftp</code>. We have not provided support for this process.<br/>
</details>
</p>
<br/>

#### 1 Quality control of raw RNA sequencing data with FastQC
<p align="justify">
Raw sequencing data is in FASTQ format. For each sequencing read (commonly 75, 100 or 150 base pairs in length), there is a name, the As, Ts, Cs and Gs of sequence, and empty line (kept as a "just in case" by the format inventors) and a quality score. When we QC our raw data, it is these quality scores that we are largely assessing - how much do we trust the data we got from the sequencing machine?<br/>
We won't spend more time on the format here, but I have produced a video with <a href="https://elixiruknode.org/">Elixir-UK</a> which I have included below.<br/>
</p>
{% include youtube.html id="tO2H3zuBouw" %} <br/>

Now to run FastQC yourself.

```sh
# enter the fastqc test directory
cd 1_fastqc_test
 
# load the required module
module load FastQC/0.12.1-Java-11

# run on one reduced file (1 million reads only) from the infected and uninfected datasets
fastqc -o ./ ../rnaseq_data/01_workshop4_BKPyV-infection/01_test_files_for_fastqc/BKPyV-fastqctest_read1.fq.gz
fastqc -o ./ ../rnaseq_data/01_workshop4_BKPyV-infection/01_test_files_for_fastqc/uninfected-fastqctest_read1.fq.gz

```

<p align="justify">
These should finish pretty quickly, creating a zip file (which we don't care about), and an html file showing the quality assessment. You can open this html file using firefox.<br/>
</p>

```sh
# open both files using a wildcard
firefox *html

```

<p align="justify">
Look at the uninfected first. Immediately you can see that the summary is all ticks (i.e. good!) except one. Take a look through all the plots, particularly the "fail" - what is this graph telling you?<br/><br/>
Now look at the BKPyV infected report. First look at the overrepresented sequence and copy the first one and <a href="https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch">check what it is using BLASTn</a>. What is it?<br/>
Now look at the GC content graphs of both. These are both RNAseq datasets from human, so why do they look different?<br/><br/>
In summary, these data were really high quality (as are all the data we have provided to you). If quality was lacking, we could trim parts of reads, or remove entire reads. No trimming is necessary here, so we can move on to pseudoalignment.<br/><br/>
If you're interested to find out more about RNAseq QC and trimming, I've made some videos with <a href="https://elixiruknode.org/">Elixir-UK</a> on these topics.<br/>
</p><br/>
{% include youtube.html id="0nFwZC6VZyQ" %} <br/>
{% include youtube.html id="wXKxVhOSVa0" %} <br/>
{% include youtube.html id="megMSTmQN7g" %} <br/><br/>


#### 2 (Pseudo)alignment of reads to the human transcriptome with kallisto
<p align="justify">
In many cases we align/map our RNAseq reads to the reference genome, but this can be slow and error prone. With a genome as well annotated as human, we can map to the reference transcriptome - a fasta file with every known and predicted transcript from the human genome. This is much much faster and much more accurate.<br/>
In this workshop we will use Gencode v44 protein-coding genes. <a href="https://www.gencodegenes.org/">Gencode</a> is a genome annotation consortium used to inform the <a href="https://www.ensembl.org/index.html">Ensembl genome browser</a>. These are all major tools used by the community. I will cover these concepts in videos in weeks 7, 8 and 9.<br/><br/>
</p>

```sh
# change into directory for part 2 of the workshop
cd ../2_kallisto_subset/

# load the required module
module load kallisto/0.51.1-gompi-2023b

# run kallisto using the index and relevant read files
kallisto quant --index ../rnaseq_data/gencode.v44.pc_transcripts.processed-kallisto --output-dir=BKPyVinfected-01 ../rnaseq_data/01_workshop4_BKPyV-infection/02_reduced_files_for_workshop_mapping/BKPyVinfected-01-2M_read1.fq.gz ../rnaseq_data/01_workshop4_BKPyV-infection/02_reduced_files_for_workshop_mapping/BKPyVinfected-01-2M_read2.fq.gz
kallisto quant --index ../rnaseq_data/gencode.v44.pc_transcripts.processed-kallisto --output-dir=uninfected-01 ../rnaseq_data/01_workshop4_BKPyV-infection/02_reduced_files_for_workshop_mapping/uninfected-01-2M_read1.fq.gz ../rnaseq_data/01_workshop4_BKPyV-infection/02_reduced_files_for_workshop_mapping/uninfected-01-2M_read2.fq.gz

```

<p align="justify">
Each kallisto quant line should take ~2 minutes to run. After you have both - look at the rates (%) of pseudoalignment. How do they differ, and can you think why?<br/><br/>
We have using kallisto to map to all the annotated protein-coding gene transcripts of the human genome. Whilst knowing how individual transcript variants behave can be very interesting, we typically focus our interpretation at the <i>gene</i> level. More is known about gene function rather than that of individual transcripts. This does lose information, particularly if gene expression is regulated using antisense transcripts or function is changed by alternative splicing... But that's for another time.<br/><br/>
We now want our data to be in TPMs - transcripts per million. This is a really good metric for RNAseq data as it gives a proportional relationship between transcripts in a population of cells, but because the denominator is big (a million) it is quite robust to changes unless they're real. If you're not sure about counts vs TPMs vs FPKMs and why we use different ones (and shouldn't use the latter anymore) - you can check out another <a href="https://elixiruknode.org/">Elixir-UK</a> video below.<br/>
</p>
{% include youtube.html id="3Pe9xcGF_Wo" %} <br/>

<p align="justify">
Now we're going to combine all the kallisto output, sum to gene-level and start to look at changes. To do this (and to save time) our coding needs to get more complex. Here we will use a <b>for</b> loop.<br/>
The logic of a for loop is to say, for each item in a list, do the same set of commands, until you run out of items in your list. The list could be numbers in a range, lines from a file, or files in a directory (and less commonly the individual letters in a string). This is a really important skill as it saves you typing the same line multiple times (with the chance for typos) for different samples (like the kallisto quant lines from above).<br/>
</p>

The standard format for a for loop is:
```sh
# don't run this - it is just an example!
for item in list
  do
    command1
    command2
done

```
<p align="justify">
The <b>for</b> is a keyword which tells the command line that a loop statement is coming. The <b>do</b> keyword indicates the start of the loop, and the <b>done</b> shows the end. This means you can hit enter after each line and the terminal waits to execute the entire command block. We will use this now:<br/>
</p>

```sh
# do run this
for dir in *
  do
    cp ${dir}/abundance.tsv ${dir}_abundance.tsv
done

```

<p align="justify">
These abundance files contain expression information for each transcript, but the gene name is not clear. We've copied the abundance data for each sample into the same directory so we can work on them together more easily.
</p>

```sh
# check what an abundance.tsv file looks like
head BKPyVinfected-01_abundance.tsv

# take the first transcript ID, missing the version (.2) at the end and check for its gene name
grep "ENST00000641515" ../rnaseq_data/gencode.v44.pc_transcripts.t2g

```

<p align="justify">
Now it's time to go into R. R is a program like any other and can be run from the terminal by entering the command <code>R</code>. This starts R, just like you would open R or RStudio on a Windows (or Linux...<i>next week...</i>) machine.<br/><br/>
We've done all the package installations you need for this workshop, so you <i>shouldn't</i> need to install any packages. If you're asked to do so, ask us for help before start installing things.<br/>
We will now use an R package called tximport to combine the different datasets together and summarise to the gene level.<br/>
</p>

```sh
# enter R by typing R in the terminal
# make a note of which version you're using
R
```

```R
# you are now using R!

# access the Genomics shared libraries with this command
.libPaths("../rnaseq_data/R_4.3.3")

# load tximport
library(tximport)

# load the transcript to gene (t2g) file
t2g <- read.csv("../rnaseq_data/gencode.v44.pc_transcripts.t2g", sep="\t", header=TRUE)

# create a list of abundance files
files <- list.files(".","abundance.tsv$")

# remove the abundance.tsv from the name to get the sample ID
names(files) <- gsub("_abundance.tsv", "", files)

# use tximport to combine the data and write gene-level TPMs to file
txi <- tximport(files, type = "kallisto", tx2gene = t2g, ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)
txiDF <- as.data.frame(txi$abundance)
txiDF <- cbind(genes = rownames(txiDF), txiDF)
rownames(txiDF) <- NULL
write.table(txiDF,file="allTPMs.tsv",sep="\t", quote=FALSE, row.names=FALSE)

```

<p align="justify">
Now we have gene-level TPM values for our samples - nice work! This is approximately where you would have started when you have done gene expression work before, but now you have experience of data QC, alignment and quantification.<br/><br/>
Think back to our hypothesis for these data. We are trying to work out if BKPyV can infect urothelium (yep, tick), what the impact is on the transcriptome (in progress), and whether this is consistent with BKPyV causing bladder cancer (can we answer this with this experiment?).<br/><br/>
Often with RNAseq data you have a few key indicator genes in mind which could let you know if your experiment has worked. Before doing differential expression it is always nice to check these and get a feel for the data. This is the first time we can check the data is biologically meaningful, not just that the data we got off the machine is technically fine (our initial QC).<br/>
We will now check a few genes: <i>APOBEC3A</i> and <i>APOBEC3B</i> (viral response genes thought to cause the mutational signatures we see in bladder cancer), <i>MKI67</i> (marker of active proliferation, when bladder urothelium is typically quiescent; <i>i.e.</i> out of cell cycle) and <i>KRT13</i> (a marker of urothelial differentiation).<br/> 
</p>

```R
# load libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# write a function which can take a list of multiple genes, and create average and standard deviation value per group, ready to plot
prepare_data_for_plot <- function(data, gene_names) {
  data %>%
    # Filter for the specified genes
    filter(genes %in% gene_names) %>%
    # Reshape from wide to long format
    pivot_longer(cols = -genes, names_to = "Group", values_to = "Value") %>%
    # Extract group prefix based on column name
    mutate(GroupPrefix = case_when(
      startsWith(Group, "BKPyV") ~ "BKPyV",
      startsWith(Group, "uninfected") ~ "uninfected",
      TRUE ~ "Other"  # This handles unexpected prefixes, if any
    )) %>%
    # Calculate mean and standard deviation for each prefix group and gene
    group_by(genes, GroupPrefix) %>%
    summarise(
      Mean = mean(Value),
      SD = sd(Value),
      .groups = "drop"
    )
}

# define the list of genes we want to plot
gene_list <- c("APOBEC3A", "APOBEC3B", "MKI67", "KRT13")

# run the function for those genes
grpd_gene_list_data <- prepare_data_for_plot(txiDF, gene_list)

# check you've only got the data you want
print(grpd_gene_list_data)

# plot
ggplot(grpd_gene_list_data, aes(x = genes, y = Mean, fill = GroupPrefix)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), 
                width = 0.2, position = position_dodge(0.9)) +
  labs(x = "Gene", y = "TPM", title = "Indicator gene expression") +
  theme_minimal() +
  scale_fill_manual(values = c("BKPyV" = "red", "uninfected" = "blue"))

# save
ggsave("indicator_genes.pdf")

# quit R
# don't save your workspace - it's usually unnecessary
q()

```

<p align="justify">
What does the (limited) data suggest to you?<br/><br/>
We used a function so we didn't have to do the same (long) block of R code for each gene we wanted to plot. It's possible you may want to mess with the y axis scale (log, perhaps?) as the genes have quite different expression levels.<br/>
The code will also generate bars with standard deviation errors (but we only had one sample in each group).<br/><br/>
Whilst this code is useful here to just check your dataset makes sense, it could also be useful code to use in your report if you want to highlight particular lists of genes you think are relevant. For now, you could try some other genes by making a new <code>gene_list</code> and re-running the code from there.<br/><br/>
Now that you've exited R, you can view your saved graph by running: <code>firefox indicator_genes.pdf</code>
</p><br/>

#### 3 Differential Expression Analysis (DEA) with Sleuth
<p align="justify">
Taking a quick look at TPMs is one thing, but we want to be unbiased and use the full statistical power of the dataset. For that we need DEA (done in R) using the full kallisto directories similar to those you made in part 2 (<i>not</i> the abundance.tsv files) and an input file showing our experimental design, which we will make now. This is another <code>for</code> loop, with each line doing text manipulation of the experiment metadata stored in our read names - this is why consistent and informative naming is so useful.<br/><br/>
A quick confession. You will now use some kallisto output which I created from the full dataset (not just the first 2 million reads as you did). This is because it would take more like 20-25 minutes per sample to run (rather than 2 minutes). If you do the RNAseq data for your report, you should run kallisto in full for all samples (I will give the full command at the end of the workshop material below). Part of this is we use a statistical method called bootstrapping to give us confidence in our expression values. Bootstrapping is where we re-run the data alignment multiple times to get consistent answers (like you would when building phylogenetic trees). In the full run you will do 20 bootstraps compared to the 0 you did in part 2. This explains why you <i>might</i> see (very) slight differences between your TPMs from part 3 and those of your classmates around you.<br/>
</p>

```sh
# enter the results directory for this section
cd ../3_DEA

# create the experimental design file for sleuth (the DEA tool)
# first write the headers
echo -e "sample\texpress\tbiorep\tpath" > infectionDEA.info

# now loop through the kallisto output directory filenames to get the information needed
for sampledir in ../rnaseq_data/01_workshop4_BKPyV-infection/03_full_kallisto_output_for_workshop_DEA/*/abundance.h5
  do
    sample=`echo $sampledir | rev | cut -d'/' -f2 | rev`
    donor=`echo $sample | cut -d'-' -f2`
    exp=`echo $sample | awk -v sample=$sample '{if (sample~/^BKPyV/) {print "0"} else {print "1"}}'`
    echo -e $sample'\t'$exp'\t'$donor'\t'$sampledir >> infectionDEA.info
done

```

<p align="justify">
The <code>for</code> loop looks very complicated. But look at each part to see what it is doing. We are looping ("iterating") through the abundance.h5 files in each sample directory (the variable <code>sampledir</code> here is out "iterable" - the variable we're using at each stage of the loop). The abundance.h5 files are made by kallisto. We then use the name of this directory to create the rest of the information we need. Some useful linux commands are used:<br/><br/>
<code>rev</code> - reverses the characters in your output. This is useful if you know you want the end of something, but the start may by of inconsistent length or format<br/>
<code>cut</code> - splits strings or data into lists based on a delimiter, this is the <code>-d</code> bit. We're splitting our string into a list on <code>/</code> or <code>-</code> in this loop, and then capturing the second (<code>-f2</code>) item in that list.<br/>
<code>awk</code> is a whole language in itself, like R, and is fantastic for manipulating files in columns. Here we're using it to check our sample names to put our samples into groups (either 0 or 1).<br/>
Then we use <code>echo</code> to print everything we want.<br/><br/>
Easy, really... The key thing is that this loop looks complex, but whenever you write your own loop you do it in stages: you don't just write that in one go. You check that you're looping through the right thing, that you're getting the information from each variable, then that your output looks right <i>before</i> you commit it to file. <br/><br/>
No we can use the output of the loop and jump back into R.
</p><br/>

```sh 
# enter R again from the terminal
R

```
```R
# access the Genomics shared libraries with this command
.libPaths("../rnaseq_data/R_4.3.3")

# load libraries we need
library(tidyverse)
library(readr)
library(sleuth)

# create an object for the experimental design and t2g
s2c <- read.table("infectionDEA.info", header=TRUE, stringsAsFactors=FALSE)
t2g <- dplyr::select(read.table("../rnaseq_data/gencode.v44.pc_transcripts.t2g", header=TRUE, stringsAsFactors=FALSE), target_id = ensembl_transcript_id, ext_gene = external_gene_name)

# build a sleuth object, aggregating to gene level
# the "summarizing bootstraps" step can take a few minutes
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
As we are running statistical tests for so many genes, we need to have a p value correction. Sleuth uses Benjamini-Hochberg, a ranking system to reduce false positives. We then use these 'q' values, again with a stringent threshold such as <0.05. Use <code>head</code>, <code>less</code>, <code>grep</code> or other tools you know to look at the significant genes - do any come up that you recognise? Do any of our "quick check" genes come up? If not, why do you think this is?<br/><br/>
The other thing to consider here is that q<0.05 genes are <i>statistically</i> different, but that doesn't necessarily equal biologically significant - the fold change in expression could be very small, but appear very significant just due to the power of the test (i.e. the ability of the test to see changes). This is particularly common when you have replicates from cell lines (as these are more like technical repeats than biological ones).<br/><br/>
Now we're going to combine our TPMs and our sleuth results and make a volcano plot.<br/>
</p>

```sh 
# create symbolic link for TPMs
ln -s ../rnaseq_data/01_workshop4_BKPyV-infection/03_full_kallisto_output_for_workshop_DEA/allTPMs.tsv allTPMs.tsv

# now get back into R
R

```

```R
# access the Genomics shared libraries with this command
.libPaths("../rnaseq_data/R_4.3.3")

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
library(ggrepel)
library(ggplot2)

# check how many significantly different genes we have, using log2FC threshold of >=1 (100% increase) or <=-1 (100% decrease) and qvalue of <0.05
sum(res$log2FC>=1 & res$qval<0.05)
sum(res$log2FC<=-1 & res$qval<0.05)

# create new column stating whether a gene is sig up or down, or not
res$DEA <- "NO"
res$DEA[res$log2FC > 1 & res$qval < 0.05] <- "UP"
res$DEA[res$log2FC < -1 & res$qval < 0.05] <- "DOWN"

# let's create our volcano plot - it should open as a separate window once rendered
# first line is the plot, second gives the colours and no legend, third defines the colours
# fourth line gives the threshold lines, fifth gives a clean white background, sixth sets the labels 
ggplot(res, aes(x=log2FC, y=-log10(qval))) + 
  geom_point(aes(colour = DEA), show.legend = FALSE) + 
  scale_colour_manual(values = c("blue", "gray", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") + geom_vline(xintercept = c(-1,1), linetype = "dotted") + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text_repel(data=subset(res, (abs(log2FC) > 1 & qval < 0.05)), aes(x=log2FC, y=-log10(qval), label=genes), max.overlaps = 1000, size=2.5)

# save it
ggsave("volcano.pdf")

# don't close R!

```

<p align="justify">
You now have a volcano plot - the hallmark of an RNAseq experiment! In a volcano you plot log2 of the fold change (x axis) against the -log10 of the significance value (y axis). These two log transformations are really important for actually understanding the data.<br/><br/>
We use log2FC so that up and down changes are symmetrical around zero. If your expression goes from 10 to 20, this is a fold change of 2, but if you go from 10 to 5 this is a fold change of 0.5. Whilst the fold change should be symmetrical, 0.5 is closer to 1 than 1 is to 2. This makes your plots wonky, and makes the up changes seem more important. If you do a log2 transform, a fold change of 1 (i.e. no change) gives a log2FC of 0. log2 of 0.5 is equal to -1 and log2 of 2 is equal to 1. Now we have a symmetrical plot.<br/><br/>
We want to have small p/q values to give us confidence differences between conditions are not just due to chance. But, we want these important points to be highlighted at the top of our graph. If we plotted fold change against p or q values, all the important stuff would be on or very near the x axis (squished at the bottom) and all the non-significant stuff would be at the top. So, we do -log10 of the stat test value, and this puts the rubbish at the bottom and the (hopefully) interesting genes at the top.<br/><br/>
One important note. We did include a +1 transformation when calculating fold change to get rid of seemingly huge changes between decimals (the fold change between 0.000005 and 0.0001 is huge!) but some of your most significant changes may come from very low expression. Is this biologically meaningful? Think about times when genes with low expression can still be incredibly important for cell function.<br/><br/>
Some comparisons will yield thousands of significant differences. This dataset is more manageable. One approach is to go through the genes which have changed and try and see patterns - do they make sense? This is a lot of work, even with just a few genes.<br/>
A more powerful and unbiased technique is to use gene set enrichment analysis (GSEA). Instead of just focusing on genes we consider significant (based on arbitrary thresholds), GSEA looks for patterns in the dataset as a whole. We'll do this now.
<br/>
</p>

```R 
# still in R

library(fgsea)

# use our pi values to rank the genes biologically from most up to most down
prerank <- res[c("genes", "pi")]
prerank <- setNames(prerank$pi, prerank$genes)
str(prerank)

# run GSEA using a list of gene sets curated by MSigDB (Broad Institute)
genesets = gmtPathways("../rnaseq_data/h.all.v2024.1.Hs.symbols.gmt")

# run GSEA - this will produce warnings, don't worry!
fgseaRes <- fgsea(pathways = genesets, stats = prerank, minSize=15, maxSize=500)

# check the top most significant hits
top10_fgseaRes <- head(fgseaRes[order(pval), ], 10)
top10_fgseaRes

# plot an enrichment plot for the top hit, and explore some of the others
# google the terms, or check the MSigDB website - do the results make sense?
plotEnrichment(genesets[["HALLMARK_E2F_TARGETS"]], prerank) + labs(title="E2F targets")
ggsave("E2Ftargets.pdf")

# create a bar chart of the best hits
ggplot(top10_fgseaRes, aes(x = NES, y=reorder(pathway, -pval), fill = factor(sign(NES)))) +
	geom_bar(stat = "identity", width = 0.8) +
	labs(title = "GSEA", x = "Normalised Enrichment Score (NES)", y = "Pathway") +
	theme_minimal(base_size = 16) +
	scale_fill_manual(values = c("#0754A2", "#B10029"), guide = "none") +
	scale_y_discrete(labels = function(x) gsub("^HALLMARK_", "", x)) +
	theme(axis.text = element_text(color = "black"),
		  axis.title = element_text(color = "black"),
		  panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank())

# save bar chart
ggsave("GSEA-hallmarks_best_hits_bar.pdf")

# save your results
write.table(res, file="GSEA-hallmarks_results_file.tsv", sep="\t", row.names=FALSE, col.names=TRUE)

# what's your interpretation?

# when you're happy, the workshop is over - you can quit R!
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
Hopefully you've found this session useful, interesting and inspirational (plus the other learning objectives...). If you would like to do RNAseq analysis for your Genomics3 assessment, you <b>should not</b> use the data from above. Instead I have provided with a similar dataset of human urothelial cell cultures grown in either normoxic (20% O<sub>2</sub>) or hypoxic (2% O<sub>2</sub>) conditions:<br/>
<code>/shared/biology/bioldata1/bl-00087h/data/rnaseq_data/02_hypoxia_assessment_dataset/</code>
<br/><br/>
This is a 4vs4 dataset where urothelial cells from 4 different people have been used for cell culture. Cells were re-differentiated to form a biomimetic tissue and then split, with one half continuing to be grown in standard normoxic conditions (20% O<sub>2</sub>) for 14 days, and the others grown in hypoxia (2% O<sub>2</sub>) for 14 days. RNA was then harvested and suvbmitted for sequencing. We are able to measure the ability of biomimetic urothelium to form a tight barrier (to urine, if it was in the body), and the barrier ability was highly compromised when the cells were taken into hypoxia.<br/>
In healthy people, the urothelium is an epithelial layer some 3-6 cells thick (depending on how full the bladder is) which sits on top of a basement membrane associated with a capillary bed (<i>i.e.</i> it has good oxygen supply). For this report consider what reductions in oxygen would mean for a tissue, its physiology and the resultant transcriptome. Consider how this could impact bladder health, or in which bladder diseases hypoxia may play a role. Are there major regulators of hypoxia which do not change at the transcript level?<br/><br/>
For your analysis, run through the same process as in this workshop: FastQC for looking at the reads, kallisto for the alignment, tximport in R for getting the gene-level TPM files, sleuth in R for running differential expression analysis, and then interpretation by plotting individual genes, making a volcano plot and using GSEA to report a little on the biology of the dataset. Look at the very bottom of this page for a couple of extra commands to help you.<br/><br/>
The way to access the higher marks is to go beyond what we show you in the workshop material. This includes biological interpretation, but also showing us some independent bioinformatics. You could consider some of the following:<br/>
<ol>
     <li>Plotting individual genes of interest from the TPMs </li>
	 <li>Specifically plotting/investigating pathways/sets of genes, such as HIF targets </li>
	 <li>Expanding your GSEA analysis using other gene lists from <a href="https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp">MSigDB</a> (you will need to create a free account) - the c2 set is in the rnaseq_data directory already</li>
	 <li>Generating informative plots which we haven't show you </li>
	 <li>Investigating the relevance of interesting genes in other urothelial settings, such as diseases where hypoxia may be important, like cancer - you could try <a href="https://www.cbioportal.org/">cbioportal</a> </li>
   </ol>
<br/>
There is more information on the assessment criteria on the VLE. Essentially we are looking for a decent introduction to the topic and relevance of the dataset, an explanation of your methods, and then your attempt to interpret the results and suggest how they address the purpose of the study.<br/>
Good luck, remember to support each other using the discussion boards, and do ask us for help and guidance if you need it.<br/><br/>
</p>

#### Some extra commands to help
<p align="justify">
Remember, you should run FastQC on each read file (all 16). You don't need to put any of the graphs into your report, but it is always good practice to check the quality of your data, and putting information on average read number is often good practice.<br/>
You then need to run kallisto on each sample (all 8) using both read 1 and read 2. You will also include bootstrapping here (which you didn't do in the workshop) to make your DEA more robust. The command for that is below. Then your Sleuth, plotting and GSEA commands should be very similar, just changing the read names and any relevant paths to data.<br/>
fastQC on these full files will take 10-20 minutes per file. kallisto will take 40-60 minutes per sample. Use loops, and set the commands running in the background, or using <code>screen</code> - google this it will be incredibly helpful. <br/><br/>
Remember, you can set jobs running, and then go off and do something else while they run. Setting something running at 5pm often means it is ready for you in the morning.<br/>
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
kallisto quant --index ../rnaseq_data/gencode.v44.pc_transcripts.processed-kallisto --output-dir=SAMPLENAME --bootstrap-samples=20 ./SAMPLENAME_read1.fq.gz ./SAMPLENAME_read2.fq.gz 
```

