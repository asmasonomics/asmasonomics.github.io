---
layout: page
permalink: /courses/BABS4_Biochem_DataWorkshop4_March2025
---

![BABS4 banner](/assets/images/babs4_banner.jpg){:class="img-responsive"}

<span style="font-size:1.6em;">**BABS4 - Data Workshop 4**</span><br/>

<p align="justify">
Welcome to Data Workshop 4!<br/>
This is the second of two RNAseq workshops as part of the BABS4 (66I) "Gene expression and biochemical interactions strand". This material carries on directly from the material in workshop 3. If you didn't attend workshop 3, go back to that material and then follow straight on into this workshop - <a href="https://asmasonomics.github.io/courses/BABS4_Biochem_DataWorkshop3_March2025">please follow this link to the workshop 3 material</a>.<br/><br/>
Remember, during this module you are interested in <b>HiGam</b> and whether it plays a role in a <b>non-canonical</b> competence response. In these two data analysis workshops you will study the <b>canonical</b> competence response - does HiGam have a role?<br/><br/>
<br/>
</p>

### Introduction
<p align="justify">
In data workshop 3, you took a full RNAseq dataset from the <a href="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0217255">Black <i>et al.</i> 2020 paper</a> and eventually created a volcano plot showing the impact on the <i>Haemophilus influenzae</i> (Hi) transcriptome after cultures were transferred to starvation conditions for 30 minutes. A lot had changed, but what did it actually mean? This workshop is all about getting you to consider <b>the biology</b> so you can see how linking in this "big data" analysis can help you contextualise and understand your wet lab practicals and (hopefully) give you a better understanding of why natural competence is important.<br/><br/>
You can download the introductory slides <a href="/assets/coursefiles/2025-03_66I_replacement_plots/66I-DW4_introductory_slides.pdf" download>here as a PDF</a>.<br/><br/>
Please do ask us <b>biological</b> questions in this workshop. You have all the coding skills to explore and graph these data - so discuss with us possible questions you could ask to extend your understanding. Consequently, if you're having issues with installations on your personal laptop, we will largely direct you to the <a href="https://asmasonomics.github.io/courses/BABS4_Biochem_DataWorkshop3_March2025#access-the-data">advice in the workshop 3 material</a> or to use a managed machine (or the virtual desktop service).<br/><br/>
Remember, the data we are using comes from RNAseq read counts. There is a video explainer on what a sequencing read actually is. You can find this at the bottom of this page.<br/>
</p>
[Introduction to transcriptomics](#introduction-to-transcriptomics)<br/>
[How has the data in this workshop been processed so far?](#data-processing-before-this-workshop)<br/>
[What is a sequencing read?](#what-is-a-sequencing-read)

### The workshop
<p align="justify">
<b>REALLY IMPORTANT THING TO READ</b><br/>
Your RStudio Project submission should treat data workshops 3 and 4 as a single entity. This means ideally you should have a single script file and variable names which are consistent between the weeks.<br/><br/>
Remember, completing just the material from these two workshops is sufficient for a pass, but not to hit the excellence criteria. You should <b>expand your analysis beyond the workshop</b> and ensure that your code and project forms a coherent, single entity. I would recommend you create a fresh project for the assessment which runs successfully.<br/><br/>
</p>

#### Setup
<p align="justify">
For this workshop you should continue on from workshop 3. If you haven't completed workshop 3, <a href="https://asmasonomics.github.io/courses/BABS4_Biochem_DataWorkshop3_March2025">please follow this link to the workshop 3 material</a>.<br/><br/>
<b>New raw data for this session</b><br/>
Download these into your <code>raw_data</code> directory, either by downloading and copying it across or doing it programmatically (see below).<br/><br/>
<a href="/assets/coursefiles/2024-03_66I/Hi_GC_1kb.bed" download>Hi_GC_1kb.bed</a>. This file splits the Hi genome into 1000 bp (1kb) bins and has the GC% for each bin. Tab separated: chromosome, start, end, GC%.<br/><br/>
I will highlight any libraries you need as we go along.<br/><br/>
</p>

#### Load your data
```R 

# load your session from last time
load(BABS4_workshop3_complete.RData)

# if you didn't save your last workshop as an .RData, but you did save the files, you can load each of the raw files again, and your processed data too
# for example: tpms <- read.table("proc_data/full_dataset_TPMs.tsv", row.names = 1, header = TRUE)

# download and load new raw data for this session (for the circos plot)
download.file("https://asmasonomics.github.io/assets/coursefiles/2024-03_66I/Hi_GC_1kb.bed", destfile = paste(getwd(),"raw_data","Hi_GC_1kb.bed", sep="/"))
gc <- read.table("raw_data/Hi_GC_1kb.bed", header=FALSE, col.names = c("chr", "start", "end", "GC"))

```
<br/>

#### Check our understanding from data workshop 3
<p align="justify">
By the end of data workshop 3 you created a volcano plot highlighting the transcriptomic changes in Hi when grown in starvation conditions (MIV). Remember that the way we induce the Hi competence response in the lab is by starving the cells.<br/>
Starving the cells obviously has a knock on impact - the cells are not just inducing the competence response, but are also turning on pathways needed to synthesise/process the metabolites they need which were previously coming from the rich media (BHI). One of the challenges is we want to know what is the competence response, and what is a general stress response. <br/>
</p>

```R 

# let's regenerate your volcano plot from last time, just to check everything is working
library(ggplot2)
library(ggrepel)

# define the list of significant genes to add as plot labels
withsymbols <- dds_tpm[- grep("[.]", dds_tpm$symbol),]
siggenes <- head(withsymbols |> arrange(desc(abs(log2FC))), 60)$symbol

# generate the plot
ggplot(dds_tpm, aes(x=log2FC, y=-log10(padj))) + 
  geom_point(aes(colour = DEA), show.legend = FALSE) + 
  scale_colour_manual(values = c("blue", "gray", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") + geom_vline(xintercept = c(-1,1), linetype = "dotted") + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text_repel(size=3, data=subset(withsymbols, abs(log2FC) > 3), aes(x=log2FC, y=-log10(padj), label=symbol), max.overlaps = Inf)

```
![MIV0 MIV2 volcano](/assets/coursefiles/2025-03_66I_replacement_plots/03_dea_002.png){:class="img-responsive"}

<p align="justify">
As we explored last time, there are geners upregulated here related to competence (<i> e.g. comA</i>) and DNA binding and recognition (<i> e.g. dprA</i>), but there are also genes involved in tryptophan metabolism, purine metabolism and carbohydrate metabolism,  among others. Are these genes simply a starvation response?<br/><br/>
We could test this (slowly) by looking up genes individually. Another (faster) way to do this is gene set enrichment analysis (GSEA). This is very similar to gene ontology (GO) analysis, but standard GO typically only considers unordered gene lists, rather than lists ranked by fold change as you could provide after RNAseq.<br/>
GSEA is not part of today's workshop, but you could include it as an optional extra in your report - there are some pointers at the end of the workshop material. Today we will instead use some of the other tested conditions from the Black <i>et al</i> paper to further explore what is really happening in canonical competence. 
<br/><br/>
</p>

#### Nutrient deficiency response without the competence response?
<p align="justify">
An unfortunate consequence of our model system is that we induce many transcriptomic changes by altering the media - not just the competence response. We want to study <i>just</i> the biology of competence, so can we control for this added noise in the system?<br/>
One way is to include another differential expression analysis where only the nutrient stress response may be present, effectively allowing you to subtract the starvation response from the competence + starvation response we already have.<br/>
From our other available data (check back to the workshop 3 introduction and the table below), we have cells growing in rich media (BHI) but where the cell density is getting high - so where nutrients may be becoming sparse (without being <i>so</i> sparse it induces the same response as being in the starvation MIV media). <br/>
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
We will now run DEA on the BHI3 (most dense) cells against the MIV0 control (not dense; same control as our original) to try and tease out the specific competence response. <b>Think about this</b> - is this a good comparison to make, is it the most informative from the data we have?<br/>
Regardless, let's jump in.<br/>
</p>

```R 

# this should follow fairly closely to workshop 3
# think about what is happening at each step, and ask us about the biology if needed!

# extract appropriate TPM subset
condcomp_tpms <- tpms[, grep("kw20.MIV0|kw20.MIV2|kw20.BHI3", colnames(tpms))]

# check sample diversity with PCA
condcomp_pca <- princomp(log10(condcomp_tpms+1))
condcomp_pcavar <- round((data.frame(condcomp_pca$sdev^2/sum(condcomp_pca$sdev^2)))*100,2)
colnames(condcomp_pcavar) <- "var_perc"
condcomp_pcavar

# create variables needed for plotting
condcomp_pcacomps <- data.frame(condcomp_pca$loadings[, c("Comp.1", "Comp.2")])
pca_x <- paste(c("PC1 ("), condcomp_pcavar["Comp.1",], c("%)"))
pca_y <- paste(c("PC2 ("), condcomp_pcavar["Comp.2",], c("%)"))

# plot the samples using PC1 and PC2 (as these explain most of the variance)
ggplot(condcomp_pcacomps, aes(x=Comp.1, y=Comp.2)) + 
  geom_point() + 
  geom_text_repel(label=rownames(data.frame(condcomp_pca$loadings[, 1:2]))) +
  labs(x = pca_x, y = pca_y)

```
![three sets PCA](/assets/coursefiles/2024-03_66I/plots/04_pca_001.png){:class="img-responsive"}
<p align="justify">
<details>
   <summary>What does the PCA tell us about our samples? (check back to last week if you need a comparison)</summary>
   MIV0 and MIV2 appear very different (as we saw last week), but the MIV0 and BHI3 samples are interspersed. This could suggest they are very similar, and therefore unlikely to show many significant differences.<br/>
   You could drop the MIV2 samples and re-do the PCA - does this show a difference now?
   <br/><br/>
</details>
</p>

```R 

# PCA doesn't tell you everything - let's carry on and run DEA for these two conditions

# install DESeq2 if you need to (skipping updates again)
#BiocManager::install("DESeq2")
library(DESeq2)

# extract count subset
MIV0BHI3_counts <- counts[, grep("kw20.MIV0|kw20.BHI3", colnames(counts))]

# generate sample info for DEA
sample_info <- data.frame(colnames(MIV0BHI3_counts))
colnames(sample_info) <- c("sample")
rownames(sample_info) <- sample_info$sample
sample_info <- sample_info |> separate(sample, c("genotype", "condition", "replicate"))

# if you run into issues with separate the following line replaces the whole block above
# sample_info <- read.table(text=gsub("[.]", ",", colnames(MIV0BHI3_counts)), sep=",", col.names=c("genotype", "condition", "replicate"), row.names = colnames(MIV0BHI3_counts))

# generate DESeq2 object
MIV0BHI3_dds <- DESeqDataSetFromMatrix(countData = MIV0BHI3_counts, colData = sample_info, design = ~ condition)
MIV0BHI3_dds$condition <- relevel(MIV0BHI3_dds$condition, ref = "MIV0")
MIV0BHI3_dds <- DESeq(MIV0BHI3_dds)
MIV0BHI3_dds_results <- results(MIV0BHI3_dds)

# merge with features
MIV0BHI3_dds_results <- merge(as.data.frame(MIV0BHI3_dds_results), featname, by=0)
rownames(MIV0BHI3_dds_results) <- MIV0BHI3_dds_results$Row.names
MIV0BHI3_dds_results <- MIV0BHI3_dds_results[,-1]

# subset TPMs, calculate log2FC, merge 
MIV0BHI3_tpms <- condcomp_tpms[, -grep("MIV2", colnames(condcomp_tpms))]

# due to a dependency issue, explictly call dplyr in the select call
MIV0BHI3_tpms <- mutate(MIV0BHI3_tpms, MIV0_avg = rowMeans(dplyr::select(MIV0BHI3_tpms, contains("MIV0"))), 
                        BHI3_avg = rowMeans(dplyr::select(MIV0BHI3_tpms, contains("BHI3"))))
MIV0BHI3_tpms <- mutate(MIV0BHI3_tpms, log2FC = log2((BHI3_avg+1)/(MIV0_avg+1)))

MIV0BHI3_dds_results <- merge(MIV0BHI3_dds_results, MIV0BHI3_tpms, by=0)
rownames(MIV0BHI3_dds_results) <- MIV0BHI3_dds_results$Row.names
MIV0BHI3_dds_results <- MIV0BHI3_dds_results[,-1]

# extract genes with symbols (and not rRNA genes), and then the most significant genes
MIV0BHI3_withsymbols <- MIV0BHI3_dds_results[- grep("[.]|HI_", MIV0BHI3_dds_results$symbol),]

# set DEA column
MIV0BHI3_dds_results$DEA <- "NO"
MIV0BHI3_dds_results$DEA[MIV0BHI3_dds_results$log2FC > 1 & MIV0BHI3_dds_results$padj < 0.05] <- "UP"
MIV0BHI3_dds_results$DEA[MIV0BHI3_dds_results$log2FC < -1 & MIV0BHI3_dds_results$padj < 0.05] <- "DOWN"

# plot with ggplot()
ggplot(MIV0BHI3_dds_results, aes(x=log2FC, y=-log10(padj))) +
    geom_point(aes(colour = DEA), show.legend = FALSE) +
    scale_colour_manual(values = c("blue", "gray", "red")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") + geom_vline(xintercept = c(-1,1), linetype = "dotted") +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    geom_text_repel(size = 3, data=subset(MIV0BHI3_withsymbols, abs(log2FC) > 1), aes(x=log2FC, y=-log10(padj), label=symbol), max.overlaps = Inf)

```
![MIV0 vs BHI3 volcano](/assets/coursefiles/2025-03_66I_replacement_plots/04_dea_002.png){:class="img-responsive"}
<p align="justify">
Another volcano! But this time far fewer significantly different genes.<br/><br/>
Those that are different, again highlight increases in carbohydrate and amino acid metabolism, but we don't see the competence or DNA recognition genes - perhaps this control has worked!<br/>
Look at the downregulated genes - mu!<br/>
<i>gam</i> doesn't come up at first glance, but muA, MuI and muL all do, and they're in the same operon. Interesting!<br/><br/>
One thing to say here is that we are resolutely ignoring some genes with the biggest changes, both in last week's work and here. We are only looking at stuff with symbols, but we could be ignoring some interesting biology. These could be quite well studied genes, but they just haven't been named properly in Hi yet...  Let's flip our focus for a minute.<br/>
</p>

```R 

# focus only on genes without a symbol, using the HI identifiers in the row names
MIV0BHI3_withoutsymbols <- MIV0BHI3_dds_results[grep("[.]", MIV0BHI3_dds_results$symbol),]
MIV0BHI3_withoutsymbols_sig <- rownames(head(MIV0BHI3_withoutsymbols |> arrange(desc(abs(log2FC))), 10))

# plot the volcano if you like, make sure to change the lab and selectLab flags (think rownames again...)
# or just view the MIV0BHI3_withoutsymbols_sig dataframe
# HI_1456 and HI_1457 have massive changes!

```

<p align="justify">
So what are these? The HI_1457 gives a bit of a hint with "opacity protein". You could google this, but some features may be less informative - HI_1456 just says "predicted coding region", for example.<br/>
What I have also made available to you is the nucleotide sequence of each feature in your counts dataframe - <a href="/assets/coursefiles/2024-03_66I/Hi_feature_sequences.fa" download> Hi_feature_sequences.fa</a> - so you could open this file in notepad, or wherever and extract the nucleotide sequence you're after. You could then use this as a query for <a href="https://tinyurl.com/ncbi-blastx">NCBI's blastx</a>, which uses a translated version of your nucleotide query against a protein database (as in, it converts the nucleotides in your query to protein to get more diverse search results). What does this tell you - could this give you more information on some of the big changing genes in your DEA results?
<br/><br/>
</p>

#### Can we extract a competence-specific response?
<p align="justify">
Quick reminder of what we're attempting to do. In workshop 3, we looked at the induced competence response in Hi using starvation medium. But. This also gave a starvation response. Today, we've compared the same control samples to cells growing in greater density, so the cells may be starting to face nutrient limits, without it causing so much stress you get the full stress (starvation + competence) response we had before. <br/>
Now we have both these comparisons, we can attempt to subtract the nutrient limit response to leave just genes involved in competence.<br/>
</p>

```R

# rename columns we want to compare from each DESeq output
MIV0BHI3_DEA <- MIV0BHI3_dds_results[, c("log2FC", "padj", "symbol")]
colnames(MIV0BHI3_DEA) <- c("MIV0BHI3_log2FC", "MIV0BHI3_padj", "symbol")
MIV0MIV2_DEA <- dds_tpm[, c("log2FC", "padj")]
colnames(MIV0MIV2_DEA) <- c("MIV0MIV2_log2FC", "MIV0MIV2_padj")

# now names are unique, merge the columns of interest, extract features with symbols
dea_comp <- merge(MIV0BHI3_DEA, MIV0MIV2_DEA, by=0)
rownames(dea_comp) <- dea_comp$Row.names
dea_comp <- dea_comp[,-1]
dea_comp_symbols <- dea_comp[- grep("[.]", dea_comp$symbol),]

# correlate the log2FC values between the two comparisons
cor(dea_comp$MIV0BHI3_log2FC, dea_comp$MIV0MIV2_log2FC, method = c("pearson"))

# plot two DEAs against each other using log2FC values
ggplot(dea_comp, aes(x=MIV0BHI3_log2FC, y=MIV0MIV2_log2FC)) + 
  geom_point() + 
  geom_text_repel(size=3, data=subset(dea_comp_symbols, abs(MIV0MIV2_log2FC) > 3 & abs(MIV0BHI3_log2FC) < 1), aes(x=MIV0BHI3_log2FC, y=MIV0MIV2_log2FC, label=symbol), max.overlaps = Inf, colour="red") +
  geom_abline() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

```
![multi DEA comparison](/assets/coursefiles/2024-03_66I/plots/04_dea_003.png){:class="img-responsive"}
<p align="justify">
We have now compared the log2FC changes in two separate differential expression analyses.<br/><br/>
If the transcriptomic changes were identical in both conditions, all genes would fall on the indicated y=x line. Changes which are specific to each DEA fall directly on their respective axes - y axis for MIV0 vs MIV2 and x axis for MIV0 vs BHI3. Genes near the origin do not change (much) in either comparison. There are more changes (and of greater magnitude) for the MIV0 vs MIV2, so the plot is not square - always check the axis scale.<br/><br/>
I've chosen to highlight those genes which would not be considered as significantly different in the BHI3 comparison, but are in the MIV2 (<i>i.e.</i> genes close to the y axis). Here we see the competence genes like <i>dprA</i> and <i>comA</i>, and the competence transcription factor <i>tfoX</i>. Mess about with these thresholds and see where the carbohydrate biosynthesis genes have gone (hint - much closer to x=y).<br/><br/>
So now we have a much reduced and more specific-to-competence list of genes. We still have the <i>pur</i> and <i>trp</i> genes - does that mean they have a more specific role in competence, rather than these nutrients being absent from the media? I recommend contextualising your results with the literature (cough, mark scheme, cough) to see whether these genes make a lot of sense. <br/>
Also, would a similar comparison with any of the other datasets be equally/more informative on the competence response?<br/><br/>
Note. You don't have to have a shared control condition (like MIV0 in our case) in order to compare. It <b>is</b> useful, as having a shared condition means you have a definite (rather than potentially assumed) shared baseline. Something to mention/highlight/control for.<br/><br/>
</p>

#### Full-genome summary visualisation with Circos plots
<p align="justify">
Now we have lots of information on which genes are changing we can try and summarise this data at the full genome level using a circos plot. This plots the Hi genome as a circle and then we can layer tracks on and around this circle to show the distribution of genes or particular features (<i>e.g.</i> regions or GC content). <br/><br/>
Circos plots can get massive, complex and strange very quickly (<a href="https://circos.ca/">check out some examples here</a>), often becoming too complex to actually be useful. But, they are used a lot in biological data visualisation, so we'd like you to have a go. Hopefully it will also help you to link with the PHASTER results from data workshop 2.<br/><br/>
A big problem I've faced when designing this material is that many circos plotting libraries have thousands of dependencies. We are therefore keeping it simple, so the plots you make are definitely customisable - but we're not expecting some of the incredibly complex images we can see on the front covers of scientific journals!
<br/>
</p>

```R 

# install and load the BioCircos library
# the vignette and further information are here - https://cran.r-project.org/web/packages/BioCircos/vignettes/BioCircos.html
install.packages('BioCircos')
library(BioCircos)

# we will now create a list of tracks, each time adding more information
# then we visualise it all in one go

# set the plot title in the middle of the circle (you may want to mess with x and y until you are happy)
tracklist <- BioCircosTextTrack("titletrack", "Hi kw20", opacity = 0.5, x = -0.2, y = 0)

# define an arc region to indicate where the mu prophage region lies
tracklist <- tracklist + BioCircosArcTrack("prophage region", "L42023.1", 1558774, 1597183, 
                                           opacities = c(1), minRadius = 1.25, maxRadius = 1.4,
                                           labels = c("mu prophage"))

# use the Hi GC content file (new raw data for this workshop) to create a GC content plot
tracklist <- tracklist + BioCircosLineTrack("GC", "L42023.1", gc$start, gc$GC, 
                                            minRadius = 0.4, maxRadius = 0.68, color = "black",
                                            labels = c("GC content"))

# we will now merge our MIV0 vs MIV2 (data workshop 3) DEA results with the feature locations, retaining useful columns
# this means we will be able to view our log2FC values as a genome-distributed heatmaps
dea <- merge(dds_tpm, featlocs, by=0)
rownames(dea) <- dea$Row.names
dea <- dea[,c("symbol", "log2FC", "chr", "start", "end", "strand")]

# use this to create a new heatmap track
# the colors form a red to blue gradient (from high pos to low neg log2FC)
# pick colours you like!
tracklist <- tracklist + BioCircosHeatmapTrack("DEA", "L42023.1", 
                                               dea$start, dea$end, dea$log2FC,
                                               minRadius = 0.8, maxRadius = 0.95,
                                               color = c("#FF0000", "#0000FF"),
                                               labels = c("MIV0vsMIV2 DEA"))

# I'm only showing one heatmap track here - this could be another good way to highlight regions of the genome 
# which are up or down in different condition comparisons...

# finally, render the circos plot
BioCircos(tracklist, genome = list("L42023.1" = 1830138), genomeLabelTextSize = 0,
          genomeTicksScale = 1e5, genomeTicksTextSize = 12)

# remember to save
ggsave("plots/circos.pdf")

```
![circos](/assets/coursefiles/2024-03_66I/plots/04_circos_001.png){:class="img-responsive"}
<p align="justify">
As so many genes were upregulated in the MIV2 condition relative to MIV0, the heatmap looks pretty warm. If you add another heatmap track with today's BHI3 vs MIV0 comparison, you should start to see those regions which are consistently up between conditions (<i>i.e.</i> carbohydrate metabolism) and those which are up only in the competence-inducing conditions (<i>i.e. comA, dprA</i>). Nice to have the GC plot too to get some genome-level descriptions. GC content always fluctuates, but the mean values are often quite species-specific. Interesting that the prophage region has a higher than Hi-average content - pretty good supporting evidence for this actually being a phage integration, rather than convergent evolution. You can include gene density or TF binding site density as line tracks in this way too.<br/><br/>
The <code>BioCircos</code> functionality is pretty limited, so you might want to label a few genes manually on top of your generated figure. Similarly, it would be good to add a little scale bar to the GC content track - you can extract the min and max values from the <code>gc</code> dataframe.<br/>
Sometimes coding a figure to exactly how we want it is very very hard, and not worth the effort when you can add labels manually. Always consider how much time you sink into a task (often compared with how many times you are likely to do that task).<br/><br/>
We can label up specific regions on a circos plot very easily, such as the prophage region. Think about how you might link figures from different parts of the practicals and data workshops. Remember, you don't <i>have</i> to report your results in the order you did them - think about the narrative. In this case the given order is pretty sensible - but always think about your structure for maximum clarity.<br/><br/><br/>
Now we've started to consider the spatial organisation of genes in the Hi genome, and we're looking at the prophage region again, we can ask whether the RNAseq data directly informs your outstanding diagnostic PCR questions.<br/><br/>
</p>

#### Can the sequencing data show us if muA, muB and gam are transcribed together?
<p align="justify">
An outstanding question from your diagnostic PCR practical is whether muA, muB and gam are explicitly transcribed together. Our short read sequencing data has the potential to answer this question by looking to see if these genes (and other in the mu prophage region) are more closely correlated with each other than with distant genes. Let's look.<br/>
</p>

```R 

library(sjmisc)

# rotate tpm dataframe to get features as columns
t_tpms <- tpms |> rotate_df()

# create a correlation matrix against the gam feature (HI_1483)
gam_tpm_cor <- cor(as.matrix(t_tpms[,c("gene-HI_1483")]), as.matrix(t_tpms))
gam_cors <- as.data.frame(gam_tpm_cor[,order(-gam_tpm_cor[1,])])
colnames(gam_cors) <- c("gam_cor")

# add in symbols where possible
gam_cors <- merge(gam_cors, featname, by=0)
rownames(gam_cors) <- gam_cors$Row.names
gam_cors <- gam_cors[,-1]
gam_cors <- gam_cors[order(-gam_cors$gam_cor),]
head(gam_cors, 25)

# look at all those mu genes, or those without symbols but numerically close HI IDs

# add in location information
gam_cors <- merge(gam_cors, featlocs, by=0)
rownames(gam_cors) <- gam_cors$Row.names
gam_cors <- gam_cors[,-1]

# create column with distance between gam st/end, accounting for circular genome
# gam is found at coordinates 1565297-1565806
# genome is 1830138 bp - features at coordinate 10 are closer if you measure across the 'top' of the circle
gam_cors <- mutate(gam_cors, gam_circ_dist = ((start+1830138)-1565806), gam_noncirc_dist = abs(1565297-start))
gam_cors <- transform(gam_cors, min_gam_dist = pmin(gam_circ_dist, gam_noncirc_dist))

# sort and subset columns before plotting
gam_cors_top <- (gam_cors[order(-gam_cors$gam_cor),])[,c("gam_cor", "symbol", "min_gam_dist")]
ggplot(gam_cors_top, aes(x=log10(min_gam_dist+1), y=gam_cor)) + 
  geom_point() + 
  geom_text_repel(size=3, data = gam_cors_top |> mutate(label = ifelse(gam_cor > 0.87, rownames(gam_cors_top), "")), aes(label = label), max.overlaps = Inf) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

```
![mu correlation by distance](/assets/coursefiles/2024-03_66I/plots/04_mu_001.png){:class="img-responsive"}
<p align="justify">
This plot shows correlation values between the expression of different genes and <i>gam</i> (y axis) and how far away from <i>gam</i> each of those genes is (x axis). Most genes are quite a long way from <i>gam</i> (>100,000 bases away; >10^5 bases) and have varying to little correlation. However, genes closest to the prophage region are well correlated with <i>gam</i>.<br/>
This suggests there is a strong positional effect to expression across the mu prophage region.<br/>
<details>
   <summary>But does this actually prove they are transcribed together?</summary>
   No! For bacteria, such a good correlation is <b>highly</b> suggestive, but again, non-conclusive.
   <br/><br/>
</details>
With short read sequencing, the reads themselves are too short to cover the span of distance. The only way would be to see read pairs spanning the region, but this (with a distance over a few hundred bases) is unlikely to occur due to size selection during library prep. <br/><br/>
<details>
   <summary>Can you think of a sequencing-based technology that could resolve this?</summary>
   Long read sequencing. In theory these methods capture the entire length of a transcript, even if that is an entire operon. For such a specific question (which could be solved with the right PCR primers...), long read sequencing is probably expensive overkill, unless you wanted to confirm other co-transcribed regions at the same time. Always think about narrow vs broad use in the context of financial decisions in experimental planning.<br/><br/>
</details>
<br/>
</p>

#### Finishing up for today
<p align="justify">
The aim of today was to reinforce and expand on data workshop 3, and start to properly integrate the biological interpretation of the data. Well done! Your coding skills have now taken you from a single dataframe of counts to multi-condition comparisons, gene set enrichment analysis and integration with the literature. This is a real expansion of your coding skills and an opportunity to have worked with big data.<br/><br/>
Before you finish today, make sure you understand <b>why</b> you have been doing this analysis - chat with the demonstrators - does it all link together in your head, particularly for making a single narrative with your wet lab work? Also: save, save, save and make lots of comments on your code!<br/><br/>
</p>

#### Expansions to the project
<p align="justify">
Some ideas for expanding your analysis to improve your understanding of Hi competence (and likely your grade simultaneously):
<ol>
  <li>Use other conditions in the provided dataset.</li>
  <li>Explore significant features which lack a gene symbol.</li>
  <li>Add in more plots which help you to explore the data - heatmaps, networks, bar graphs <i>etc.</i></li>
  <li>Further biological exploration - GSEA (see below), pathways, KEGG <i>etc.</i></li>
  <li>Anything else cool you can think of!</li>
</ol>
<br/><br/><br/>
GSEA (gene set enrichment analysis) is a really good way to explore what is happening in a transcriptome dataset. It stops you from looking up each gene one by one, but instead grouping genes together into pathways, functions or locations within a cell. The problem is that this only works really well with a model organism, and Hi is not a model organism. We can make use of databases from model organisms, such as <i>E. coli</i>, but there are caveats - <i>E. coli</i> is not naturally competent for example, so it simply doesn't have that pathway - so you can't find it in Hi using those databases.<br/>
A big caveat, but you can still find some interesting stuff. Major biochemical pathways will be covered well, and <i>E. coli</i> is studied a lot, so we should get the main players, but we have to be careful on our interpretation. Also, not all Hi genes have symbols. If we had a specific Hi gene set database, that's fine. But we need to convert to E. coli IDs, so we lose lots of information.<br/><br/>
The first time we ran this workshop we did some GSEA, but we found most of the class focused too much on the code and not the actual biology. The GSEA now forms an optional extra, but you can still make use of the code. The are some packages to install using <code>BiocManager</code>, then we create our curated and ranked gene list, then we run the GSEA algoirthm.
</p>

```R 

# install BiocManager if you need to, otherwise just load it
# install.packages("BiocManager")
library(BiocManager)

# install and load the following packages
# remember, skip any updates, usually by typing "n" when prompted
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)

# install and load the E. coli gene set information
organism  <- "org.EcK12.eg.db"
BiocManager::install(organism)
library(organism, character.only = TRUE)

## for clusterProfiler, we need an array with the log2FC values in descending order
# we want symbol names, and we need to get rid of genes where there is no symbol
# and get rid of tRNA and rRNA genes
# and get rid of any genes where the symbol appears more than once

# remove genes without a symbol
symbols_only <- dds_tpm[- grep("[.]", dds_tpm$symbol),]

# order the genes by symbol and by log2FC before removing duplicate lines
symbols_only <- symbols_only[order(symbols_only[,"symbol"],-symbols_only[,"log2FC"]),]
symbols_only <- symbols_only[!duplicated(symbols_only$symbol),]

# remove rRNA and tRNA
symbols_only <- symbols_only[- grep("ribosomal_RNA", symbols_only$description),]
symbols_only <- symbols_only[- grep("tRNA", symbols_only$symbol),]

# use this filtered set to create the gene_list for GSEA
gene_list <- symbols_only$log2FC
names(gene_list) <- symbols_only$symbol
gene_list <- sort(gene_list, decreasing = TRUE)

# run GSEA
# this will give warnings, but don't worry (unless there are errors, rather than warnings)
# run GSEA
gse <- gseGO(geneList = gene_list, ont = "ALL", keyType = "SYMBOL",
             minGSSize = 3, maxGSSize = 800, verbose = FALSE,
             pvalueCutoff = 0.05, OrgDb = organism, pAdjustMethod = "none")

# store useful results
gsearesults <- gse[,c("Description","setSize","qvalue","NES","core_enrichment")]

# create a list of the most significant (by q value)
gsearesults_top <- head(gsearesults[order(gsearesults$qvalue),], 40)
gsearesults_top <- gsearesults_top[order(gsearesults_top$NES, decreasing = TRUE),]

# take some time to view and explore these results - what do they tell you?

```

<p align="justify">
A few points on interpretation here. 
<ol>
  <li>GSEA works by ordering the genes (from biggest to smallest log2FC in our case) and then seeing if particular lists of genes ("gene sets") occur closer to either end of the ranked list than you would expect by chance.</li>
  <li>A gene set found mostly at the high log2FC end is an "enriched" gene set in our data, and the opposite for suppressed pathways.</li>
  <li>You could spend ages refining the gene sets - better to use this to get the global picture of the transcriptomic changes.</li>
  <li>The enrichment score (ES) is a score generated by how many genes appear close together in the rankings. This means you get higher scores if many genes are close together, but also if there are lots of genes in a set (say 200 vs a set of 8). This is why you have a normalised enrichment score (NES), which accounts for the differences in set size.</li>
  <li>The description can be a bit vague. Sorry! You can get more information from the "leading edge" (aka "core_enrichment") genes - as these are the ones which contribute to the ES. You'll see that many overlap between different gene sets.</li>
  <li>When you report GSEA results, give the gene set description, NES and q value (the adjusted p value).</li>
</ol>
<br/>
Let's summarise these GSEA results with some plots.<br/>
</p>

```R 

# create summary dotplot
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

# remember to save plots of interest
ggsave("plots/MIV0MIV2_GSEA_dotplot.pdf")

# create GSEA enrichment plots for individual processes
# these are based on gsearesults dataframe order, not gsearesults_top
# plot the most enriched gene set
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)

```
![MIV0 MIV2 gsea dotplot](/assets/coursefiles/2024-03_66I/plots/04_gsea_001.png){:class="img-responsive"}

<p align="justify">
Here we have the top 10 terms for enriched gene sets (activated) and suppressed ones. Size of the circle represents how many genes are in each set, with the gene ratio saying how many contributed to the ES in our data. The colour indicates the significance. There's lots of customisation available here.<br/><br/>
</p>
![carbhydrate metabolic process](/assets/coursefiles/2024-03_66I/plots/04_gsea_002.png){:class="img-responsive"}
<p align="justify">
This two panel plot is how the individual gene set enrichment scores are calculated. On the x axis at the bottom, the vertical lines show the position of each gene in the ranked log2FC order. The green line is the ES building up and up based on the spacing of gene set genes in the rank list, eventually peaking (this is the ES). The top panel shows the rank position against the metric for ranking (the fold change).<br/><br/>
The GSEA really highlights the starvation response - increases in sugar, amino acid and small molecule metabolism, and decreases in protein synthesis. Nothing about competence? But is that a surprise?
</p>
```R

# check for gam in any leading edge gene set
gsearesults[grep("gam", gsearesults$core_enrichment),]

# check for any pathway labelled competence
gsearesults[grep("Competence|competence", gsearesults$Description),]

# It wasn't there to be found...
# But, as predicted, GSEA (using E. coli genes...) does show the big metabolic shifts

```
<p align="justify">
Remember, GSEA is an optional extra!
<br/>
</p>



#### Tips for the RStudio Project
<p align="justify">
<ol>
  <li>One single cohesive project - likely with a single .R script file and sensible directories</li>
  <li>Clear commenting which is relevant to understanding the code</li>
  <li>Include plots and printing you did to check the data and make your next decisions - don't include the 6 times you misspelt something and got errors</li>
  <li>Include your <code>library()</code> lines, but comment out any installation lines</li>
  <li>Don't load in the 2 processed data files I gave you at the start of this workshop - work with the raw data</li>
  <li>Sensible naming, including for plots</li>
  <li>Clear and consistent style (spaces, indentations <i>etc.</i>) - keepi it legible.</li>
  <li>I won't penalise you for not using complex loops, only if code was repeated multiple times without need (efficiency and conciseness)</li>
</ol>
<br/>
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

#### What is a sequencing read
<p align="justify">
This video was included in the week 6 teaching materials on the VLE to support the workshop material. We had some confusion as to what a sequencing read actually is - so hopefully this explainer helps.<br/>
</p>
<iframe src="https://york.cloud.panopto.eu/Panopto/Pages/Embed.aspx?id=b2f4da0d-6a2c-4e94-a633-b13300fce3cc&autoplay=false&offerviewer=true&showtitle=true&showbrand=true&captions=true&interactivity=all" height="405" width="100%" style="border: 1px solid #464646;" allowfullscreen allow="autoplay" aria-label="Panopto Embedded Video Player"></iframe>
<br/>
[Return to tutorial introduction.](#introduction)<br/>
[Begin the workshop.](#the-workshop)<br/><br/>
