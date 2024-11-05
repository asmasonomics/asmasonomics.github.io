---
layout: page
permalink: /courses/Genomics3_Workshop5_scRNAseq_Nov2024
---

![Genomics3 banner](/assets/images/genomics_banner.jpeg){:class="img-responsive"}

<span style="font-size:1.6em;">**Genomics 3 - Workshop 5: scRNAseq**</span><br/>

<p align="justify">
Welcome to Workshop 5! In the last workshop we worked on <b>bulk</b> RNAseq data - <i>i.e.</i> we derived a single, averaged transcriptome for entire populations of cells, and then compared different conditions to see what was different. In a cell culture system a homogenous response may be biologically accurate (although not always...), but in reality cells are independent entities, even highly similar cells. Differences may be minor, almost technical, based on cell cycle stage for example, or it could be that your population has completely different cell types in it. The latter context is what we will consider in this workshop as we work with single cell RNA sequencing (scRNAseq) data derived from human bladder cancers - "urothelial carcinoma".<br/><br/>
</p>
{% include image_cols.html 
	file = "/assets/images/labelled_urinary_tract_wb_smaller.jpg"
	content = "Cancers are complex ecosystems. For bladder cancer, the actual cancer consists of dysregulated urothelial cells - the epithelial cells which line the bladder wall (and the upper urinary tract up the ureters and into the renal pelvis of the kidney). But surgical samples do not just contain these urothelial carcinoma cells. We also capture white and red blood cells, fibroblasts, muscle <i>etc. etc.</i> Some of these cells may be genuinely part of the tumour microenvironment (TME) or it could be that the surgeon also sampled some normal(ish) tissue which was next to the main tumour (paracancerous). This is really important as patients with tumours which are heavily infiltrated (lots of immune cells within the TME) tend to do better than those where the tumour is evading immune surveillance. Infilitration is often used as a marker of treatment success (or resistance)."
%}
<br/>
Whilst a bulk transcriptome averages out the transcriptomes of all cells in a population (<i>e.g.</i> making a very murky cancer plus immune plus RBCs plus muscle picture, where the murkiness may be different between patients, but only for technical reasons, not biological), scRNAseq allows you to identify different cell types in your data. You can work with these independently (effectively filtering out the noise) or you can study interactions, changes in function <i>etc. etc.</i> - we're still only in the infancy of working out how scRNAseq can aid our research questions. <br/><br/>
One thing to remember is that scRNAseq still has some drawbacks:<br/>
<ol>
     <li>scRNAseq still loses all positional information from a dataset, for that you need spatial transcriptomics (and lots of £££)</li>
	 <li>You can look at heterogeneity, but you can only sample a (relatively) small number of cells (typically ~10k per sequencing run - bulk will have RNA from a population often in excess of 1M cells)</li>
	 <li>There is a limit of detection - genes with low expression will be poorly represented</li>
   </ol>
<br/>
These are important points to consider as you do your analysis (and when you're interpreting the data).<br/><br/>
</p>

### Introduction to the material
<p align="justify">
scRNAseq can be massive, so we're starting with data which has already mapped to the human transcriptome. In this workshop you will work exclusively in <code>R</code>, but using RStudio (nice and familiar), but still on the Linux system. The workshop is aimed towards the dual-boot Linux machines in G/N/169, but in theory you <i>could</i> do this analysis on a Windows machine. We <b>do not</b> recommend this however, as working on the managed machines means you have access to all the data, any workshop-required R libraries have been installed for you already, and you have increased computational power working on Linux rather than Windows.<br/><br/>
In the workshop you will start with genome-mapped data. You will perform QC and filtering to keep only high quality and informative cells. You will perform dimension reduction, clustering, community annotation, differential expression, gene set enrichment analysis and functional annotation - lots of graphs.<br/><br/>
If you choose to base your final report on <i>this</i> workshop, you will need to expand/adapt the analysis in the workshop to some related, but different data we have provided. The bioinformatic approach will be very similar, but you will need to address an appropriate question for your chosen dataset, and bring in the relevant biology. More details on these options are at the end of the workshop material.<br/>
</p>

### The workshop
#### Research Question
<p align="justify">
Cancer "T" stage is determined by how invasive the cancer has become. If a growth of cells is contained to its layer of origin (the urothelium in our case, including the bladder lumen as this is acellular space, not another tissue), it isn't <i>really</i> a cancer. These T<sub>a</sub> "tumours" are considered worth removing in bladder cancer because having a T<sub>a</sub> tumour is a major risk factor for developing subsequent invasive disease. T1 tumours invade through the basement membrane into the lamina propria (the stroma, where blood cells, fibroblasts, neurones <i>etc.</i> are). T2 and T3 tumours are different extents of muscle-invasive disease, and T4 tumours have escaped the bounds of the bladder, with metastasis highly likely if not evident at bladder cancer diagnosis.<br/><br/>
This is summarised nicely in this <a href="https://doi.org/10.1038/nrc3817">2015 bladder cancer review</a> with Figure 1 included below. 
</p><br/>
![Stages of bladder cancer from Knowles & Hurst, 2015](/assets/images/Knowles2015_PMID25533674_Fig1_BLCA_stages.jpg){:class="img-responsive"}
<br/>
<p align="justify">
Around 80% of people present at non-muscle-invasive stage (NMIBC) and 5-year survival is pretty good. However muscle-invasive (MIBC) disease has a 5-year survival lower than 50%, even with the standard-of-care radical treatment of removing the bladder. In this workshop you will work with scRNAseq data from a T1 tumour and a T3 tumour to identify whether there are differences between the malignant cells and tumour microenvironment between these tumour stages. These data were derived from Lai <i>et al.</i> (2021) published in the <a href="https://doi.org/10.1002/ijc.33794"><i>International Journal of Cancer</i></a>.
<br/><br/>

#### Workshop Aims
1. Quality check your scRNAseq data, removing uninformative or low quality cells
2. Apply appropriate techniques for identifying distinct clusters of similar cells
3. Explore the biology of these communities through annotation, differential expression and functional characterisation
4. Compare alike cells between two different samples
<br/><br/>

#### 0 Set up your RStudio project
<p align="justify">
<b>Remember, you should be working on a Linux machine in G/N/169 for this workshop.</b><br/>
Assuming that you are, open the Apps menu and open RStudio. This should feel familiar. You need to use all the good practices your learned in first and second year - do revisit Emma Rand's first year material and your BABS module VLEs if you need a refresher.<br/><br/>
####Give specific advice on RStudio on G/N/169 machines - paths etc. <br/><br/>
Finally, open a new script file called <code>Workshop5.R</code> or similar, and add the following lines to the top:<br/>
</p>
```R
# access the Genomics shared libraries with this command
.libPaths("/shared/biology/bioldata1/bl-00087h/data/singlecell_rnaseq_data/R_4.3.3")

# load relevant libraries
library(tidyverse)
library(seurat)
library(clustree)


```
<br/><br/>

#### 1 Relevant Subheading
<p align="justify">
Placeholder text.
<br/><br/>
</p>

#### 2 Relevant Subheading
<p align="justify">
Placeholder text.
<br/><br/>
</p>

### Concluding remarks
<p align="justify">
Placeholder text.
<br/><br/>
</p>

### What to do if you want to do RNAseq for your report
<p align="justify">
Placeholder text.
<br/><br/>
The way to access the higher marks is to go beyond what we show you in the workshop material. This includes biological interpretation, but also showing us some independent bioinformatics. You could consider some of the following:<br/>
<ol>
     <li> </li>
	 <li> </li>
   </ol>
<br/>
There is more information on the assessment criteria on the VLE. Essentially we are looking for a decent introduction to the topic and relevance of the dataset, an explanation of your methods, and then your attempt to interpret the results and suggest how they address the purpose of the study.<br/>
Good luck, remember to support each other using the discussion boards, and do ask us for help and guidance if you need it.<br/><br/>
</p>

#### Some extra comands to help
<p align="justify">
Placeholder text.
<br/><br/>
</p>