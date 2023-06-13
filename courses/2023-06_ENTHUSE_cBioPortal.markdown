---
layout: page
permalink: /courses/ENTHUSE_cBioPortal_Jun2023
---
<span style="font-size:1.6em;">**ENTHUSE bioinformatics session: using cBioPortal to explore *TP53* mutations**</span><br/>

<p align="justify">This teaching material was designed as the second part of Dr Simon Baker's ENTHUSE project <i>TP53</i> practical, providing hands-on bioinformatic analysis of real cancer sequencing data. This material can be used directly as a web broswer-based tutorial, or adapted with acknowledgement into a specific lesson plan.<br/></p>
[The video embedded below](#session-walkthrough) is a walkthrough of the exercises by Dr Andrew Mason.

### Introduction
<p align="justify">In your laboratory session you identified a mutation in <i>TP53</i> (R273H) using a restriction enzyme site created as a result of the mutation. Excellent work! By finding this specific mutation in the P53 protein ("the guardian of the genome"), you will be able to recommend <a href="https://pubmed.ncbi.nlm.nih.gov/33449813/">Eprenetapopt</a>, a drug which targets mutant versions of P53, reducing the tumour's ability to grow. Not bad for a simple lab-based assay!<br/><br/>
Using mutations to subtype different types of cancer improves our ability to treat cancer properly in the clinic. This is a move towards personalised medicine. But, how did researchers know this mutation existed, and how often do we find it in people with cancer? Commonly occuring mutations, known as <b>hotspots</b>, help us prioritise efforts for drug development, where one drug can help the most people.<br/><br/>
In this session you will explore <i>TP53</i> mutations in over 10,000 cancer samples using <a href="https://www.cbioportal.org/">cBioPortal</a>, a publicly available online resource for cancer genomics. Data mining is now a huge part of biomedical research, and the discipline of <b>bioinformatics</b> combines computational analysis with the interpretation of biological data. You will use the results of tumour DNA sequencing data from 32 different cancer types within <b>The Cancer Genome Atlas</b> (TCGA) study, an enormous international study to better understand the diversity of genetic changes in cancer.<br/></p>

### Learning objectives
1. placeholder text
2. placeholder text
<br/>

### The Session
<p align="justify">If you follow the guided steps below, you will explore <a href="https://www.cbioportal.org/">cBioPortal</a>, investigate the prevalence of <i>TP53</i> mutations across different cancers, and begin to see the complexity of cancer biology by seeing how different cancers exploit different parts of normal cell biology. Look out for 'explainer' links - these give more information and background to the concepts and data used in the session.<br/><br/>
I would recommend using a computer or a tablet for this session, rather than a phone. cBioPortal hasn't really been designed with phones in mind, and you'll find everything a bit small and tricky to navigate.<br/>
Take your time and explore the website and data - you could spend days on this website and still find brand new data. If you get lost, go back to the homepage and you're only a few clicks from where you need to be!<br/><br/></p>

#### 1 Open the cBioPortal website
<p align="justify">
Open up any web browser and head to <a href="https://www.cbioportal.org/">www.cbioportal.org</a><br/><br/>
The home screen (pictured below) gives immediate access to results from sequencing data generated from thousands of cancer samples from all over the world. Scroll down the list. Each study represents potentially years of work from doctors, surgeons and nurses in hospitals to recruit patients, work with families and collect samples, and then from biomedical researchers in universities to process the samples, then generate, analyse, interpret and publish the results.<br/></p>
![cBioPortal homepage www.cbioportal.org](/assets/images/ENTHUSE-01_cBioPortal_home.jpg){:class="img-responsive"}
<br/>
**EXPLAINER** [What is cBioPortal and where does the data come from?](#cbioportal-explained)<br/>
**EXPLAINER** [What is sequencing?](#sequencing-explained)<br/>
<br/><br/>

#### 2 Select the pancancer TCGA dataset
<p align="justify">
Return to the top of the cBioPortal homepage. We're going to work with The Cancer Genome Atlas (TCGA) pancancer study. "Pan" in this context just means "across lots of different cancers".<br/>
Click the Quick select link for <b>TCGA PanCancer Atlas Studies</b>. You have just loaded data for 32 different cancer studies and 10967 samples coming from 10528 different people! Next, click on the blue <b>Query By Gene</b> button, scroll to the bottom and type <b>TP53</b> into the Enter Genes box, then <b>Submit Query</b>.<br/>
</p>
**EXPLAINER** [What is The Cancer Genome Atlas?](#the-cancer-genome-atlas-explained)<br/>
**EXPLAINER** [Why are gene and protein names different?](#gene-and-protein-naming-explained)<br/>
<br/><br/>

#### 3 Explore the *TP53* OncoPrint plot
<p align="justify">
Lots of colours! Don't panic. An OncoPrint is a big summary of mutation data - is there a mutation in a particular sample, and what type of mutation is it? In an OncoPrint each column represents a single sample. You have loads of samples, so everything is very squished.<br/>
You can largely ignore the three "Profiled for..." rows - these just indicate which sequencing technologies were used on those samples. For example, you'll see a large number of glioblastoma samples were only profiled for copy number changes, not mutations or structural changes.<br/></p>
![TP53 OncoPrint from TCGA PanCancer Atlas](/assets/images/ENTHUSE-02_TP53_OncoPrint.jpg){:class="img-responsive"}
<br/>
<p align="justify">
Use the zoom slider to go to 1%, then you can see the whole cohort, and see that <i>TP53</i> is mutated in 36% of the queried samples. Some patients have multiple samples, hence the disparity in patient/sample numbers (top right).<br/>
Use the zoom slider or the mouse (click and hold to draw a small box over the OncoPrint) to zoom right in on a small number of samples. Hovering over individual samples gives more information such as the number of samples per patient (usually 1), the tumour study (top row) and the specific mutation(s) in a patient (bottom row). Some mutations have symbols attached to show any information we might have on what that mutation does. A blue target suggests it is <b>oncogenic</b> (<i>i.e.</i> important in cancer). Flames show it is a hotspot mutation (<i>i.e.</i> mutated in lots of people). Sometimes mutations can be tolerated by cells, or we simply don't know yet if that specific mutation impacts how the protein works.<br/><br/>
A potentially new concept for you will be <b>copy number changes</b>. This literally means how many copies of each gene a person has in every cell in their tumour. In healthy cells we should have 2 copies of every gene (as we are diploid organisms), with some exceptions when genes are found on the X/Y chromosomes. In cancers the genome can be very unstable and this can lead to deletions of some parts of the DNA, or amplifications where you get more than 2 copies of genes. This can have a huge impact on tumour biology!<br/>
With <i>TP53</i>, <b>nonsense</b> mutations are important (where cancers break the P53 proliferation brakes), but the same effect happens from a <b>deep deletion</b> of <i>TP53</i> (where the brakes get totally removed). This is a really important concept in cancer biology: different types of mutations can have the same impact on tumour biology. This can mean that the same drug can be given to people with different mutations, such as P53 mutations R175H, R248Q, R273C, <b>R273H</b>, R273L and R282W, as these mutations all change how P53 binds to DNA.<br/>
</p>
**EXPLAINER** [Which mutations are likely to have an impact in cancer?](#mutations-explained)<br/>
<br/><br/>

#### 4 *TP53* hotspots
<p align="justify">
Click on the <b>Mutations</b> tab above the plot. This will take you to a "Lollipop plot" for all mutations in <i>TP53</i> across TCGA. Lollipop plots where in the protein mutations occur (x axis) and how common they are (y axis - tallest lollipop sticks are most common).<br/>
You will see that mutations at R273 are the most common across the cohort. In green, red and blue blocks the plot shows the functional domains of P53 - the parts of the protein which perform its function. Just at a glance you'll see that there are more mutations towards the right hand side of the red domain and far fewer after the blue domain.<br/>
</p>
![TP53 Lollipop from TCGA PanCancer Atlas](/assets/images/ENTHUSE-03_TP53_Lollipop.jpg){:class="img-responsive"}
<br/>
<p align="justify">
Mutations happen all over <i>TP53</i>, but there are some patterns which allow us to understand the biology of what P53 is doing in tumours. Dark green points represent missense mutations, where the amino acid is changing. Black points are nonsense mutations, where the rest of the protein after the mutation is truncated (creates an early stop codon).<br/>
Use the Missense, Trucating, Inframe <i>etc.</i> table to the right of the plot to select only Driver mutations (hover your mouse to the right of 4213 and the word 'ONLY' will appear), then select only truncating mutations. What is the distribution of these mutations? Why are truncating mutations much less common at the end of the protein (after the blue domain)?<br/>
Remove that filter and do the same for Missense mutations - what are the major differences here?<br/>
</p>
![TP53 Lollipop from TCGA PanCancer Atlas](/assets/images/ENTHUSE-04_TP53_Lollipop_Missense.jpg){:class="img-responsive"}
<br/>
<p align="justify">
The red domain is where P53 binds to DNA. The blue domain is where P53 binds to other P53 proteins (working as a team). Missense mutations occur almost exclusively in these domains, impacting how P53 does its job. Nonsense/truncating mutations can happen pretty much anywhere, as long as they disrupt the function - that's why there are fewer after the blue domain, where a mutation is less likely to impact how the protein can function.<br/>
There is a lot of information here. Explore the plot (and the table below) further by applying different filters, hovering over mutations, annotations, domains <i>etc.</i> and following links to explore the known biology. After playing with the data (or if you get lost!), <a href="https://tinyurl.com/TP53-cBioPortal-TCGA-Lollipop">use this link to refresh the page back to the original plot and filters</a>.<br/>
</p>
<br/><br/>

#### 5 *TP53* mutations and prognosis
<p align="justify">
When trying to understand the impact of mutations, we can associate the presence of a mutation with patient prognosis (<i>i.e.</i> how long they are likely to survive after diagnosis). Again, this helps us better prioritise research and development funding.<br/><br/>
Click on the <b>Comparison/Survival</b> tab and then the small <b>Survival</b> tab which appears. This is called a Kaplan-Meier plot, and it shows the survival time of patients after diagnosis. The way to understand this plot is that patients with tumours with a <i>TP53</i> mutation (red, altered group) died more quickly than patients without a <i>TP53</i> mutation. If you're confused looking at the plot, visualise a line from 50% survival and see at how many months this line crosses the red and blue datasets. This is the median survival time for each group.<br/>
The x axis here goes to 30 years, but often cancer survival statistics are measured at 5 or 10 years. Try using the slider to see how survival rates change over time.<br/><br/>
cBioPortal helpfully gives us some warnings (blue and yellow boxes). Any statistical test has assumptions and confounding variables - features which could explain differences in the data which we are not being shown here. Before scrolling past the screenshot below (where <i>some</i> answers are), write down 4 confounding variables which could influence how we interpret this graph. Think about how those variables could influence your conclusions.<br/>
</p>
![TP53 survival from TCGA PanCancer Atlas](/assets/images/ENTHUSE-05_TP53_Survival.jpg){:class="img-responsive"}
<br/><br/>
You could have come up with some of the following ideas for confounding factors:
- Different responses in different cancers
- Different survival times in different cancers (<i>i.e.</i> treatment success)
- Difference in how advanced a tumour was when diagnosed
- Biological sex differences
- Severity of the *TP53* mutation in that person compared to other people in the altered group
- Other mutations in *TP53* mutated tumours (or in the unaltered group)
- Age at diagnosis (see how the plot shape changes with "Disease-specific" survival)
- Lifestyle factors such as weight, smoking status, economic status, employment history (think exposures)<br/>

<p align="justify">
Biology gives <b>noisy</b> data (unlike physics and chemistry) because there is natural variation between individuals, so there are always more confounding variables to consider. In statistics, if something is a big confounder, we can try and "control" for it during analysis. 
</p>
**EXPLAINER** [What is "big data" and is it useful?](#the-problem-with-big-data-explained)<br/>
<br/><br/>

#### 6 *TP53* mutations across different cancers
<p align="justify">
As we identified above, cancer type is a likely confounder in our understanding of <i>TP53</i> mutations. <i>TP53</i> is the most commonly mutated gene in cancer generally, but do you think it's mutated in all cancer types equally?<br/><br/>
Click on the <b>Cancer Types Summary</b> tab. As you can see the proportion of tumours with <i>TP53</i> mutations does vary a lot. In the plot, green is any kind of single nucleotide mutation, and the other colours represent bigger changes.<br/>Have a play with the y-axis value (does <b>Counts</b> change your interpretation? Is frequency or count more informative?), and changing from Cancer Study to Cancer Type to Cancer Type Detailed. What patterns can you see?<br/>
</p>
![TP53 cancer type split from TCGA PanCancer Atlas](/assets/images/ENTHUSE-06_TP53_CancerTypes.jpg){:class="img-responsive"}
<br/><br/>
<p align="justify">
We are now going to explore mutations in <i>TP53</i> and other genes in three cancer case studies. Return to the Cancer Study level summary - <a href="https://tinyurl.com/TP53-cBioPortal-TCGA-Summary">use this link to get back</a>.
</p>
<br/><br/>

#### 7A Sarcoma - multiple hits to the P53 pathway
<p align="justify">
On the Cancer Types Summary graph, hover over the Sarcoma bar on the chart and click on <b>Query this study for <i>TP53</i></b>. Now focused on this one study ,look at the OncoPrint, Mutations and Comparison/Survival tabs - what do you observe in Sarcoma compared to all cancers from TCGA?<br/>
Interestingly, eventhough <i>TP53</i> mutations are very common in Sarcoma, and the mutations are spread across the protein as before, there is no significant difference in survival.<br/><br/>
In biology, proteins rarely act on their own, instead forming <b>pathways</b> with other proteins and molecules to complete functions. Select the <b>Pathways</b> tab.<br/>
The top pathway is the eponymous P53 pathway, with <i>TP53</i> at the centre promoting (arrow) senescence (ageing) and apoptosis (cell death), and inhibiting (flat arrowhead) growth and proliferation.<br/>
</p>
![TP53 pathway](/assets/images/ENTHUSE-07_TP53_Pathway.jpg){:class="img-responsive"}
<br/>
<p align="justify">
The percentages show how frequently each gene in the pathway is mutated in Sarcoma. <i>TP53</i> mutations are common, but so are those in <i>CDKN2A</i> and <i>MDM2</i> - both of which negatively regulate P53. The question is, why mutate other genes in this pathway, particularly if they help to limit the regulatory role of P53? Is this what is actually happening? Let's find out.<br/><br/>
At the top of the page, under Sarcoma (TCGA, PanCancer Atlas), click on the pencil symbol next to TP53. Add MDM2 and CDKN2A to the list and <b>Submit Query</b>.
</p>
![TP53 pathway](/assets/images/ENTHUSE-08_TP53_MDM2_CDKN2A.png){:class="img-responsive"}
<br/>
<p align="justify">
Remember that <i>TP53</i> was mutated in 47% of Sarcoma cases, but the percentage of tumours with a mutation in at least 1 of <i>TP53</i>, <i>CDKN2A</i> and <i>MDM2</i> is now 74%. What does this suggest? Head to the OncoPrint tab again - what does this show?<br/>
</p>
![TP53 MDM2 CDKN2A mutual exclusivity](/assets/images/ENTHUSE-09_TP53_MDM2_CDKN2A_mutexcl.jpg){:class="img-responsive"}
<br/>
<p align="justify">
<i>MDM2</i> is almost never mutated in tumours where <i>TP53</i> is mutated. This is called <b>mutual exclusivity</b> - simply put, if you have a mutation in one, you don't in the other. Remember in the pathway that <i>MDM2</i> was seen to inhibit <i>TP53</i>, and the mutations we get in <i>MDM2</i> are all amplifications - generating <b>more</b> MDM2 protein, and having a bigger inhibtory effect on P53 activity.<br/>
The picture is more complicated with <i>CDKN2A</i> but there is little overlap with <i>MDM2</i> mutations, and most of the <i>TP53</i> overlap is when there is a missense, change/gain of function in P53 (dark green), rather than a truncating mutation (black).<br/>
The importance in tumour biology here is that in 74% of Sarcoma samples, the P53 pathway is being broken. This is happening in multiple ways with the cancer using different mutations. This means that in Sarcoma, we can't just use <i>TP53</i> mutations as a marker of altered P53 biology - we should use all three genes in tandem.<br/>
Now, go back to the Comparison/Survival tab, then the Survival tab. The altered group (a mutation in any one of the three genes) now has much worse survival than the unaltered group, reflecting the negative impact of breaking the P53 pathway, not just breaking P53 itself. <br/><br/>
<a href="https://tinyurl.com/TP53-cBioPortal-TCGA-Summary">Return to the Cancer Study summary graph</a> for the next case study.
</p>
<br/><br/>


#### Next sections (not done)
Cancer Types Summary
Cancer Type case studies:
- Sarcoma (TP53, MDM2, CDKN2A; mut excl - consequence of hitting same pathway in different ways)
- AML (TP53, NPM1; mut excl; very different mechanisms - road map idea)
- Head & Neck (differential clinical association with HPV - different disease initiations)<br/>
Finish on points about data use, and that these are real people.










<br/><br/><br/><br/><br/><br/><br/><br/><br/><br/><br/><br/><br/><br/><br/>

### Session walkthrough
<p align="justify">
The video below is a walkthrough of the exercises outlined above for the exploration of <i>TP53</i> mutations in over 10,000 samples of The Cancer Genome Atlas.<br/>
CURRENTLY A PLACEHOLDER VIDEO<br/>
</p>
{% include youtube.html id="1z03uo_kpQI" %}
<br/>
[Return to tutorial introduction.](#introduction)
<br/><br/><br/><br/><br/><br/><br/><br/><br/>

### Explainers
#### cBioPortal explained
<p align="justify">
Text holder.
</p>
[Return to session](#1-open-the-cbioportal-website)
<br/><br/><br/>

#### Sequencing explained
<p align="justify">
Text holder.
</p>
[Return to session](#1-open-the-cbioportal-website)
<br/><br/><br/>

#### The Cancer Genome Atlas explained
<p align="justify">
Text holder.
</p>
[Return to session](#2-select-the-pancancer-tcga-dataset)
<br/><br/><br/>

#### Gene and Protein naming explained
<p align="justify">
Text holder.
</p>
[Return to session](#2-select-the-pancancer-tcga-dataset)
<br/><br/><br/>

#### Mutations explained
<p align="justify">
Text holder.
</p>
[Return to session](#3-explore-the-tp53-oncoprint-plot)
<br/><br/><br/>

#### The Problem with Big Data explained
<p align="justify">
Text holder.
</p>
[Return to session](#5-TP53-mutations-and-prognosis)
<br/><br/><br/>
