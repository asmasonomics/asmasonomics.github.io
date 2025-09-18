---
layout: page
permalink: /courses/BSc_dissertation_2025
---

![project banner](/assets/coursefiles/2025_BSc_dissertation_project/BSc_multiomic_BLCA_banner.png){:class="img-responsive"}

<span style="font-size:1.6em;">**Mason Lab BSc project - Multiomic analysis of bladder cancer**</span><br/>

<p align="justify">
This site contains helpful information for completing your 28H final year BSc dissertation project in the Mason Lab. <a href="https://vle.york.ac.uk/ultra/courses/_114373_1/outline">General information on the module is available on the VLE</a> - check this regularly.<br/>
</p>

### Expectations, meetings and assessments
<p align="justify">
Across the year you are expected to spend 1.5-2 days per week on your project, balancing this 40 credit module with your other module topics. As this is a data project, you can be very flexible with when you are working. Remember, this is an individual project (and must be completed and submitted independently), but it is likely the other students who also selected this topic will face similar issues, so support each other! You may enjoy arranging to meet in a computer lab or in the library at the same time in small groups. Whatever works for you.<br/><br/>
28H does have some workshops, so make sure you check the Module Planner on the <a href="https://vle.york.ac.uk/ultra/courses/_114373_1/outline">VLE site</a> and look out for any announcements from the module organiser. There is a <a href="https://docs.google.com/document/d/1I23Xv6QPa_hg3XFhiO3tt6Ghl3S379HwJx9lddAYUko">useful checklist of expectations available here</a>.<br/>
During semester 1, we will have our group supervision meetings on <b>Wednesdays 2-3 in B/M/049</b>. Some particular dates to note are tabled below:<br/><br/>
</p>

| Week | Date | |
| --- | --- | --- |
| S1 W2 | W 01/10/25 | First group supervision meeting |
| S1 Consolidation Week | W 29/10/25 | No supervision meeting |
| S1 Week 6 | T 04/11/25 | Meeting is on Tuesday at 1pm (not Wednesday) in B/M/049 | 
| S1 Week 8 | W 19/11/25 | Extended session to include formative presentations | 
| S1 Week 11 | W 10/12/25 | Last meeting of semester |
| S2 Week 1 | w/c 09/02/26 | First supervision meeting of semester 2 (day/time TBC) | 
| S2 Week 6 | w/c 16/03/26 | Final supervision meeting of project (day/time TBC) | 

<p align="justify">
As your project develops I will set aside an open office hour where I will be available in case you want to drop in to discuss your project. This will be on a first come, first served basis.<br/><br/>

The <b>main assessment is your 4000 word scientific report</b>. Deadline: Tuesday 5th May at 12 noon (S2 W11).<br/><br/>
There is also the formative oral presentation in semester 1 (S1 W8). There are multiple opportunities for feedback on your writing (S2 W0, S2 W8) - these are for you to make the most of, and are not compulsory.
<br/><br/>
</p>

### Scope and aims of your project
<p align="justify">
Bladder cancer is an incredibly diverse disease, characterised by a very high mutational burden and the highest number of different driver genes among solid cancers. The challenge is to understand the drivers in each patient’s cancer so that we can select, repurpose or develop more personalised treatment.<br/>
In this project you will use different types of sequencing data generated from 408 advanced bladder cancers included in The Cancer Genome Atlas (TCGA), an international programme built to understand over 30 different cancer types from over 10,000 patients. Whilst the power of TCGA was having multiple data types, most studies have only looked at one. In your project you will combine datasets to better understand the impact of mutations and gene expression changes on bladder cancer biology. The overall aim of my lab is to better understand existing bladder cancer subtypes, or to find new subtypes which could have a route to clinic.<br/><br/>
So, how about <i>your</i> project? You will do the following:
<ol>
<li>Familiarise yourself with muscle-invasive bladder cancer literature</li>
<li>Devise a relevant, sensible and testable investigation where you can compare different types of muscle-invasive bladder cancers using the multiple data types you have available to you</li>
<li>Perform a systematic data investigation to answer your question</li>
<li>Critically review your results and determine how this fits within the wider literature - what have we learned and is it clinically relevant?</li>
<li>Consider and plan the future work you would propose to further understand your results and their relevance (incl. wet lab work, further data analysis, clinical follow up <i>etc.</i>)</li>
</ol>
We would expect you to complete all your analysis using the skills you have developed in R. This means you can use RStudio on a personal laptop, or using university-managed machines - whatever works best for you. Just make sure you keep saving and backing up your data and scripts.<br/><br/>
Coming up with your own study might feel daunting just now, but don't worry - there are some useful starting points below including some key papers. If you are doing the Cancer cell & molecular biology module, or the Genomics module, you should get some useful context and help from lectures and workshops - don't worry if you don't take either of these. A minimal project would be to identify tumours with/without a mutation in a relevant cancer gene, and then to use RNA sequencing data to compare the gene expression profiles of these cancers, and to try and work out what that means. We have lots of different data (and there are >400 cancers in the dataset) so the project can grow, develop, extend in lots of different ways - the best approach is to try and find a question where you'd be interested in finding out the answer.<br/><br/>
Something you won't have maybe thought about too much before is the ethics of working with human genomic data. We'll speak about this during our sessions, but your project work is covered under <a href="/assets/coursefiles/2025_BSc_dissertation_project/MoA_AM202410.pdf" download>an existing ethics approval from the Biology Ethics Committee.</a>
<br/><br/>
</p>

### Introduction to bladder cancer
<p align="justify">
To get you started I've recorded a mini lecture on muscle invasive bladder cancer. This is enough to get you going, but you'll want to dig deeper in areas relevant to your project (including looking at studies from other cancer types). Remember - what I've said is my scientific assessment based on the literature and my own research experience - your critique is welcome and expected!<br/>
</p>
{% include youtube.html id="e3pWczidol0" %}
Download slides in <a href="/assets/coursefiles/2025_BSc_dissertation_project/2025_BSc_Intro-to-BLCA.pptx" download>pptx</a> or <a href="/assets/coursefiles/2025_BSc_dissertation_project/2025_BSc_Intro-to-BLCA.pdf" download>PDF</a> formats.
<br/><br/>

#### Relevant (big) starting papers
<p align="justify">
<a href="https://doi.org/10.1016/j.cell.2017.09.007">Robertson <i>et al</i>, 2017, Cell</a> - the paper which describes the full TCGA MIBC cohort<br/>
<a href="https://doi.org/10.1016/j.eururo.2019.09.006">Kamoun <i>et al</i>, 2020, European Urology</a> - the paper which introduces the MIBC "consensus classifier"<br/>
<a href="https://doi.org/10.1186/s13073-020-00781-y">Shi <i>et al</i>, 2020, Genome Medicine</a> - the paper which describes bladder cancer hotspot mutations, particularly those created by APOBEC mutagenesis<br/>
<a href="https://doi.org/10.1038/s41388-022-02235-8">Baker <i>et al</i>, 2022, Oncogene</a> - paper from our unit which attempts to link viral infection to bladder cancer<br/>
<a href="https://doi.org/10.1016/j.cell.2017.09.042">Martincorena <i>et al</i>, 2017, Cell</a> - paper identifies recurrent driver genes across all cancers in TCGA - see Figure 2<br/>
<a href="https://doi.org/10.1016/j.jmoldx.2024.08.005">Cotillas <i>et al</i>, 2024, The Journal of Molecular Diagnostics</a> - an alternative view to bladder cancer classification (all disease stages)
<br/><br/>
</p>

### The data
<p align="justify">
<a href="https://drive.google.com/drive/folders/1NW0C4x3VF04U5HYUrsQrNHocIMPWQJ_q?usp=sharing">This google drive </a>contains major data files from The Cancer Genome Atlas muscle invasive bladder cancer cohort. We have pre-processed some of these to make them less difficult to work with, but some of them are still very large! You can and should work with these in R. You can also open these in notepad (but use notepad++ because this has some actual functionality).<br/>
Most of the data files are tab spaced value (tsv) files which means you can open and process them easily (even in Excel if you want a quick look). The .maf and .bed format files are tsv format files as well, but they are organised in a partciular way with particular columns in a particular order. A description of each file is given in the table below:<br/><br/>
</p>

| Filename | Technology | Description |
| --- | --- | --- |
| CNA_gene-level-CN.tsv | copy number array  |  number of copies of each gene in each tumour (could represent segmental or casette alterations) |
| gc47_gene_locations.bed | annotation file | GRCh38 gene locations across the genome | 
| methylation_CpG-array.tsv |  5mC methylation array | methylation intensity beta values between 0 (low) and 1 (high) | 
| miRNA_counts.tsv | miRNA sequencing (short RNAseq) | raw read counts per miRNA gene |
| miRNA_CPMs.tsv | miRNA sequencing (short RNAseq) | counts per million to control for differences in sequencing depth | 
| mRNA_gc47-counts.tsv | mRNA sequencing (polyA RNAseq) | raw read counts per gene feature |
| mRNA_gc47-TPMs.tsv | mRNA sequencing (polyA RNAseq) | transcripts per million to control for depth and feature length |
| mRNA_gc47-TPMs_ConsensusClassifier.tsv | sample annotation | Consensus classification of each tumour according to Kamoun et al (2020) |
| mRNA_gc47-TPMs_LundTax2023Classifier.tsv | sample annotation | Lund classification of each tumour according to Cotillas et al (2024) |
| TCGA-BLCA-clinical-metadata.tsv | patient metadata | relevant metadata for patients and their cancer |
| WGS_coding-regions-only.maf | whole genome sequencing | mutations in tumours filtered to just include mutations in coding regions |
| WGS_mutation-in-at-least-4-patients.maf | whole genome sequencing | mutations in tumours filtered to only include mutations found in 4 or more patients |
| WXS.maf | whole exome sequencing | mutations in tumours identified by only profiling coding regions |

<p align="justify">
I've provided these datafiles so you can do your own specific analysis, including making your own versions of figures. There are places to get data summaries for this cohort, even some top-level analysis and figures. But these alone will not be enough for a good project mark. By coding it yourself you can ask good questions and get to grips with the data - allowing you to critique how and when it can (or should) be used (and when it really shouldn't).<br/><br/>
Before jumping into the data analysis, do some reading to get a better feel for the cancer biology and the technologies I've talked about. Then you should use <a href="https://www.cbioportal.org/">cBioPortal</a> to do some initial playing with the data (in a browser). I've written a course (~2hours) on how to use cBioPortal, including with demos, problem solving tasks to get used to using the website, and then a worked coding example in R to demonstrate how to download and work with this type of data. This will be very (very) useful for your project - mix up your reading by <a href="https://asmasonomics.github.io/courses/Intro_cBioPortal_Mar2024">completing my Intro to cBioPortal course in your own time.</a>
<br/><br/>
</p>

### Getting started with your analysis
<p align="justify">

<br/><br/>
</p>

