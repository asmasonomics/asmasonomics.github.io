---
layout: page
permalink: /courses/Intro_cBioPortal_Jan2023
---
<span style="font-size:1.6em;">**Introduction to cBioPortal**</span><br/>
![cBioPortal and Elixir logos](/assets/images/cBioPortal_Elixir.jpg){:class="img-responsive"}
<br/>
<p align="justify"><a href="https://www.cbioportal.org/">cBioPortal</a> is a publicly available online resource for cancer genomics, with omics and clinical data available for many cancer types. In this 1 hour session we will introduce cBioPortal - what's there, how to use it, and what you can learn from it. We will then have some led problem solving activities to get to grips with the website, and finish with tips and advice for how you can use the data, and how student projects could make use of it.<br/></p>

<span style="font-size:1.2em;">**Learning objectives**</span><br/>
1. Recognise the applications and utility of cBioPortal for cancer research
2. Operate and explore the cBioPortal website to identify cancer data of interest
3. Complete two cancer biology problem-solving tasks using cBioPortal
4. Recognise the process for accessing and analysing cBioPortal data
<br/>

<span style="font-size:1.2em;">**Session outline**</span><br/>
12.00&nbsp;&nbsp;&nbsp;Introduction and Learning Objectives<br/>
12.10&nbsp;&nbsp;&nbsp;cBioPortal website demo<br/>
12.20&nbsp;&nbsp;&nbsp;Problem-solving tasks<br/>
12.45&nbsp;&nbsp;&nbsp;Recap and Further Resources<br/>
12.50&nbsp;&nbsp;&nbsp;Accessing and using underlying cBioPortal data<br/>
13.00&nbsp;&nbsp;&nbsp;Close
<br/>

<span style="font-size:1.2em;">**The Session**</span><br/>
**Introduction and Learning Objectives**<br/>
<a href="/assets/files/2023-01-11_cBioPortal_01_Introduction-Learning-Objectives.pdf" download>PDF of the introductory slides available to download here</a>.

**cBioPortal website demo**<br/>
Live demonstration. Representative video to follow.

**Problem-solving tasks**<br/>
<p align="justify">You probably have time to complete 2 of the 3 tasks below, so pick those most interesting to you! This web page will remain live following the course if you want to complete all the taks, and for your future reference.<br/>The point here is to explore the data, not just rush through the questions. The questions are to lead your searches, but can you think of other data to investigate?<br/>You could take a look at TCGA's pan-cancer analysis of 32 tumour types, but focus on summary comparisons as cBioPortal can be slow >2000 samples...!<br/></p>
[Task 1](#task-1) - Exploration of the METABRIC breast cancer dataset<br/>
[Task 2](#task-2) - Exploration of two AML datasets<br/>
[Task 3](#task-3) - Exploration and comparison of two kidney cancers<br/>

#### Task 1
<p align="justify">In this task you will explore the METABRIC study, one of the largest cancer cohorts. Breast cancer is one of the best served cancers in terms of genomic resources. Explore METABRIC and consider the following questions, but also take the time to look at other sets including TCGA, a cohort built on large, aggressive tumours (across all tumour types). Can you find any breast cancer cell lines which well represent patient data?<br/></p>
1. Are there differences in survival or other clinical data between the 3-gene classifier subtypes?
2. How does ER status survival change over time (5-year vs 10-year vs 20-year)?
3. What gene expression or mutation differences can you see between ER+ and ER- BRCA?
4. Which are the most commonly mutated, likely cancer-related genes? Of the top 10, which are consisteny with being tumour supressor genes (TSGs), and which as oncogenes? (Think about the mutation types). Are there hotspots in these with clinical relevance?
5. Consider *MUC16*. It is not included in the 'likely cancerous' list, but is heavily mutated. What data from the lollipop mutation plots make it unlikely related to oncogenesis?
6. How does tumour mutational burden (TMB) compare to the demo example from the TCGA bladder cancer (2017) cohort? How could this impact our study of mutations in METABRIC?
7. Looking across BRCA studies, are there any male samples? If so, how do they compare? Can you compare them confidently?

#### Task 2
<p align="justify">Despite being a rare cancer, Acute Myeloid Leukemia is widely studied, due to its broad age range and very poor survival. Choose the TCGA (NEJM 2013) and OHSU datasets, and select the samples with mutation information.<br/></p>
1. Look at some of the clinical data - how do survival, other clinical features and mutation status vary between male and female, or depending on age of diagnosis?
2. How does the mutational burden data compare to solid cancers such as bladder and colorectal? How is this reflected in the most mutated genes? 
3. Select the cancer-likely genes (filtered) with more than 8% frequency and submit the query. Do the 2 cohorts look similar, or are there major differences?
4. Are any of these top 11 most mutated genes co-occuring or mutually exclusive?
5. Do the lollipop diagrams suggest these genes are tumour suppressors, or oncogenes? Are there hotspots, any with clinical actions?
6. I would say the *NPM1* plot is a bit strange. Highly mutated, but do you think its likely these mutations are impacting cancer development? Take a look at the UniProt link (NPM_HUMAN) and scroll to the protein domains section. Does this change your opinion?
7. Do any of these most mutated genes exhibit differences in survival, or other clinical features?

#### Task 3
<p align="justify">Often the public (and us researchers) conflate multiple cancers together from the same organ. cBioPortal includes TCGA data from both kidney renal papillary cell carcinoma, and clear cell carcinoma. In this task we will explore both, and use some of the limited group functionality comparisons - this can be a bit clunky.<br/></p>
1. Select TCGA's Firehose Legacy Kidney Renal Paillary Cell Carcinoma. As with the previous tasks, explore the most mutated genes, sex imbalances.
2. Given the male imbalance, the mutations in *AR* (androgen receptor) could be interesting. Filter for male. Filter for *AR* mutations and calculate whether expression of *AR* differs, how *AR* is mutated and what other genes are co-occuring or mutually exclusive.
3. Change your query, and have a similar exploration of TCGA's Kidney clear clear cell data. Do you notice any immediate differences?
4. Return to the cBioPortal homepage (or alter the query) and select both datasets. In the cancer types, use the compare groups functionality to compare the two cohorts. What are the first obvious differences? Are there clinical differences (beyond tissue and other coding classifiers)?
5. Both types of kidney cancer appear more common in males. Select only female patient samples and redo the comparison. What has changed?
6. Still with the female samples, take the 5 most mutated likely cancer genes. Do the 2 cohorts match well? Are there co-occurences or mutually exclusive genes?
7. The lollipop plot for *SETD2* is curious, consistent with tumour suppressor mutations, but also potentially oncogenic activation of the SET domain. Modify the query to just *SETD2* and just the renal clear cell carcinoma dataset - can't compare expression data between datasets readily. Explore the plots of expression when split by copy number or mutation status. Are there many differences?

<br/>
**Recap and Further Resources, Accessing and using underlying cBioPortal data**<br/>
<a href="/assets/files/2023-01-11_cBioPortal_0304_LO-Recap_Resources_Downloads.pdf" download>PDF of final section slides available to download here</a>.
<br/><br/>

<span style="font-size:1.2em;">**Feedback form**</span><br/>
<p align="justify">This was the first time running the Introduction to cBioPortal course. We are hoping to repeat the course in future, and to develop it into an Elixir online resource for broader use.
Hopefully you have enjoyed the workshop structure and content, but there is always room to improve. We would really appreciate if you could complete <a href="https://forms.gle/oWjnkAU6tGP2W91H7">this brief feedback form</a>. This link will be active until the end of the week.<br/></p>

<br/><br/>
<p align="justify">This course was first delivered over 1 hour at The University of York on 11th January 2023. The material was written and delivered by Dr Andrew Mason, with support from <a href="https://elixiruknode.org/">Elixir-UK</a>.<br/></p>