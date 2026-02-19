---
layout: page
permalink: /courses/MSc_cBioPortal_Feb2026
---
<span style="font-size:1.6em;">**Introduction to cBioPortal (73M tutorial)**</span><br/>
![cBioPortal and Elixir logos](/assets/coursefiles/2024-03-21_cBioPortal/cBioPortal_Elixir.jpg){:class="img-responsive"}
<br/>
<p align="justify"><a href="https://www.cbioportal.org/">cBioPortal</a> is a publicly available online resource for cancer genomics, with omics and clinical data available for many cancer types. In this tutorial we will introduce cBioPortal - what's there, how to use it, and what you can learn from it.<br/></p>

<span style="font-size:1.4em;">**Learning objectives**</span><br/>
1. Recognise the applications and utility of cBioPortal for cancer research
2. Operate and explore the cBioPortal website to identify cancer data of interest
3. Complete two cancer biology problem-solving tasks using cBioPortal
4. Recognise the process for accessing and analysing cBioPortal data
<br/>

<span style="font-size:1.4em;">**Session outline**</span><br/>
15.00&nbsp;&nbsp;&nbsp;Introduction and Learning Objectives<br/>
15.05&nbsp;&nbsp;&nbsp;cBioPortal website demo<br/>
15.15&nbsp;&nbsp;&nbsp;Problem-solving tasks<br/>
15.40&nbsp;&nbsp;&nbsp;Recap and Further Resources<br/>
15.45&nbsp;&nbsp;&nbsp;Accessing and using underlying cBioPortal data<br/>
15.55&nbsp;&nbsp;&nbsp;Close
<br/>

<span style="font-size:1.4em;">**The Session**</span><br/>
<span style="font-size:1.2em;">**Introduction and Learning Objectives**</span><br/>
<object data="/assets/files/2026-02-24_cBioPortal_0102_Intro.pdf" width="700" height="450" type='application/pdf'></object>
<br/>

<span style="font-size:1.2em;">**cBioPortal website demo**</span><br/>
<p align="justify">During the session this is a live demonstration of the cBioPortal functionality, covering the available datasets, then using TCGA Bladder Cancer (Cell 2017) data to investigate: the summary and clinical dashboard tabs, a single- and then multi-gene query. Both gene queries are used to epxlore oncoprints, lollipop plots, plotting functionality and survival.<br/></p>
The video below is an indicative recording of a cBioPortal demo.
{% include youtube.html id="1z03uo_kpQI" %}
<br/>

<span style="font-size:1.2em;">**Problem solving tasks**</span><br/>
<p align="justify">Now you've seen what cBioPortal is and the kinds of data available, it's time to explore! In the guided examples below, follow the instructions and think about the questions posed - how is the data able to support your conclusions, and what caveats must you consider?<br/></p>

**Task 1 - the pan-cancer importance of *TP53* mutations**<br/>
<p align="justify">This first task covers <i>TP53</i> mutations across the 32 different cancer studies of The Cancer Genome Atlas, an enormous international study to better understand the diversity of genetic changes in over 10,000 cancer patients. <i>TP53</i> is one of the most commonly studied genes in cancer as it is known as the "Guardian of the Genome". It is also very commonly mutated, as you will see.<br/></p>

<p align="justify">
Open up any web browser and head to <a href="https://www.cbioportal.org/">www.cbioportal.org</a><br/><br/>
The home screen (pictured below) gives immediate access to results from sequencing data generated from thousands of cancer samples from all over the world. Scroll down the list. Each study represents potentially years of work from doctors, surgeons and nurses in hospitals to recruit patients, work with families and collect samples, and then from biomedical researchers in universities to process the samples, then generate, analyse, interpret and publish the results.<br/></p>
![cBioPortal homepage www.cbioportal.org](/assets/coursefiles/2023_ENTHUSE/ENTHUSE-01_cBioPortal_home.jpg){:class="img-responsive"}
<br/>
<p align="justify">
Return to the top of the cBioPortal homepage. We're going to work with The Cancer Genome Atlas (TCGA) pancancer study. "Pan" in this context just means "across lots of different cancers".<br/>
Click the Quick select link for <b>TCGA PanCancer Atlas Studies</b>. You have just loaded data for 32 different cancer studies and 10967 samples coming from 10528 different people! Next, click on the blue <b>Query By Gene</b> button, scroll to the bottom and type <b>TP53</b> into the Enter Genes box, then <b>Submit Query</b>.<br/>
</p>
<p align="justify">
Lots of colours! Don't panic. An OncoPrint is a big summary of mutation data - is there a mutation in a particular sample, and what type of mutation is it? In an OncoPrint each column represents a single sample. You have loads of samples, so everything is very squished.<br/>
You can largely ignore the three "Profiled for..." rows - these just indicate which sequencing technologies were used on those samples. For example, you'll see a large number of glioblastoma samples were only profiled for copy number changes, not mutations or structural changes.<br/></p>
![TP53 OncoPrint from TCGA PanCancer Atlas](/assets/coursefiles/2023_ENTHUSE/ENTHUSE-02_TP53_OncoPrint.jpg){:class="img-responsive"}
<br/>
<p align="justify">
Use the zoom slider to go to 1%, then you can see the whole cohort, and see that <i>TP53</i> is mutated in 36% of the queried samples. Some patients have multiple samples, hence the disparity in patient/sample numbers (top right).<br/>
Use the zoom slider or the mouse (click and hold to draw a small box over the OncoPrint) to zoom right in on a small number of samples. Hovering over individual samples gives more information such as the number of samples per patient (usually 1), the tumour study (top row) and the specific mutation(s) in a patient (bottom row). Some mutations have symbols attached to show any information we might have on what that mutation does. A blue target suggests it is <b>oncogenic</b> (<i>i.e.</i> important in cancer). Flames show it is a hotspot mutation (<i>i.e.</i> mutated in lots of people). Sometimes mutations can be tolerated by cells, or we simply don't know yet if that specific mutation impacts how the protein works.<br/><br/>
A potentially new concept for you will be <b>copy number changes</b>. This literally means how many copies of each gene a person has in every cell in their tumour. In healthy cells we should have 2 copies of every gene (as we are diploid organisms), with some exceptions when genes are found on the X/Y chromosomes. In cancers the genome can be very unstable and this can lead to deletions of some parts of the DNA, or amplifications where you get more than 2 copies of genes. This can have a huge impact on tumour biology!<br/>
With <i>TP53</i>, <b>nonsense</b> mutations are important (where cancers break the P53 proliferation brakes), but the same effect happens from a <b>deep deletion</b> of <i>TP53</i> (where the brakes get totally removed). This is a really important concept in cancer biology: different types of mutations can have the same impact on tumour biology. This can mean that the same drug can be given to people with different mutations, such as P53 mutations R175H, R248Q, R273C, R273H, R273L and R282W, as these mutations all change how P53 binds to DNA.<br/>
</p>
<p align="justify">
Click on the <b>Mutations</b> tab above the plot. This will take you to a "Lollipop plot" for all mutations in <i>TP53</i> across TCGA. Lollipop plots where in the protein mutations occur (x axis) and how common they are (y axis - tallest lollipop sticks are most common).<br/>
You will see that mutations at R273 are the most common across the cohort. In green, red and blue blocks the plot shows the functional domains of P53 - the parts of the protein which perform its function. Just at a glance you'll see that there are more mutations towards the right hand side of the red domain and far fewer after the blue domain.<br/>
</p>
![TP53 Lollipop from TCGA PanCancer Atlas](/assets/coursefiles/2023_ENTHUSE/ENTHUSE-03_TP53_Lollipop.jpg){:class="img-responsive"}
<br/>
<p align="justify">
Mutations happen all over <i>TP53</i>, but there are some patterns which allow us to understand the biology of what P53 is doing in tumours. Dark green points represent missense mutations, where the amino acid is changing. Black points are nonsense mutations, where the rest of the protein after the mutation is truncated (creates an early stop codon).<br/>
Use the Missense, Trucating, Inframe <i>etc.</i> table to the right of the plot to select only Driver mutations (hover your mouse to the right of 4213 and the word 'ONLY' will appear), then select only truncating mutations. What is the distribution of these mutations? Why are truncating mutations much less common at the end of the protein (after the blue domain)?<br/>
</p>
<details>
   <summary>Remove that filter and do the same for Missense mutations - what are the major differences here?</summary>
   Driver missense mutations are exclusively found in the functional domains, and not elsewhere in the protein.
</details><br/>
![TP53 Lollipop from TCGA PanCancer Atlas](/assets/coursefiles/2023_ENTHUSE/ENTHUSE-04_TP53_Lollipop_Missense.jpg){:class="img-responsive"}
<br/>
<p align="justify">
The red domain is where P53 binds to DNA. The blue domain is where P53 binds to other P53 proteins (working as a team). Missense mutations occur almost exclusively in these domains, impacting how P53 does its job. Nonsense/truncating mutations can happen pretty much anywhere, as long as they disrupt the function - that's why there are fewer after the blue domain, where a mutation is less likely to impact how the protein can function.<br/>
There is a lot of information here. Explore the plot (and the table below) further by applying different filters, hovering over mutations, annotations, domains <i>etc.</i> and following links to explore the known biology. After playing with the data (or if you get lost!), <a href="https://tinyurl.com/TP53-cBioPortal-TCGA-Lollipop">use this link to refresh the page back to the original plot and filters</a>.<br/>
</p>
<p align="justify">
When trying to understand the impact of mutations, we can associate the presence of a mutation with patient prognosis (<i>i.e.</i> how long they are likely to survive after diagnosis). Again, this helps us better prioritise research and development funding.<br/><br/>
Click on the <b>Comparison/Survival</b> tab and then the small <b>Survival</b> tab which appears. This is called a Kaplan-Meier plot, and it shows the survival time of patients after diagnosis. The way to understand this plot is that patients with tumours with a <i>TP53</i> mutation (red, altered group) died more quickly than patients without a <i>TP53</i> mutation. If you're confused looking at the plot, visualise a line from 50% survival and see at how many months this line crosses the red and blue datasets. This is the median survival time for each group.<br/>
The x axis here goes to 30 years, but often cancer survival statistics are measured at 5 or 10 years. Try using the slider to see how survival rates change over time.<br/><br/>
cBioPortal helpfully gives us some warnings (blue and yellow boxes). Any statistical test has assumptions and confounding variables - features which could explain differences in the data which we are not being shown here. Before scrolling past the screenshot below (where <i>some</i> answers are), write down 4 confounding variables which could influence how we interpret this graph. Think about how those variables could influence your conclusions.<br/>
</p>
<details>
   <summary>Are your confounding variables in this list? Did you think of others?</summary>
   Different responses in different cancers<br/>
   Different survival times in different cancers (<i>i.e.</i> treatment success)<br/>
   Difference in how advanced a tumour was when diagnosed<br/>
   Biological sex differences<br/>
   Severity of the <i>TP53</i> mutation in that person compared to other people in the altered group<br/>
   Other mutations in <i>TP53</i> mutated tumours (or in the unaltered group)<br/>
   Age at diagnosis (see how the plot shape changes with "Disease-specific" survival)<br/>
   Lifestyle factors such as weight, smoking status, economic status, employment history (think exposures)<br/>
</details><br/>
![TP53 survival from TCGA PanCancer Atlas](/assets/coursefiles/2023_ENTHUSE/ENTHUSE-05_TP53_Survival.jpg){:class="img-responsive"}
<br/><br/>

<p align="justify">
Biology gives <b>noisy</b> data because there is natural variation between individuals, so there are always more confounding variables to consider. In statistics, if something is a big confounder, we can try and "control" for it during analysis. <br/>
As we identified above, cancer type is a likely confounder in our understanding of <i>TP53</i> mutations. <i>TP53</i> is the most commonly mutated gene in cancer generally, but do you think it's mutated in all cancer types equally?<br/><br/>
Click on the <b>Cancer Types Summary</b> tab. As you can see the proportion of tumours with <i>TP53</i> mutations does vary a lot. In the plot, green is any kind of single nucleotide mutation, and the other colours represent copy number changes changes such as deletions.<br/>Have a play with the y-axis value (does <b>Counts</b> change your interpretation? Is frequency or count more informative?), and changing from Cancer Study to Cancer Type to Cancer Type Detailed. What patterns can you see?<br/>
</p>
![TP53 cancer type split from TCGA PanCancer Atlas](/assets/coursefiles/2023_ENTHUSE/ENTHUSE-06_TP53_CancerTypes.jpg){:class="img-responsive"}
<br/><br/>
<p align="justify">
In general terms, <b>solid</b> cancers (rather than blood cancers) with a high mutational burden tend to have <i>TP53</i> mutations - see the lung, esophageal, colorectal and bladder datasets, compared to Acute Myeloid Leukaemia or B cell lymphoma. As this is biology, there are some that bend this observation - particularly the two types of kidney cancer found in TCGA. Often, general observations hide features which are cancer-specific (and even specific to subtypes within those cancers).<br/><br/>
This tutorial has given you an introduction to some of the power cBioPortal has for cancer research. It makes large data (usually the sole focus of bioinformaticians) accessible to laboratory scientists, students, clinicians and the general public. But, there is so much more you can find out.<br/>
Try hovering your mouse over one of the cancer type bars in the Cancer Types Summary plot and then select the study name from the box which appears. This takes you to the cBioPortal dashboard for this study - there is so much data here! Take a look about (<a href="https://tinyurl.com/TP53-cBioPortal-TCGA-Summary">using this link to return to the Cancer Study summary graph if needed</a>) and see how much information there is!<br/><br/></p>

**Task 1A - the pan-cancer importance of *TP53* mutations**<br/>
<p align="justify">
We are now going to explore mutations in <i>TP53</i> and other genes in a sarcoma case study. By looking in multiple cancers we can see whether or not <i>TP53</i> mutations are common in all cancers, or whether there is variability. As we're working with biological data, there is a <b>lot</b> of variability!<br/><br/>
Start to think about why that might be. P53 (the protein name for the product from the <i>TP53</i> gene) regulates cell proliferation and cell death. Tumours need to escape P53 activity. But is mutating P53 the only way for a tumour to do this? And is regulating proliferation and cell death the only way for a tumour to form? Clearly not, because not every tumour has a <i>TP53</i> mutation, and tumours have a large number of ways to evade detection and eradication in the body.<br/><br/>
In his epic (third) review of the <i>Hallmarks of Cancer</i> published in 2022 (<a href="/assets/files/2023-06_ENTHUSE_Hallmarks_Hanahan2022.pdf" download>download the PDF here</a>), Professor Douglas Hanahan outlines the different physiological features a cancer must overcome (image below). These are all areas where cancers can effect change on the body in order to survive. Remember, a cancer is under the pressures of evolution in just the same way as an infectious virus, a herd of giraffes, or us - just on a different scale.<br/>
</p>
![Hanahan hallmarks 2022](/assets/coursefiles/2023_ENTHUSE/ENTHUSE-12_2022_Hallmarks_of_Cancer.jpg){:class="img-responsive"}
<br/><br/>

<p align="justify">
<a href="https://tinyurl.com/TP53-cBioPortal-TCGA-Summary">Return to the Cancer Study summary graph.</a><br/><br/>
On the Cancer Types Summary graph, hover over the Sarcoma bar on the chart and click on <b>Query this study for <i>TP53</i></b>. Now focused on this one study, look at the OncoPrint, Mutations and Comparison/Survival tabs.
<details>
   <summary>What do you observe in Sarcoma compared to all cancers from TCGA?</summary>
   Interestingly, eventhough <i>TP53</i> mutations are very common in Sarcoma, and the mutations are spread across the protein as before, there is no significant difference in survival.
   </details><br/>
In biology, proteins rarely act on their own, instead forming <b>pathways</b> with other proteins and molecules to complete functions. Select the <b>Pathways</b> tab.<br/>
The top pathway is the eponymous P53 pathway, with <i>TP53</i> at the centre promoting (arrow) senescence (ageing) and apoptosis (cell death), and inhibiting (flat arrowhead) growth and proliferation.<br/>
</p>
![TP53 pathway](/assets/coursefiles/2023_ENTHUSE/ENTHUSE-07_TP53_Pathway.jpg){:class="img-responsive"}
<br/>
<p align="justify">
The percentages show how frequently each gene in the pathway is mutated in Sarcoma. <i>TP53</i> mutations are common, but so are those in <i>CDKN2A</i> and <i>MDM2</i> - both of which negatively regulate P53. The question is, why mutate other genes in this pathway, particularly if they help to limit the regulatory role of P53? Is this what is actually happening? Let's find out.<br/><br/>
At the top of the page, under Sarcoma (TCGA, PanCancer Atlas), click on the pencil symbol next to TP53. Add MDM2 and CDKN2A to the list and <b>Submit Query</b>.
</p>
![TP53 pathway](/assets/coursefiles/2023_ENTHUSE/ENTHUSE-08_TP53_MDM2_CDKN2A.png){:class="img-responsive"}
<br/>
<p align="justify">
Remember that <i>TP53</i> was mutated in 47% of Sarcoma cases, but the percentage of tumours with a mutation in at least 1 of <i>TP53</i>, <i>CDKN2A</i> and <i>MDM2</i> is now 74%. 
<details>
   <summary>What does this suggest?</summary>
   The increased percentage of affected cases suggests people are not commonly mutated in all three genes at the same time.
   </details><br/> 
<details>
   <summary>Head to the OncoPrint tab again - what does this show?</summary>
   As we thought - cases with a mutation in <i>TP53</i> very rarely have mutations in the other two genes.
   </details><br/>
</p>
![TP53 MDM2 CDKN2A mutual exclusivity](/assets/coursefiles/2023_ENTHUSE/ENTHUSE-09_TP53_MDM2_CDKN2A_mutexcl.jpg){:class="img-responsive"}
<br/>
<p align="justify">
<i>MDM2</i> is almost never mutated in tumours where <i>TP53</i> is mutated. This is called <b>mutual exclusivity</b> - simply put, if you have a mutation in one, you don't in the other. Remember in the pathway that <i>MDM2</i> was seen to inhibit <i>TP53</i>, and the mutations we get in <i>MDM2</i> are all amplifications - generating <b>more</b> MDM2 protein, and having a bigger inhibtory effect on P53 activity.<br/>
The picture is more complicated with <i>CDKN2A</i> but there is little overlap with <i>MDM2</i> mutations, and most of the <i>TP53</i> overlap is when there is a missense, change/gain of function in P53 (dark green), rather than a truncating mutation (black).<br/><br/>
The importance in tumour biology here is that in 74% of Sarcoma samples, the P53 pathway is being broken. This is happening in multiple ways with the cancer using different mutations. This means that in Sarcoma, we can't just use <i>TP53</i> mutations as a marker of altered P53 biology - we should use all three genes in tandem.<br/>
Now, go back to the Comparison/Survival tab, then the Survival tab. The altered group (a mutation in any one of the three genes) now has much worse survival than the unaltered group, reflecting the negative impact of breaking the P53 pathway, not just breaking P53 itself. <br/><br/>
The video below is a walkthrough and explainer for this case study.<br/>
</p>
{% include youtube.html id="vP5pLJmf5io" %}
<br/>

**Task sum up**<br/>
<p align="justify">
This <i>TP53</i> study should have given you a large exposure to what cBioPortal can show you. Think about how you can use this in your own research.<br/>
Perhaps you have done an experiment in the lab and you want to know whether a particularly gene is important in the relevant cancer - cBioPortal can help you do this.<br/>
You can also start with the cBioPortal data for a cancer (or cancers) of interest and approach it as a data question.
<br/><br/></p>

<span style="font-size:1.2em;">**Recap and Further Resources, Accessing and using underlying cBioPortal data**</span><br/>
<object data="/assets/files/2025-02-25_cBioPortal_0304_LO-Recap_Resources_Downloads.pdf" width="700" height="450" type='application/pdf'></object>
<br/><br/>

<p align="justify">This tutorial was first delivered at The University of York on 27th February 2024 for the 59M and 68M MSc students. The material was written and delivered by Dr Andrew Mason, adapting previous material supported by <a href="https://elixiruknode.org/">Elixir-UK</a> and <a href="https://www.yorkagainstcancer.org.uk/">York Against Cancer</a>.<br/></p>