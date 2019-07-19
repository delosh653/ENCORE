# Welcome to ENCORE!

<p align="center">
<img src="ENCORE Shiny App/www/ont_map.PNG" width="275" />  <img src="ENCORE Shiny App/www/ont_nav.PNG" width="275" />  <img src="ENCORE Shiny App/www/group_comp.PNG" width="275" />
</p>

This is the second step in the PAICE (Pipeline for Amplitude Integration of Circadian Exploration) Suite! This suite of tools provides high-throughput applications for circadian, ultradian, and infradian rhythms. The first step, ECHO, can be found [here](https://github.com/delosh653/ECHO), and must be run before using this application.

## README Outline

* Overview
* Use and First-Time Set-Up Instructions
* ENCORE Features
* Contact Information and Bug Reporting
* FAQ

## Overview

Note: currently in beta testing. Since this app is currently in development, please check here for updates, especially for an ENCORE citation.

ENCORE (ECHO Native Circadian Ontological Ryhtmicity Explorer) is an R-powered Shiny application to navigate and understand the gene function of amplitude change categories for circadian rhythms, using [ECHO](https://github.com/delosh653/ECHO) output.  

All images created by ENCORE using data from [*Harmonics of Circadian Gene Transcription in Mammals*](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000442) by M.E.Hughes, et al. (2009).

## Use and First-Time Set-Up Instructions

Thank you for downloading ENCORE (ECHO Native Circadian Ontological Rhythmic Explorer)! ENCORE is an app designed to help you navigate and understand the function of amplitude change categories for circadian rhythms with gene ontologies. This guide will lead you in first time set-up and use. Pictures have been provided for ease of use, using Windows 10, in the files ENCORE_README.docx and ENCORE_README.pdf, found above. A double asterisk indicates the step has an explanation below, and a tilde indicates the step is first-time set up only.

Steps: 
1.	** ~ Download [Firefox](https://www.mozilla.org/en-US/firefox/new/) or [Chrome](https://www.google.com/chrome/browser/desktop/index.html) and make it your default browser.

2.	~ Download the protein links for your organism from [STRING](https://string-db.org/cgi/download.pl). The organisms that are available in this application are:

| Organism  | Taxonomy Number |
| ------------- | ------------- |
| Mus musculus  | 10090  |
| Homo sapiens  | 9606  |
| Neurospora crassa  | 5141  |
| Drosophila melanogaster  | 7227  |
| Anopheles gambiae  | 7165  |
| Saccharomyces cerevisiae  | 4932  |
| Escherichia coli K12 substr MG1655  | 511145  |

Enter your organism’s name in the ‘choose an organism’ box and click update. You should then see the taxonomy number preceding the file names. Download the file called TAX#.protein.links.v11.0.txt.gz.

3.	~ Once you’ve downloaded that file, extract the text file. (You may need to download an alternative unzip program that can unzip .gz files, such as [7Zip](https://www.7-zip.org/).) Place that text file, unaltered, in the ‘links’ folder of the ‘ENCORE Shiny App’ Folder.

4.	~ [Download R](https://www.r-project.org/), if you do not already have it. 

5.	~ [Download RStudio](https://www.rstudio.com/products/rstudio/download/), if you do not already have it (RStudio Desktop is sufficient).

6.	Plug in your computer, if it is not plugged in already, and open RStudio.

7.	~ Copy and paste the following text into the console window (bottom right window of the RStudio Session), then press enter:

```r
install.packages("rstudioapi")
install.packages("shiny")
install.packages("ggplot2")
install.packages("r2d3")
install.packages("data.table")
install.packages("jsonlite")
install.packages("stringr")
install.packages("igraph")
```

Copy each of the following lines in the console separately and press enter for each:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("STRINGdb")
BiocManager::install("AnnotationHub")
BiocManager::install("mygene")
BiocManager::install("topGO")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Ag.eg.db")
BiocManager::install("org.Dm.eg.db")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("org.EcK12.eg.db")
BiocManager::install("org.Sc.sgd.db")
```

This will install these packages (a set of functions that this application uses) onto your computer. This may ask for your input, so just say no to the questions asked. If you run into errors saying “no,” just say yes instead. Note: this may take some time.

8.	Open encore_app.R, which should be included in the .zip file you downloaded and also contained this README. It should open in the top left window of your RStudio session.

9.	In the top right corner of the encore_app.R window, you should see the button, “Run App”. Click on the small downwards arrow next to it and choose “Run External”. 

10.	Now click “Run App”. This should open the ENCORE application in your now default browser window (either Firefox or Chrome). The picture in the READMEs is a representation in Firefox.

11.	Start by choosing the ‘Create ENCORE File’ section to calculate your [ECHO](https://github.com/delosh653/ECHO) output for ENCORE and download resulting file. ENCORE downloads appear in the ‘downloads’ folder of the ENCORE Shiny App folder that you’ve downloaded.

12.	Once you’ve downloaded the ENCORE file, upload it to the ‘Explore’ section and proceed, following the navigation instructions.

13.	Have fun!

** Why do I have to install either Firefox or Chrome, you ask? Why not Internet Explorer, or some other browser? Well, it is known there are problems downloading files when viewing shiny apps in Internet Explorer, as well as some aspects of d3 not working, so we definitely want to avoid that. However, I have not tested this app in browsers like Microsoft Edge, Safari, etc. If you can verify that these work, please let me know at delosh@rpi.edu.

## ENCORE Features

ENCORE's interface is divided into two sections: **Explore** and **Create Encore File**. 

Within the **Create ENCORE File** tab, you can upload your ECHO visualization file results (.csv) and enter ontology enrichment information, such as the organism, gene ID type, desired significance level and adjustment, and period range restrictions. Note that this step takes some time. You can then download your results as an .RData file, which appears in the 'downloads' folder of the ENCORE Shiny App folder.

<p align="center">
<img src="ENCORE Shiny App/www/ont_map.PNG" width="400" />

<p align="center">
<img src="ENCORE Shiny App/www/hover_ont_map_node.png" width="200px">  <img src="ENCORE Shiny App/www/hover_ont_map_link.png" width="200px"> 
</p>

In the **Explore** tab, simply upload the ENCORE .RData file from your results and select desired ontologies and categories. You can then see a sankey, or flow, diagram showing the path to a selected ontological term for a specific AC category. Possible paths from the ontological type to the selected term flow from left to right. Not significant terms are grey and significant terms are colored. To jump to a specific ontological category's children, click its bar node. Hover and click on bars and links for more information.

<p align="center">
<img src="ENCORE Shiny App/www/ont_nav.PNG" width="400" />

<p align="center">
<img src="ENCORE Shiny App/www/hover_fc_bar.png" width="200px">  <img src="ENCORE Shiny App/www/hover_ont_bar.png" width="200px"> 
</p>

You can then explore a layered bar graph in the Ontology Explorer, representing the ontologies significantly enriched for the categories you chose. You can then hover and click for more information on enirchments and ontological terms. 

<p align="center">
<img src="ENCORE Shiny App/www/group_comp.PNG" width="400" />

<p align="center">
<img src="ENCORE Shiny App/www/hover_fc_group.png" width="200px"> <img src="ENCORE Shiny App/www/hover_gene_group.png" width="200px"> <img src="ENCORE Shiny App/www/hover_gene_connect.png" width="200px"> <img src="ENCORE Shiny App/www/hover_hs.png" width="200px">
</p>

Once you've selected a GO Term in the Ontology Explorer, you can explore protein-protein ineractions between different genes belonging to that category, as well as their mean-centered, normalized heat maps of gene expression sorted by phase. We can also find more information about these connections through different hover and click interactions.

<p align="center">
<img src="ENCORE Shiny App/www/gene_expr.PNG" width="400" />
</p>

When you click on an ontological term in the Ontology Explorer or a gene in the Group Comparison, information about the term appears here. If you click on a gene from the Group Comparison Tool, the expression with summary information will appear on the left, with a link connecting to gene information on above. Genes appearing in the chord diagram also appear here, along with their amplitude change coefficient categories and hours shifted (phase) values.

Note: further descriptions and instructions appear within the app.

## Contact Information and Bug Reporting

As you may have noticed, this is still in beta testing! Therefore, we would love to hear your feedback on the program. For general feedback, email delosh@rpi.edu with the subject line "ENCORE Feedback".

If you run into any errors, please email delosh@rpi.edu with the following (subject line: "ENCORE Error"): 
- a short desciption of your problem
- ENCORE version number 
- your ENCORE file 
  - Feel free to compress, as it is a large file. If it's too large after compression, simply send the ECHO .RData file and your exact inputs to get the ENCORE file.
- your exact settings and click path (this appears in the top right of any image)
- your exact error from the console window (a screenshot will do)

However, please read the FAQ below before sending error reports.

Contact:
Hannah De los Santos /
email: delosh@rpi.edu /
Rensselaer Polytechnic Institute

## FAQ

**Q:** My dataset isn't working for some reason and I'm not using the latest ECHO version! Why?

**A:** Please update to the current ECHO version (you can see this at the top of the github repository by checking out what the latest commit was), since this may have been corrected in an update. If this problem persists, feel free to send another email!
