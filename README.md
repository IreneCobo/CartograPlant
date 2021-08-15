# Imputation-based meta-analysis workflow for GEA and GWAS analyses in Galaxy (Tripal Galaxy) for CartograTree

## Table of Contents
1. [Introduction](#intro)
2. [Material and methods](#mat)
    - [Workflow Structure](#next)
    - [Study panels and reference panel](#panels)
    - [Phenotypic variables](#phen)
    - [Environmental variables](#env)
    - [Data quality filtering](#dat)


<a name="intro"></a>
## Introduction  
This is the reference documentation for the imputation-based meta-analysis workflow development for GEA and GWAS analyses in Galaxy for CartograTree. 

Over the last decades, there has been an increase in high resolution data due to the emergence of new technologies: high throughput sequencing for genotypic data, phenomics for phenotypic data, and fine resolution environmental data. These technologies allow us to address adaptive questions impossible to answer to date. Nevertheless, they require intelligent workflows and advanced data integration. 

CartograTree is the first web-based application which integrates genotypic, phenotypic and environmental data, and associated meta-data, from georeferenced plants, connected by analytic workflows implemented in Galaxy. The development of this robust and integrative resource is relevant and timely since consistent data collection for plant population studies is lacking. The functionalities of CartograTree include data import, data storage, data selection/filtering, data analysis and data integration. All organisms represented in CartograTree must be georeferenced to provide integration with environmental layers, one of the main powers of the platform. Another advantage is derived from the ability to analyze the stored data directly through CartograTree thanks to the development of Tripal Galaxy, which helps end-users easily integrate data from the site within workflows. 

This imputation-based meta-analysis workflow will improve the ability to analyze multiple datasets from independent studies in integrative platforms such as CartograTree. On the one hand, meta-analyses will allow the integration of genomic, phenotypic and environmental data from different studies available in CartograTree. On the other hand, imputation methods are useful when working with low-density SNP assays. They work by combining a reference panel of individuals genotyped at a dense set of polymorphic sites (usually single-nucleotide polymorphisms, or “SNPs”) with a study sample collected from a genetically similar population and genotyped at a subset of these sites [(Howie *et al*. 2009)](https://doi.org/10.1371/journal.pgen.1000529). Hence, they predict unobserved genotypes in the study sample by using a population genetic model to extrapolate allelic correlations measured in the reference panel [(Howie *et al*. 2009)](https://doi.org/10.1371/journal.pgen.1000529). The imputed genotypes expand the set of SNPs that can be tested for association, and this more comprehensive view of the genetic variation in a study can enhance true association signals and facilitate meta-analysis [(Zegginni *et al*. 2008](https://doi.org/10.1038/ng.120), [Barret *et al*. 2008)](https://doi.org/10.1038/ng.175).

*Populus trichocarpa* is used as model species to develop this workflow given its compact genome size (~500 Mb, 19 chromosomes) and the high number of genomic and phenotypic resources available, since it was the first sequenced tree genome [(Tuskan *et al*. 2006)](https://science.sciencemag.org/content/313/5793/1596). It is a deciduous broadleaf tree that is native to Western North America. It is an economically important source of timber. Its rapid growth, together with its compact genome size, has lead to its use as a model organism for tree species. 

<a name="mat"></a>
## Material and methods

The study panels, reference panel, as well as their phenotypic and environmental information, are stored in: 

/labs/Wegrzyn/IreneCobo/Ptrichocarpa

<a name="next"></a>
### Workflow Structure

Our meta-GWAS/GEA approach was based on [Zhao et al. 2019](https://doi.org/10.1038/s41467-019-09462-w) imputation-based meta-GWAS workflow for tomato, with some modifications for meta-GEA analysis.


![Workflow](galaxy_workflows/Imputation-based.meta-analysis/Workflow-Imputationbasedmetaanalysis.png "Workflow")


<a name="panels"></a>
### Study panels and reference panel

We used three different GWAS panels already published and genotyped using different technologies (Table 1). 

We imputed SNP data for panels 1 and 3 from a reference panel (Table 1). Panel 3 was genotyped by using a low-density assay (34K). Panel 1, although it was sequenced using a whole genome techology (exome capture), given the high number of missing genotypes in chromosome 1 and 19, these chromosomes were imputed as well.

Raw reads from panel 1 and reference panel were downloaded from SRA in fastq format. For panel 2, high quality SNPs (table 1), ready to use, were downloaded in .vcf format. For panel 3, a .csv file containing the SNPs from a 34K assay, was downloaded. The references for the assay design can be found in table 1. 


**Table 1**. Study panels and reference panel used for the workflow development.

| **Variables** | **Panel 1** | **Panel 2** | **Panel 3** | **Reference panel** |
| ------ | ------ | ------ | ------ | ------ |
| *Number of individuals* | 411  (451 in article) | 544 | 448 | 461 clones, 101 provenances |
| *Number of loci (SNPs)* | Raw reads (1,311,373 in the article) | 17,902,740 | 29,355 | Raw reads (813,280 in the article) |
| *DNA sequencing technique* | Exome capture (design of oligonucleotide baits based on P.trichocarpa genome v2.0) | High-coverage whole-genome sequencing (WGS) | 34 K Populus SNP array (Assay design information in: [Geraldes et al. 2011](https://doi.org/10.1111/j.1755-0998.2010.02960.x) and  [2013](https://doi.org/10.1111/1755-0998.12056)| Exome capture (design of oligonucleotide baits based on *P. trichocarpa* genome v2.0) |
| *Reference genome* | v.2.0 | v.3.0 | v.2.2 | v.2.0 |
| *Latitude range (N)* | Latitudinal transect: 37.537433 to 61.01; Altitudinal transect: 49.27878 to 50.38812 | 38.882-54.25 | 44-59.62 | 48°54′ N  (Nooksack River, Whatcom County, Washington) - 43°47′ N (Middle Fork, Willamette River, Lane County, Oregon) |
| *Longitude range (W)* | Latitudinal transect: -118.70285 to -153.582; Altitudinal transect: -121.00119 to -123.67191 | 120.087-128.683 | 121.17-137.92 | |
| *Elevation range (m)* | Latitudinal transect: 0-2363; Altitudinal transect: 7-1161 | 3.048-2144 | 0-900 | |
| *Phenotypes* | Two common garden locations: Crizt (Virginia): Bud flush, bud set and height; British Columbia (Canada): Bud flush, bud set, height (among others in the article) | Bud flush, bud set and height | Bud flush, bud set and height (among others in the article) |
| *Reference article* | [Zhang et al. 2019](https://doi.org/10.1093/gbe/evz151), [Oubida et al. 2015](https://doi.org/10.3389/fpls.2015.00181) | [Evans et al. 2014](https://doi.org/10.1038/ng.3075) | [McKown et al. 2014a](https://doi.org/10.1111/nph.12815), [b](https://doi.org/10.1111/nph.12601) | [Guerra et al. 2019](https://doi.org/10.1186/s12864-019-6160-9) |

The individual locations for the individuals of each study panel, as well as their sequencing technology, are represented in the following figure.

![Panellocations](galaxy_workflows/Imputation-based.meta-analysis/Panellocations+sequencing.png "Panellocations")

<a name="phen"></a>
### Phenotypic variables
Regarding phenotypic data, studies containing bud set, bud flush and height phenotypic data for *P. trichocarpa* were selected. Phenology traits, such as bud set and bud flush are highly related to fitness and show high heritability. On the other hand, height in trees and growth-related phenotypes in general can be considered as a proxy for fitness and they are also traits associated with the genotype. Thus, a strong selection over all these phenotypic traits is expected [(Aitken and Bemmels 2015)](https://doi.org/10.1111/eva.12293).

The phenotypic information for each study panel is represented in table 2

**Table 2**. Phenotypic information for each study panel

| **Variables** | **Panel 1 replicated only in space** |  | **Panel 2 replicated in space (3 replicates per phenotype) and in time (measured in different years)** |  |  | **Panel 3 (replicated in space (4–20 clonal ramets of similar age and condition) and in time (repeated measurements across years))** |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| *Common garden locations* | British Columbia CA | Crizt VA | Clatskanie OR | Corvallis OR | Placerville CA | Totem field common garden (UBC) |
| *Bud flush (Budflush day CO_357:1000007, budflush date CO_357:0000226, Poplar bud flush scoring scale CO_357:3000010)* | Number of days to the first fully unfolded leaf | Number of days to the first fully unfolded leaf | Bud flush score (0 buds tightly closed - 5 leaves fully expanded), on 29 March, 2012, 26 March 2013, 1 April 2013 | Bud flush score (0 buds tightly closed - 5 leaves fully expanded), on 29 March, 2012, 26 March 2013, 1 April 2013 | Bud flush score (0 buds tightly closed - 5 leaves fully expanded), on 29 March, 2012, 26 March 2013, 1 April 2013 | Julian date (2010-2011) Phenology events were marked using visual observations of the terminal bud on the main bole or canopy as a whole |
| *Bud set (Budset date CO_357:1000009, CO_357:0000016, budset scoring CO_357:1000010, CO_357:0000017, Poplar budset scoring scale CO_357:3000016)* | Days elapsed for trees to reach a fully developed bud from January 1 | Days elapsed for trees to reach a fully developed bud from January 1 | Bud set score (1 actively growing - 6 large buds fully formed), on Fall 2009, 26-31 August, 2010, 15-17 September, 2010 and 28-30 September 2010 | Bud set score (1 actively growing - 6 large buds fully formed), on Fall 2009, 26-31 August, 2010, 15-17 September, 2010 and 28-30 September 2010 | Bud set score (1 actively growing - 6 large buds fully formed), on Fall 2009, 26-31 August, 2010, 15-17 September, 2010 and 28-30 September 2010 | Julian date (2008-2010/2009-2010: Bud set dates occurring before the summer solstice (day 186) removed) Phenology events were marked using visual observations of the terminal bud on the main bole or canopy as a whole |
| *Height (Tree total height CO_357:1000037, CO_357:0000048)* | cm (after all trees set bud in the fall of 2013) | cm (after all trees set bud in the fall of 2013) | cm (2009-2012) | cm (2009-2012) | cm (2009-2012) | cm (2008-2012, estimated at the end of each season) |

<a name="env"></a>
### Environmental variables

Regarding environmental data, we used temperature and photoperiod, since they are the environmental variables associated with the fitness-related phenotypic data analyzed in our study. However, we also took into consideration latitude and altitude. Since temperature and photoperiod co-vary with latitude, it is difficult to differentiate the effect of temperature and photoperiod due to latitude. Thus, we subsetted samples for the same latitude located through altitudinal gradients to perform our GEA analyses. Since they are located at the same latitude, they are also subjected to the same photoperiod. Thus, this design allows us to test the effect of temperature alone [(Aitken and Bemmels 2015)](https://doi.org/10.1111/eva.12293).

<a name="dat"></a>
### Data quality filtering

Raw sequences for panel 1 and reference panel were downloaded in .fastq format from SRA (NCBI) (accessions SRP055524 and SRA058855, respectively) running the following command in SRA Toolkit:

```
for (( i = 206; i <=616; i++))
do 
./fastq-dump --split-files SRR1819$i
done

```
(This is the code for panel 1. Same for the reference panel, changing the accessions)

Then, they had to be filtered by data quality and mapped against reference genomes. Afterwards, a SNP calling was performed.

First, sequences were trimmed for adapters and bases for quality <35 (based on the Illumina filter flag). In addition, very short reads (length <35) were also removed.

To accomplish this, FASTQC v0.11.7 was used to detect the adapters.

Afterwards, trimmomatic was used to remove the adapters and filter bases with quality <35, as well as very short reads (length<35). Read having poor local quality scores (score<25 and window size =5) were also removed using trimmomatic.



