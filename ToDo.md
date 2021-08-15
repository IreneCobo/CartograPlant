### To Do's

**21 July 2021**
- [ ] Meta-analysis with METASOFT: Reading the manual
- [ ] Landscape genomics with LFMM 2: Writing the R code to be wrapped in Galaxy
- [ ] Remapping SNPs against another reference genome
- [ ] CartograPlant article: Including last updates

**14 June 2021**

Fraxinus: https://gitlab.com/PlantGenomicsLab/greenash_radseq/-/blob/master/EMMAX/EMMAXfixed.md 
- [x] EMMAX performed in the imputed set of SNPs using LinkImputeR (only correcting for kinship, no for population structure): 35941 SNPs for HKP phenotype and 33741 SNPs for MLW phenotype were statistically significant (P<0.05) before the multiple testing correction. But, after multiple testing correction, I obtained again 0 statistically significant SNPs (P<0.05) for all the methods (FDR, Bonferroni and qvalue) and for both phenotypes (HKP and MLW).
- [ ] Finding a way to perform LD calculations (LD prune) using the genetic map

Brainstorm with Gabe:
- [ ] Finite-sample genome-wide regression p-values (GWRPV) with a non-normally distributed phenotype: https://doi.org/10.1101/204727, https://cran.r-project.org/web/packages/gwrpvr/gwrpvr.pdf
- [ ] Bayesian approach (instead of frequentist):https://books.google.com/books?id=A18ICgAAQBAJ&pg=PA243&lpg=PA243&dq=lddesign+r+package&source=bl&ots=AGfyl7bvnK&sig=ACfU3U1Z65Gq-4VqpgDT9MRPL_vXFoHbew&hl=es&sa=X&ved=2ahUKEwitvZeLwJbxAhXIWM0KHakNAoUQ6AEwCXoECAUQAw#v=onepage&q=lddesign%20r%20package&f=false (page 242), https://doi.org/10.1007/978-1-62703-447-0_3, ldDesign in R https://pubmed.ncbi.nlm.nih.gov/23756887/, BayPass http://www1.montpellier.inra.fr/CBGP/software/baypass/files/BayPass_manual_2.1.pdf . Problem: Correction by kinship/population structure

**24 May 2021**

Fraxinus: https://gitlab.com/PlantGenomicsLab/greenash_radseq/-/blob/master/EMMAX/EMMAX.md
- [x] EMMAX performed in the imputed set of SNPs using LinkImputeR. I obtained 35188 and 35584 statistically significant SNPs (P<0.05) for HKP and MLW, respectively. After multiple testing correction (FDR and Bonferroni), in both cases and phenotypes, no significant SNPs.
- [x] Since I obtain significant SNPs before the multiple testing correction, I have been reading bibliography and apparently, since Linear Mixed Models (such as EMMAX) correct by population structure and kinship for false positives, FDR and Bonferroni can be too restrictive. I tried using MultiTrans, a multitesting correction program for LMM. However, I obtained an error when trying to run it: Error: cannot allocate vector of size 6453.2 Gb
- [x] I have also tried to perform QTL for half-sib families. I found the GridQTL program, which perform this type of analysis (https://www.researchgate.net/publication/266499255_gridqtl_A_Grid_Portal_for_QTL_Mapping_of_Compute_Intensive_Datasets). There are some articles published that used this program in half-sib families, but I do not manage to find the program online, nor the users manual. I have written the GridQTL developers some days ago but they have not answered yet. I am thinking about writing the authors of the several published articles using this program in both full-sib and half-sib families, in case they have the software and the manual. 

**10 May 2021**
- [x] 2021 Forest genetics poster and poster talk submitted
- [x] GCC2021 abstract submitted
- [x] ISG registration, poster and poster talk submitted 

- Workflows Galaxy: 
  - [ ] Remapping SNPs: Working on the last step: Include the new SNP positions in the vcf file (so far so good)
- Fraxinus: 
  - I need clarification for some steps
 
**26 April 2021**
- Poster Forest Genetics 2021: In operation
- Article CartograPlant: In operation 
- Workflows Galaxy: 
  - Alternatives to DAPC for accounting for Population Structure: PCA seems to be a good option (although not perfect): https://www.broadinstitute.org/talks/controlling-stratification-meta-gwas-pca-theory-applications-and-implications and https://faculty.washington.edu/tathornt/SISG2015/lectures/assoc2015session06.pdf
  - Quality filtering: done
  - Remapping SNPs: Filtering the sam file obtained after mapping the flanking regions against the most recent reference genome with bwa (because there are some flanking regions that map agains multiple sites): https://gist.github.com/davfre/8596159 and https://wikis.utexas.edu/display/CoreNGSTools/Filtering+with+SAMTools#FilteringwithSAMTools-Filteringhigh-qualityreads 

 1. To keep only reads that map without any mismatches (or it is too much restrictive?):

 ```
bamtools filter -tag XM:0 -in reads.bam -out reads.noMismatch.bam
```
2. Keep only the duplicate with the highest quality. If two or more duplicates have the same quality, remove both.

Remove unmapped (0x4) and secondary reads (0x100), and retain the highest quality alignments (q 20)
```
samtools view -S -F 0x104 -q 20 flankchr1.sam #201837 out of 210878 original SNPs
```
I tried with q 10
```
samtools view -S -F 0x104 -q 20 flankchr1.sam #204284 out of 210878 original SNPs
```


2. 

**22 March 2021**
- Fraxinus: https://gitlab.com/PlantGenomicsLab/greenash_radseq/-/blob/master/EMMAX.md
Success! I tried performing a kinship matrix using the LD0.2 SNPs and the association analysis run without erros. Now I am trying to figure out a threshold for the P-value (correct for multiple testing): FDR, Bonferroni, but I have read about P = 1/n, being "n" the effective number of independent SNPs. Reading bibliography to try to understand this last method. 
- DAPC analysis in Panel2. vcf file too big (WGS), need to run R in an interactive session using screen. I asked for help in the CBC slack and I followed their guidelines but still having the same issue (my screen session is aborted when I close my terminal), plus a new issue (not only my screen session is aborted when I close my terminal, but also, now, my srun interactive session is also aborted). I talked again with them (via Zoom call) and they said that, by default, the screen session stop running authomatically after some hours. They told me that they will let me know if they manage to find a solution. Still waiting


**5 March 2021**
- DAPC analysis in Panel2. vcf file too big (WGS), need to run R in an interactive session using screen. I asked for help in the CBC slack and I followed their guidelines but still having the same issue (my screen session is aborted when I close my terminal), plus a new issue (not only my screen session is aborted when I close my terminal, but also, now, my srun interactive session is also aborted)
- Galaxy workflow: Remapping SNPs against new reference genome almost finished. Three issues:
    - Genome indexing is not working (neither with bwa nor samtools)
    - Doubt with the sam file quality filtering,resulting from the mapping of the SNP flanking regions against the new reference genome (is this code correct? I filter by mapping quality (20) but it is specific of the aligner/it varies with the aligner, would it be a flag-based filtering more accurate?: samtools view -S -F 0x104 -q 20 SNP.sam )
    - Figuring out how to include the resulting SNP positions in the vcf file (almost there, but I am having problems in "reheading" the new vcf file with the v4 positions)
- Reference panel (Genetic map code): Pending
- Fraxinus: Trying to use EMMAX. Error. I have been investigating but I do not figure out the cause
```
./emmax-intel64 -v -d 10 -t genotype -p phenoEMMAXorderHKP.txt -c assignedPoporder.txt -k output.aBN.kinf -o HKPassoc

Reading TFAM file genotype.tfam ....


Reading kinship file output.aBN.kinf...

  98 rows and 98 columns were observed with 0 missing values.


Reading the phenotype file phenoEMMAXorderHKP.txt...

  98 rows and 1 columns were observed with 7 missing values

nmiss = 7 , mphenoflag = 1

Reading covariate file assignedPoporder.txt...

  98 rows and 2 columns were observed with 0 missing values. Make sure that the intercept is included in the covariates

File reading - elapsed CPU time is 0.000000

evals[0] = nan, evals[1] = nan, evals[n-1] = nan
FATAL ERROR : Minimum q eigenvalues of SKS is supposed to be close to zero, but actually 179769313486231570814527423731704356798070567525844996598917476803157260780028538760589558632766878171540458953514382464234321326889464182768467546703537516986049910576551282076245490090389328944075868508455133942304583236903222948165808559332123348274797826204144723168738177180919299881250404026184124858368.000000
FATAL ERROR : Minimum q eigenvalues of SKS is supposed to be close to zero, but actually 179769313486231570814527423731704356798070567525844996598917476803157260780028538760589558632766878171540458953514382464234321326889464182768467546703537516986049910576551282076245490090389328944075868508455133942304583236903222948165808559332123348274797826204144723168738177180919299881250404026184124858368.000000
FATAL ERROR : k = 90 > n-q = 89
Aborted
```

**22 February 2021**
- **Reference panel**: Still working on the reference panel's code, no relevant updates

**8 February 2021**
- **Reference panel**: Apparent still running in himem partition (500 GB memory, 12 cores). I am running two files, the chr3 for all individuals (.out looks fine, 21 % run so far, no errors in the .err), and all the individuals/all chromosomes (no errors in the .err file but anything in the .out, as usual). I was wondering, however, if it makes sense to run apparent in the Reference panel, since it is no family-based data. In the reference article I am using for the workflow (https://www.nature.com/articles/s41467-019-09462-w#Sec8), the authors mention that the developed a Python code to perform the genetic map necessary for imputation (they do not mention any pedigree information). I have been looking for this code in the article/supplemental information but I found nothing. I wrote the authors asking them for the code and they sent me the code this morning. I am reading it and trying to understand it.

**25 January 2021**
- **Reference panel**: Tried running Apparent in 10 individuals to test if the memory error was in fact a input format error: It run perfectly (thanks to Jeremy advice of loading module load R/3.5.1 in my .sh code). Two outputs: "Triad_allchr3_10ind.csv" and "Triad_sigchr3_10ind.csv", three significant tests showing the most probable parents for the three significant individuals (p-value=0). Now running all the individuals on chr3 using the "himem" partition (Apparentcode2moremem.sh). It has been running for 2 days and a half approx. (still running) Working directory: /labs/Wegrzyn/IreneCobo/Ptrichocarpa/Panel2/Referencepanel/Jeremytest/Testchr01  

```
> library(outliers)
> InputFile <- read.table(file="chr310ind.txt",sep="\t",h=F)
> apparentOUT <- apparent(InputFile, MaxIdent=0.10, alpha=0.01, nloci=300, self=FALSE, plot=TRUE, Dyad=FALSE)
  |======================================================================| 100%
There were 20 warnings (use warnings() to see them)
> warnings()
Warning messages:
1: In min(genoOut$GD) : no non-missing arguments to min; returning Inf
2: In max(genoOut$GD) : no non-missing arguments to max; returning -Inf
3: In min(genoOut$GD) : no non-missing arguments to min; returning Inf
4: In max(genoOut$GD) : no non-missing arguments to max; returning -Inf
5: In min(genoOut$GD) : no non-missing arguments to min; returning Inf
6: In max(genoOut$GD) : no non-missing arguments to max; returning -Inf
7: In min(genoOut$GD) : no non-missing arguments to min; returning Inf
8: In max(genoOut$GD) : no non-missing arguments to max; returning -Inf
9: In min(genoOut$GD) : no non-missing arguments to min; returning Inf
10: In max(genoOut$GD) : no non-missing arguments to max; returning -Inf
11: In min(genoOut$GD) : no non-missing arguments to min; returning Inf
12: In max(genoOut$GD) : no non-missing arguments to max; returning -Inf
13: In min(genoOut$GD) : no non-missing arguments to min; returning Inf
14: In max(genoOut$GD) : no non-missing arguments to max; returning -Inf
15: In min(genoOut$GD) : no non-missing arguments to min; returning Inf
16: In max(genoOut$GD) : no non-missing arguments to max; returning -Inf
17: In min(genoOut$GD) : no non-missing arguments to min; returning Inf
18: In max(genoOut$GD) : no non-missing arguments to max; returning -Inf
19: In min(genoOut$GD) : no non-missing arguments to min; returning Inf
20: In max(genoOut$GD) : no non-missing arguments to max; returning -Inf
> write.csv(apparentOUT$Triad_all,"Triad_allchr3_10ind.csv")
> write.csv(apparentOUT$Triad_sig,"Triad_sigchr3_10ind.csv")
```
- **EAB (green ash)**: I checked and corrected the families and individuals for v5 (I included some extra families in error), but in spite of that, Mendelian errors continued. Jeremy developed a mendelian filter script to fix the Mendelian error and performed new .ped, .dat and .map files. We both run Merlin independently using these corrected input files (QTL and association analysis) and we got statistically significant results. In addition, I performed a family informativeness analysis and also, I redid the QTL analysis and family informativeness analysis after removing unlikely genotypes using --error flag in Merlin. Detailed information about all the analysis can be found here: https://gitlab.com/PlantGenomicsLab/greenash_radseq/-/blob/master/Merlin_analysis/MerlinAnalysisv1.4.md . An Excel file with the summary of the results of all the analyses can be found here: https://drive.google.com/file/d/1X7o9p2XAq2mCPUwrzcleNHiYPsypn7Em/view?usp=sharing 

All the information about the analyses and results can be found in Gitlab: https://gitlab.com/PlantGenomicsLab/greenash_radseq/-/blob/master/Merlin_analysis/MerlinAnalysisv1.4.md 

In summary:
QTL results:
- Host kill percentage: 121 significant SNPs (P<0.05) out of 3650
- Mean_larval_height: 135 significant SNPs (P<0.05) out of 3651

After removing the unlikely genotypes:
- HKP_nounlikelygenotypes: 111 significant SNPs (P<0.05) out of 3654
- MLH_nounlikelygenotypes: 205 significant SNPs (P<0.05) out of 3645

Family informativeness results were the same in both cases:
======================
         Family           Trait  People  Phenos   Pairs    Info  ELOD20
           C__1         Model 1      11      10      45   4.535   0.039
           C__1         Model 2      11      10      45   1.958   0.017
           HA_1         Model 1       7       7      21   4.453   0.039
           HA_1         Model 2       7       7      21   0.473   0.004
           I__1 *** 36-bit family skipped ***
           J__1         Model 1       9       8      28   1.515   0.013
           J__1         Model 2       9       8      28   1.658   0.014
           K__1 *** 32-bit family skipped ***
          TOTAL         Model 1      27      25      94  10.502   0.091
          TOTAL         Model 2      27      25      94   4.089   0.036


**Doubt**: What concerns me the most right now is the possibility of having a lot of false positives due to the high number of p-values obtained (resulted from multiple testing at multiple loci). I was wondering if it would be advisable to perform a False Discovery Rate (FDR) or Bonferroni correction in the p-values. The QTL analysis we used in Merlin is based on this method http://csg.sph.umich.edu/abecasis/publications/pdf/Am.J.Hum.Genet.vol.71-pp.238.pdf. The authors tested the type I errors (fase positives) of the method and they concluded that "For normally distributed traits and **in large samples**, the method is found to give the correct type I error rate and an unbiased estimate of the proportion of trait variance accounted for by the additive effects of the locus—although, **in cases where asymptotic theory is doubtful, significance levels should be checked by simulations**.". Our sample is not large at all, so based on this, it should be advisable for us to check the significance levels by simulations. I am not very familiar with this way to check significance levels (using simulations) but seems like the --simulate flag in Merlin could do the trick? http://csg.sph.umich.edu/abecasis/merlin/tour/simulation.html.

**6 January 2021**
- Reference panel: Error on December 30th (after 7 days running): There is not package called "outliers" (in himem I suppose). I have tried to installed the package including "install.packages("outliers")" in the .R document, but it gives me the following error: 

```
Installing package into ‘/usr/lib64/R/library’
(as ‘lib’ is unspecified)
Warning in install.packages("outliers") :
  'lib = "/usr/lib64/R/library"' is not writable
Error in install.packages("outliers") : unable to install packages
Execution halted
```
I have tried running the same code (400GB Memory and 20 CPUs) but in "general": /labs/Wegrzyn/IreneCobo/Ptrichocarpa/Panel2/Referencepanel/Jeremytest/Testchr01/Apparentcode2morememgeneral.sh. Another error:

```
sbatch: error: Memory specification can not be satisfied
sbatch: error: Batch job submission failed: Requested node configuration is not available
```
I have been investigating how to install R packages in Xanadu, but I have just found a tutorial about how to install them locally (home directory). Not sure how to install it in himem partition

**23 December 2020**
- Rerun Apparent using 400GB Memory and 20 CPUs in "himem" instead of "general". Working directory: /labs/Wegrzyn/IreneCobo/Ptrichocarpa/Panel2/Referencepanel/Jeremytest/Testchr01. Job name: Apparentcode2moremem.sh Code:
```
#!/bin/bash
#SBATCH --job-name=apparent
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mail-type=END
#SBATCH --mem=400G
#SBATCH --mail-user=irene.cobo_simon@uconn.edu
#SBATCH -o apparentmoremem_%j.out
#SBATCH -e apparentmoremem_%j.err

Rscript Apparentcode2.R
```

**9 December 2020**

- Poster for ISG's event
- Fraxinus: 
    - Rerun the analysis using the v3 of the genotype file and the new map file (removed families whose parents are not sure, and filtered by Mendelian errors/Quality filtered): https://gitlab.com/PlantGenomicsLab/greenash_radseq/-/tree/master/Filtering/v3: Adapting the input file formats to Merlin
        - .dat file: done
        - .map file: done
        - .ped file: in operation
- Reference panel: Rerun apparent using one chromosome to check whether it is a memory issue or not: in operation

**1st December 2020**
 
- Fraxinus: 
    - Rerun the analysis using the v3 of the genotype file and the new map file (removed families whose parents are not sure, and filtered by Mendelian errors/Quality filtered): https://gitlab.com/PlantGenomicsLab/greenash_radseq/-/tree/master/Filtering/v3
    - Performed a GWAS analysis in the Quality filtered data
- CartograTree/Plant workflows: 
    - SNP data quality filtering workflow in the Galaxy instance: Problems with vcftools, I cannot find the flags needed for the filtering
    - Population structure analysis workflow: Performed in Xanadu, pending to be performed in the Galaxy instance
    - Imputation of SNPs: I need a pedigree file to peform a High density recombination map in [LepMap3 program](https://avikarn.com/2019-04-17-Genetic-Mapping-in-Lep-MAP3/) to be used in IMPUTEv2: To do so, I am using Apparent program. Memory issues. Looking for a more efficient program to perform this analysis (pedigree file) in larger samples sizes. 
    - Use the --IMPUTE output flag in vcftools to transform the Reference panel in vcf format into IMPUTE format reference panel (http://vcftools.sourceforge.net/man_latest.html , OUTPUT OTHER FORMATS selection)

**25 November 2020**

- Reference Panel, Apparent: Redo the analysis using the apparentTest data in a separate directory (/labs/Wegrzyn/IreneCobo/Ptrichocarpa/Panel2/Referencepanel/Jeremytest/ApparentTest). 

1. Tried using Dyad=TRUE, self=TRUE it gave me an error (Error in order(Out3a$DOf, Out3a$DCumPv) : argument 1 is not a vector
Calls: apparent -> [ -> [.data.frame -> order
Execution halted). 
2. Tried using Dyad=FALSE, self=TRUE, any errors (Triad_all3.txt, Triad_sig3.txt). 
3. Tried using Dyad=FALSE, self=FALSE, any errors neither (Triad_all4.txt, Triad_sig4.txt). 

Running Apparent for Reference Panel in a separate directory (/labs/Wegrzyn/IreneCobo/Ptrichocarpa/Panel2/Referencepanel/Jeremytest/ApparentReferencepanel) and with Dyad=FALSE instead of TRUE, self=FALSE, 200GB memory again

I canceled the job because .out and .err files did not show any errors but it has been running for more than a week. Investigating where can be the problem. 

**11 November 2020**

- [ ] Panel 1: ValidateVariants to check the validity of the g.vcf files (same error in different positions): 

- 206 sample: A USER ERROR has occurred: Input 206.sort.mkdup.g.vcf fails strict validation of type ALL: one or more of the ALT allele(s) for the record at position Chr01:11562 are not observed at all in the sample genotypes
- 207 sample: A USER ERROR has occurred: Input 207.sort.mkdup.g.vcf fails strict validation of type ALL: one or more of the ALT allele(s) for the record at position Chr01:9983 are not observed at all in the sample genotypes
- 208 sample: A USER ERROR has occurred: Input 208.sort.mkdup.g.vcf fails strict validation of type ALL: one or more of the ALT allele(s) for the record at position Chr01:8580 are not observed at all in the sample genotypes

I looked the error on the GATK page and found [this answer](https://gatk.broadinstitute.org/hc/en-us/community/posts/360061452132-GATK4-RNAseq-short-variant-discovery-SNPs-Indels-), so I am following the GATK team guidelines: Redoing the SNP calling using the newest version of GATK (GATK/4.1.8.1), since I used the previous version (GATK/4.1.3.0). Aftewards, I will check whether the error persist.

I run it in GATK/4.1.8.1 and same error: A USER ERROR has occurred: Input 206.sort.mkdupnew.g.vcf fails strict validation of type ALL: one or more of the ALT allele(s) for the record at position Chr01:11562 are not observed at all in the sample genotypes

- [ ] Reference panel *P. trichocarpa*: Running Apparent in the Reference panel to obtain pedigree data, needed as an input for the Genetic map building(/labs/Wegrzyn/IreneCobo/Ptrichocarpa/Panel2/Referencepanel/Jeremytest): 
1. Jeremy helped me to transform the Reference Panel vcf file (ReferencepanelSNPposv4fixed2sorted.vcf) into Apparent format (converted2.txt) using his python code (convert_vcf_to_format.py). I run it in Xanadu using the following job (Convertvcftoapparent.sh). Output (Reference panel in Apparent format): converted2.txt
2. Running Apparent. To do so, it is necessary, first of all, to run the [Apparent R Script](https://github.com/halelab/apparent/blob/master/apparent.R) and afterwards, the Apparent code

```
library(outliers)
InputFile <- read.table(file="converted2.txt",sep="\t",h=F)
apparentOUT <- apparent(InputFile, MaxIdent=0.10, alpha=0.01, nloci=300, self=FALSE, plot=TRUE, Dyad=TRUE)
write.csv(apparentOUT$Triad_all,"Triad_allRefPan.csv")
write.csv(apparentOUT$Triad_sig,"Triad_sigRefPan.csv")
```
- self: Logical value for instructing 'apparent' whether or not to consider self-crossing (parent i = parent j). The default value is TRUE. I changed it into "FALSE" since P.trichocarpa is a dioic species with male and female flowers located in different individuals (the same for Fraxinus pennsylvanica). 
- Dyad: Logical value for instructing 'apparent' to perform a Dyad analysis, following the Triad analysis. The default value is FALSE. I changed it to "TRUE", since Dyad analysis can be employed to identify a likely single parent for a given offspring.

The rest of flags were kept with the default values

I put the Apparent R Script + the Apparent code together in the same file (first the Apparent R script and the the Apparent code): Apparentcode2.R

I run it as a R batch job: Apparentcode2.sh

```
#!/bin/bash
#SBATCH --job-name=apparent
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=100G
#SBATCH --mail-user=irene.cobo_simon@uconn.edu
#SBATCH -o apparent_%j.out
#SBATCH -e apparent_%j.err

Rscript Apparentcode2.R
```

It has been running for more than 2 days (I had to cancel it). It gave the following error
```
slurmstepd: error: get_exit_code task 0 died by signal
```
I have been investigating on the Internet and it could be a memory problem. I am redoing the analysis requesting 200G of memory instead of 100GB (#SBATCH --mem=100G). I have also ask a question in the CBC Slack group about the meaning of this error and how to fix it. Waiting for the new job results and the answer for my question


- [x] *F. pennsylvanica*: Run the more accurate association test in merlin (likelihood-ratio test, ./merlin -d pedstats.dat -p pedstats.ped -m geneticmaptab.map --assoc --pdf --tabulate) and I have also obtained a high number of significantly associated (pvalue<0.05) SNPs in both traits (5500 in total, 488 HKP (merlinassocHKP0.05.csv) and 5012 MLW (merlinassocMLW0.05.csv), /labs/Wegrzyn/IreneCobo/Fraxinus/Merlinassocsig0.05tab.txt ). PDF with figures: merlin.pdf. However, the linkage disequilibrium analysis looked suspicious (/labs/Wegrzyn/IreneCobo/Fraxinus/merlinassocout.txt), since the program did not manage to calculate the LD for all the SNPs and the LOD, H2 and pvalues were exactly the same in all the available LD results (LOD=0, H2=33%, pvalue=0.5). Thus, I have tried running a QTL analysis (./merlin-regress -d pedstats.dat -p pedstats.ped -m geneticmaptab.map -t modelfile.txt) but I have obtained the following output for all the SNPs:

```
Pedigree-Wide Regression Analysis (Model 1)
======================================================
       Position      H2   Stdev    Info     LOD  pvalue
          0.000      na      na      na      na      na
          1.299      na      na      na      na      na
          2.051      na      na      na      na      na
          7.971      na      na      na      na      na
          8.592      na      na      na      na      na
```

I have checked the input files and all look nice and run in merlin without giving any errors. It seems like the issue is not format-related but maybe related to the pedigree information provided.

I performed a genotype inference analysis (--infer) and another QTL with the infered SNPs. Same results. More detailed information about the performed analysis [here](https://gitlab.com/PlantGenomicsLab/greenash_radseq/-/blob/master/Merlin_analysis/README.md)

Trying to assign the parents to their families. Error obtained: Mendelian inheritance errors detected
FATAL ERROR - 
Mendelian inheritance errors detected

Checking the genotype alignment in the pedigree files (no errors found so far). Redoing all the input format changing into Merlin format from the geno/pheno files Jeremy sent me to make sure I did not introduce any errors in the genotype alignment during the input format changes (phenogenopedigree3.ped)

**27 October 2020**
- [ ] Run Merlin in Fraxinus pennsylvanica 
- [x] Suspicious flanking regions mapping: Tested, positions were correct
- [ ] Working on buiding genetic map using Jeremy's code and the Reference panel vcf file (successfuly obtain in tassel from the hapmap file).
    - [x] Successfully built the apparent format file from the vcf file (thanks to Jeremy's help)
    - [ ] Running Apparent to obtain the pedigree (in progress)
    - [ ] Use this [tutorial](https://avikarn.com/2019-04-17-Genetic-Mapping-in-Lep-MAP3/) to build the Genetic map (pedigree, Referencepanel in vcf)
- [ ] Perform SNP imputation of Panel 3 using the Reference panel and the Genetic map in Imputev2
- [x] Problems with hapmap file: Fixed, there was an extra column I removed and some positions had rare values (_1, -1 etc). I removed them and the code now works fine. 
Remove duplicated based on first column name (SNP names):
```
sort -t ' ' -k 1,1 -u Referencepanelposv4namescopy.txt > Referencepanelposv4namesnorep.txt
```
-t ‘ ‘ = space-delimited input file
-k 1,1 = filtering by first column
-u = only print unique values on first column

Remove the value “0” in forth column (quality filtering)
```
awk '$4 != "0"' Referencepanelposv4namesnorep.txt > Referencepanelposv4namesnorepsubset.txt
```
Remove also "_1" values in forth column 
```
awk '$4 != "_1"' Referencepanelposv4namesnorepsubset.txt > Referencepanelposv4namesnorepfixed.txt
```
Put header (I removed them to remove the first column)
```
awk 'NF' Refpanelheader.txt ReferencepanelSNPposv4noheaderfixed.hmp.txt > ReferencepanelSNPposv4fixed.hmp.txt
```

Fixed hmp.txt, Panel 4 (Panel4SNPforGWASposv4.csv), Change Chr01 by 01, AA SNPs by A, commas by tabs. Then, sort in tassel

Sort
```
./tassel-5-standalone/run_pipeline.pl -SortGenotypeFilePlugin -inputFile Panel4SNPforGWASposv4tab.hmp.txt -outputFile Panel4SNPforGWASposv4tabsorted.hmp.txt -fileType Hapmap
```
Error: there are positions in the forth column with value "-1". I removed them:
```
awk '$4 != "-1"' Panel4SNPforGWASposv4tabsorted.hmp.txt  > Panel4SNPforGWASposv4tabsortedfixed.hmp.txt
```

Finally, I run the following code in tassel to tranform hmp into vcf

Transform hmp into vcf
```
./tassel5-standalone/run_pipeline.pl -Xmx5g -fork1 -h Panel4SNPforGWASposv4tabsortedfixed.hmp.txt -export -exportType VCF -runfork1
```

It worked fine. I used the same tassel code (sort and transform into vcf) for Reference panel too and also worked fine.


**22 October 2020**

- [ ] Problems with the .bam file used to remap de SNP flanking regions for the Reference panel in Ptrichocarpa_533_v4.fa . Running the following code for finding possible errors, following this [pipeline](https://gatk.broadinstitute.org/hc/en-us/articles/360035891231-Errors-in-SAM-or-BAM-files-can-be-diagnosed-with-ValidateSamFile)
```
module load picard/2.9.2
export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch
java -jar /isg/shared/apps/picard/picard-tools-2.9.2/picard.jar ValidateSamFile I=Referencepanelflank200pbparessorted.bam MODE=SUMMARY
```
Error Type	Count

ERROR:INVALID_VERSION_NUMBER	1

ERROR:MISSING_READ_GROUP	1

WARNING:RECORD_MISSING_READ_GROUP	1698948

```
java -jar /isg/shared/apps/picard/picard-tools-2.9.2/picard.jar ValidateSamFile I=Referencepanelflank200pbparessorted.bam IGNORE_WARNINGS=true MODE=VERBOSE
```
ERROR: Header version: 1.6 does not match any of the acceptable versions: 1.0, 1.3, 1.4, 1.5

ERROR: Read groups is empty

I tried to track back possible errors. Looked at the .sam file 

```
java -jar /isg/shared/apps/picard/picard-tools-2.9.2/picard.jar ValidateSamFile I=Referencepanelflank200pbpares.sam IGNORE_WARNINGS=true MODE=VERBOSE
```
ERROR: Read groups is empty

I tried including the read groups names in the .sam file using the following code
```
java -jar /isg/shared/apps/picard/picard-tools-2.9.2/picard.jar AddOrReplaceReadGroups I=Referencepanelflank200pbpares.sam O=Referencepanelflank200paresreadgroups.sam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=1
```
I run again the code: 

```
java -jar /isg/shared/apps/picard/picard-tools-2.9.2/picard.jar ValidateSamFile I=Referencepanelflank200pbparesreadgroups.sam IGNORE_WARNINGS=true MODE=VERBOSE
```
No errors found


**21 October 2020**
- [ ] Dealing with weid things in my Reference panel file. Revising all the previous steps. Trying to use directly the hmp.txt files that Dr Guerra sent me but I just realized that the new positions after mapping against the P.trichocarpav4 do not make sense. Looking for errors in previous steps. Tried to use IGV to visualize the mapping but it does not manage to open the bam files. 
- [ ] [CartograPlant article](https://docs.google.com/document/d/14vmeAs_xEYJacCB_ocSMp-GqGEVrnme6/edit) Fully written. Working on the figures. Some doubts in comments

**14 October 2020**
- [ ] Still working on the tassel error. I tried running example files and other hmp.txt files and the program did not recognize the plugin. So I decided to download the program again. It runs now but it gave a  new error: ERROR net.maizegenetics.plugindef.AbstractPlugin - java.lang.IllegalStateException: BuilderFromHapMap: Error Parsing Line: 1	Chr01_10000376	C/T	1	10...
BuilderFromHapMap: SNP Named: 1 has illegal value: 

- I think that it is because hmp.txt file contains column names (after importing from R). Fixed. Run again. New error:

ERROR net.maizegenetics.plugindef.AbstractPlugin - java.lang.IllegalStateException: BuilderFromHapMap: Error Parsing Line: 	Chr01_10007115	C/T	1	101...
BuilderFromHapMap: SNP Named:  has illegal value:

- I compared the .hmp.txt file that run correctly (Poplar_Genotype_10.hmp.txt) with mine (Referencepanelfixed.hmp.txt). The chromosome names are different. In the good file, chr=01, in the file that gives an error: chr=1. Working on changing all the chromosome names.

- [ ] [CartograPlant article](https://docs.google.com/document/d/1GBrQOGN44o4mUAEpZZipJP52q-gX-J2B/edit) (in progress)

**7 October 2020**
- [ ] Imputation of Panel3 SNPs using Reference panel. High density recombination map: (famlily assignment to obtain the pedigree data, Jeremy's code). This information will allow me to use [Lep-Map3 program](https://avikarn.com/2019-04-17-Genetic-Mapping-in-Lep-MAP3/) to perform the recombination map. Problem: Referencepanel has to be in .vcf format and it is in hapmap (.hmp). I tried to convert it using tassel (using a code that I have worked previously in other hmp to vcf conversions). It gave me an error (Error loading: Referencepanelfixed.hmp.txt
BuilderFromHapMap: Ordering incorrect HapMap must be ordered by position). It gives this error because the .hmp.txt file has to be sorted following the tassel sorting requirements: TASSEL5 has strict requirements for the sites in a genotype file. Each site must be unique (as defined by its locus/chromosome, position, and name) and they must be in order in the file. Genotype files produced by other programs (and also earlier versions of TASSEL) often do not meet this second requirement and throw an error when TASSEL tries to load them. It can be difficult to recreate TASSEL’s internal sort order by hand, so this plugin allows the user to sort an input genotype file according to TASSEL’s rules and output it to a new file ready for further analysis. To sort a file from the command line, use the following command:
```
run_pipeline.pl -SortGenotypeFilePlugin -inputFile [filename] -outputFile [filename] -fileType [Hapmap or VCF]
```
The -fileType flag is optional and is only needed if the input file’s extension doesn’t match a known file extension (“.hmp.txt”, “.vcf”, etc.).

It gives me another error:

ERROR net.maizegenetics.pipeline.TasselPipeline - java.lang.IllegalArgumentException: TasselPipeline: parseArgs: Unknown parameter: -SortGenotypeFilePlugin


- [ ] Panel 2: waiting for the R packages installation
- [ ] Panel 1: ValidateVariants to check the validity of the gcvcfs

**9 September 2020**
- [x] Being updated about the ash project (Meeting with Jeremy last Friday)
- [ ] Writing Research plan for Smith postdoctoral fellowship (in operation, deadline 2 October)
- [ ] Recombination map (famlily assignment to obtain the pedigree data, Jeremy's code). This information will allow me to use [Lep-Map3 program](https://avikarn.com/2019-04-17-Genetic-Mapping-in-Lep-MAP3/) to perform the recombination map.


**2 September 2020**
- [ ] Imputation of Panel3 SNPs using Reference panel. High density recombination map: 

| Position | Combined_rate(cM/Mb) | Genetic_Map(cM) |
| ----------- | ----------- | ----------- |        
| Physical position (base pairs) | Recombination rate between current position and the next position in map (cM/Mb) | Genetic Map position (in cM) |

- Thinking about using MareyMap [R package](https://cran.r-project.org/web/packages/MareyMap/) to calculate the recombination rate: Problems with R package installation, unable to install the require "tcltk" package
- Trying using the MareyMap [online version](http://lbbe-shiny.univ-lyon1.fr/MareyMapOnline/). I need the genetic distances between markers (figuring out how to calculate them)
- Using [this article](https://www.g3journal.org/content/10/1/299#sec-1) as reference, which used two different methods (LD-based, **LDhelmet program**; and genetic-map based, **MareyMap**) to buid a recombination map in *Populus tremula*. They **compared both methods** and concluded that recombination rates obtained from the two methods largely agree, **although the LD-based method identifies a number of genomic regions with very high recombination rates that the map-based method fails to detect**. Genetic-map and LD-based estimates of recombination rates are positively correlated and show similar correlations with other genomic features, showing that **both methods can accurately infer recombination rate variation across the genome**. 

- [ ] Panel 2: waiting for the R packages installation
- [ ] Panel 1: ValidateVariants to check the validity of the gcvcfs

**10 August 2020**
- [x] Panel3: V4 positions obtained. Duplicates removed. After removing duplicates and merging in R: 24,702 SNPs, hmp format (Panel4SNPforGWASposv4.hmp.txt)
- [X] Reference panel: V4 positions obtained, duplicates removed: 1,653,460 SNPs, hmp format (ReferencepanelSNPposv4.hmp.txt)
- [ ] Imputation of Panel3 SNPs using Reference panel. High density recombination map: In operation, reading Imputev2 best practices to understand the input formats for study panel, reference panel and recombination map, and using [Jeremy's tutorial](https://avikarn.com/2019-04-17-Genetic-Mapping-in-Lep-MAP3/) to perform the recombination map
- [ ] Panel 2: waiting for the R packages installation
- [ ] Panel 1: ValidateVariants to check the validity of the gcvcfs (in operation)

**5 August 2020**
- [x] Panel3: V4 positions obtained. Duplicates removed. After removing duplicates: 31,171 SNPs. Genotyped SNPs for GWAS 29,355. After merging in R: 24,702 SNPs
- [x] Reference panel: Mapped Reference panel flanking regions, v4 positions obtained, removing duplicates (in operation)
- [ ] Imputation of Panel3 SNPs using Reference panel (High density recombination map needed, learning how to build it, avaible a "high density" recombination map for *P. trichocarpa*, 3568 SNPs [Muchero et al. 2015](https://doi.org/10.1186/s12864-015-1215-z))
- [ ] Panel 2: waiting for the R packages installation
- [ ] Panel 1: ValidateVariants to check the validity of the gcvcfs (in operation)

**3 August 2020**
- [x] Panel3: Incongruence between the number of positions in v3 (provided by the authors) and the number of positions after mapping the flanking regions against the v4 (some scaffolds duplicated, two different qualities, using "unique" in R, same number of unique scaffolds in both files, v3 and v4 positions: 31346): Selecting the ones with the higher quality. Doubt: 38 among them were just duplicated, same quality, same position. The other 192, differences in quality and position (I selected the ones showing higher quality 73M127S, 144H56M)
- [ ] Reference panel: Mapping Reference panel flanking regions
- [ ] Panel 2: waiting for the R packages installation
- [ ] Panel 1: ValidateVariants to check the validity of the gcvcfs (in operation)

**22 July 2020**
- [x] Reference panel: Problem solved (I talked with the authors), the reference genome used was v3 instead of v2
- [x] Reference panel/Panel3: Flanking regions from Reference panel, using v3; and panel 3, using the SNP possitions sent by Athena McKnown. I calculated 100 pb and 200 pb reads around the SNP locations in fasta format, to contrast (verify) the results. If my further calculations are correct, the results from both fasta files (SNP positions in v4.0) should be the same. 
- [X] Reads in fasta format mapped against the reference genome v4.0 using bwa.

Index the reference genome
```
module load bwa/0.7.17
bwa index Ptrichocarpa_533_v4.0.f
```
Mapping
```
bwa mem Ptrichocarpa_533_v4.0.fa Panel4flank100pb.fasta>Panel4flank100pb.sam
```
Sam to bam
```
samtools view -bhS Panel4flank100pb.sam > Panel4flank100pb.bam
```
- [x] Reference panel/panel3: Sam files after mapping against the v4.0 reference genome. Code to extract the positions:
Remove header from sam file (header lines start with the at sign @)
```
sed '/@/d' ./Panel4flank200pbprueba.sam>Panel4flank200pbprueba2.sam
```
In the alignment, 3rd column is the chromosome, 4th column is the position in which the flanking region start mapping against the reference genome. I got these two columns using the following code
```
awk '{print $3,$4}' Panel4flank200pbprueba2.sam > Panel4posv4.txt
```
How the result looks like
Chr01 5160045
Chr01 5160247
Chr01 5159789
Chr01 5159991

I kept only the even lines
```
awk 'NR%2==0' Panel4posv4.txt > Panel4posv4pares.txt
```
And finally, I substract "-1" from the even lines positions (second column). This way I obtain the v4 positions
```
awk '{print $1,$2 = $2 - 1}' Panel4posv4pares.txt > Panel4posv4paresdef.txt
```
I obtained the same positions using the two fasta files (100pb and 200pb flanking regions). Problem: The number of positions is longer that the positions of the initial file (v3 genome). I checked all the files and the increase in the number lines in the files happened in the sam file (after bwa mapping). In addition, the sam files resulting from the 100 pb fasta file and from the 200 pb fasta file, have different number of lines. I am investigating what happened. Idea: Change the headers of the fasta file by their position names in the original SNP panel to track what happen with each flanking region during the mapping

- [ ] Panel1: ValidateVariants to check the validity of the gcvcfs (in operation)
- [ ] Panel 2: no successful installation of dependencies from Computational Biology Core. I have asked them what it can be done then to have the R packages installed. They said that I have to drop an email to cbcsupport@uconn.edu for installation of the package specifying all the details, but at present Mike is busy with installation of new cores.

**21 July 2020**
- [x] Reference panel: Problem solved (I talked with the authors), the reference genome used was v3 instead of v2
- [x] Reference panel/Panel3: Flanking regions from Reference panel, using v3; and panel 3, using the SNP possitions sent by Athena McKnown. I calculated 25pb and 100 pb reads around the SNP locations in fasta format, to contrast (verify) the results. If my further calculations are correct, the results from both fasta files (SNP positions in v4.0) should be the same. 
- [X] Reads in fasta format mapped against the reference genome v4.0 using bwa.

Index the reference genome
```
module load bwa/0.7.17
bwa index Ptrichocarpa_533_v4.0.f
```
Mapping
```
bwa mem Ptrichocarpa_533_v4.0.fa Panel4flank100pb.fasta>Panel4flank100pb.sam
```
Sam to bam
```
samtools view -bhS Panel4flank100pb.sam > Panel4flank100pb.bam
```
- [ ] Panel1: ValidateVariants to check the validity of the gcvcfs (in operation)
- [ ] Panel 2: no successful installation of dependencies from Computational Biology Core. I have asked them what it can be done then to have the R packages installed. They said that I have to drop an email to cbcsupport@uconn.edu for installation of the package specifying all the details, but at present Mike is busy with installation of new cores.


**7 July 2020**
- [ ] Panel 1: ValidateVariants to check the validity of the gcvfs". Learning how to apply this code to check the validity of the gvcfs.
- [ ] Panel 2: Still waiting for the Computational Biology Core Slack answer regarding the R packages installation for DAPC analysis
- [x] Reference panel: Problem identified regarding the error I got from the SNP flanking regions code. I have checked the code with other files. The code is correct. The problem is that the chromosome 1 in P.trichocarpav2 has a lenght of 48367220. However, the last SNP in chromosome 1 from the Reference panel has its location in the position 50494826, which do not exist in the Ptrichocarpav2 genome. Dr Guerra told me that they used Ptrichocarpav2 to perform the SNP calling but I do not understand why these differences in lenght. I have written Dr Guerra and he said that Dr Jason Holliday and his team did the genotyping. Waiting for his answer. 
- [ ] Panel 3: Flanking regions. Doubts about the SNP positions. In the .csv file, the SNP names look like this (scaffold_1_10003096). I suppose that the SNP position is 10003096, located in the scaffold 1. But I need confirmation and I did not find this information in the reference article nor in the article describing the assay design (Geraldes et al. 2013). In addition, not sure if these SNP positions correspond to the v2.2, which was used for the assay design, or v3, since they say in the article that they remapped the SNP flanking regions against the Ptrichocarpa v3. I wrote the author (Dr. Athena McKnow) and she sent me the positions. I will try to use the code for the Reference panel to obtain the flanking regions and map them against the Ptrichocarpa v4
- [ ] Building a recombination map for *P. trichocarpa*

**29 June 2020**
- [x] Panel1: Still dealing with errors, now related to the GVCF files format. I have written to GATK community asking for help. Their answer: "The gvcfs are probably malformed. Use ValidateVariants to check the validity of the gcvfs". Learning how to apply this code to check the validity of the gvcfs.
- [x] Panel2: Discriminant analysis of principal components (DAPC) and kindship (using EMMAX) calculation to account for population structure and kindship during GWAS analysis. Problems during the package installation in R to perform DAPC analysis (packages "adegenet", "vcfR" and "poppr". I have reported my problem in Computational Biology Core Slack, they are helping me to install the packages. Apparently, “libudunits2.so”  is required for installation of some of dependent packages. They have submitted a request for it, still waiting for the answer)
- Reference panel and panel 3 (assay). Applying the code to get the flanking regions using the SNP positions: 
- [x] Bed file containing the Reference Panel SNP positions in v2.0 P.trichocarpa genome (successfully created from .hmp files, since the .hmp contain the SNPs positions)
- [x] Applying the following code to get the flanking regions. I am using bedtools flankBed to add 25 bases to the interval around each snp:
```
flankBed  -i Referencepanel.bed -g Populus_trichocarpa.v2.fa.gz -b 25 > flank.bed
```
- It is necessary a .genome file for the -g flag
- [x] .genome file succesfully created using the following code:
```
samtools faidx Populus_trichocarpa.v2.fa
cut -f1,2 Populus_trichocarpa.v2.fa.fai > sizesv2.genome
```
- Problem:  I have realized that the scaffold number from the Ptrichocarpav2 (3134) are different from the scaffold number from the Reference panel (3595). I have sent an email to Dr Guerra for clarifying the reason of this difference (waiting for the answer)
- In the meantime, I have removed the scaffolds and only keep the 19 chromosomes for further analysis in both .genome file and .bed file
- Head of Referencepanelonlychr.bed
```
chr01	52	52
chr01	117	117
chr01	133	133
chr01	2387	2387
chr01	3111	3111
chr01	3112	3112
chr01	3126	3126
chr01	3211	3211
chr01	3241	3241
```
- Head of sizesv2chr.genome
```
chr01	48367220
chr02	23562801
chr03	20216605
chr04	23188140
chr05	25802683
chr06	26894541
chr07	15101417
chr08	18835763
chr09	12942059
chr10	21538349
```
- I run the following code

```
flankBed  -i Referencepanelonlychr.bed -g sizesv2chr.genome -b 25 > flankchr.bed
```
- I obtained the following flankchr.bed file
```
chr01	26	51
chr01	53	78
chr01	91	116
chr01	118	143
chr01	107	132
chr01	134	159
chr01	2361	2386
chr01	2388	2413
chr01	3085	3110
chr01	3112	3137
```
- Finally I run the following code to retrieve these regions from the Populus_trichocarparv2.fa reference genome in fasta format 

```
bedtools getfasta -fi Populus_trichocarpav2.fa -bed flankchr.bed -fo flankchr.fasta
```
Error: Feature (scaffold_1:48367428-48367453) beyond the length of scaffold_1 size (48367220 bp).  Skipping.
Error: malformed BED entry at line 404232. Start was greater than end. Exiting.

- I have checked the Reference panel article (Guerra et al. 2019) just in case they used another version of the reference genome, but this is not the problem apparently (Phytozome version 7.0 annotation and assembly files for P. trichocarpa (corresponding to assembly version 2.0 of the black cottonwood genome))

- [ ] Building a recombination map for *P. trichocarpa*

**24 June 2020**
- [x] GATK Haplotype caller in the error sample: Successfully rerun, now consolidating the resulting .vcf files. Error due to the format of the interval flag: Invalid interval record contains 1 fields: chr1, for input source: chromosomes.interval_list
 ```
 LIST=chromosomes.interval_list
gatk --java-options "-Xmx10g -Xms4g" GenomicsDBImport \
    -V 206.sort.mkdup.g.vcf \
    -V 207.sort.mkdup.g.vcf \
    -V 208.sort.mkdup.g.vcf \
    -V 209.sort.mkdup.g.vcf \
    -V 210.sort.mkdup.g.vcf \
    --genomicsdb-workspace-path Panel1 \
    --overwrite-existing-genomicsdb-workspace true \
    --intervals $LIST
```
```
chr1
chr2
chr3
chr4
chr5
chr6
chr7
chr8
chr9
chr10
chr11
chr12
chr13
chr14
chr15
chr16
chr17
chr18
chr19
```

```
gatk --java-options "-Xmx10g -Xms4g" GenomicsDBImport \
    -V 206.sort.mkdup.g.vcf \
    -V 207.sort.mkdup.g.vcf \
    -V 208.sort.mkdup.g.vcf \
    -V 209.sort.mkdup.g.vcf \
    -V 210.sort.mkdup.g.vcf \
    --genomicsdb-workspace-path Panel1 \
    --overwrite-existing-genomicsdb-workspace true \
    --L chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19
```
A USER ERROR has occurred: Badly formed genome unclippedLoc: Query interval "chr1" is not valid for this input.

- [ ] Panel2: Discriminant analysis of principal components (DAPC)(problem, it can be only run in R, I tried in my computer, lack of memory for running the .vcf file), I decided to use the PCs from the PCA as covariates in EMMAX), and kindship (using EMMAX) calculation to account for population structure and kindship during GWAS analysis (in operation, understanding the code). PC1-PC4 cumulatively explain the 49.9% of the variance, should I choose them to account for population structure in EMMAX? Or better follow the Evans et al. 2014 method: To account for population structure, for each trait we included as covariates the principal-component axes that were significant predictors of the trait, chosen using stepwise regression (model selection based on Akaike's Information Criterion in the step function in the R package). 
```
calculate the cumulative sum of the percentage variance explained
> cumsum(pve$pve)
 [1]  21.82813  32.67095  42.67034  49.91037  55.42560  60.33279  63.93155  67.24105  70.42995
[10]  73.51992  76.46584  79.26046  82.02138  84.72062  87.40212  90.04829  92.58906  95.11887
[19]  97.59661 100.00000
```

- [ ] Reference panel and panel 3 (assay). Figuring out the code to get the flanking regions using the SNP positions (in operation)
- [ ] Build a recombination map for *P. trichocarpa*

**22 June 2020**
- [x] GATK Haplotype caller in the error sample: Successfully rerun, now consolidating the resulting .vcf files. 
- [x] GWAS in Panel 2: Filtered by LD in plink (threshold r2<0.2) to retain only independent SNPs for population structure and kindship analyses.
- [x] PCA performed [here](https://gitlab.com/TreeGenes/CartograTree/-/blob/cartogratrees_master/galaxy_workflows/Imputation-based.meta-analysis/PCAPanel2location.png). Very similar to the obtained in the reference article [Evans et al. 2014](https://doi.org/10.1038/ng.3075). They exclude the highly differentiated Tahoe, Willamette Valley, and far northern British Columbia samples because strong stratification can lead to spurious associations (even if I account for population structure in the further GWAS analysis???). Ideas: Maybe I can perform the GWAS with and without these samples and see how the look like or perform the GWAS separatedly for Tahoe, Willamette Valley, far northern British Columbia samples and the rest of samples, since I think it is a shame to skip these highly differentiated populations, they might be, precisely because of their strong genetic differentiation, a good source of genetic diversity/adaptative genetic variants?
- [ ] Discriminant analysis of principal components (DAPC) and kindship (using EMMAX) calculation to account for population structure and kindship during GWAS analysis (in operation)
- [ ] Reference panel and panel 3 (assay). Figuring out the code to get the flanking regions using the SNP positions (in operation)
- [x] Instructions manual for the workflow (with all the checks, input/output formats, programs and flags to include in TreeGenes Galaxy) [here](https://gitlab.com/TreeGenes/CartograTree/-/blob/cartogratrees_master/galaxy_workflows/Imputation-based.meta-analysis/Imputation-based_meta-analysis_workflow_for_GWAS_and_GEA_analysis_in_CartograTree.md)
- [ ] Build a recombination map for *P. trichocarpa*
- [x] Wild relatives added to this [page](https://docs.google.com/spreadsheets/d/1ZMcREgdogpOPbW_aWdD2TL7QkdszAdqSjwNBH-DRXzA/edit#gid=1963599508) (CartograPlant/Wildrelatives) 


**17 June 2020**
- [x] Running GATK Haplotype caller in the error sample

- [x] GWAS in Panel2 (in operation using this [tutorial](http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/vcf-activity/)). Quality filtering is being performed sequencially (1. Missing data, 2. Frequency, 3. Quality). Missing data 10% threshold already run (21,694,149 out of a possible 28,342,758 sites), MAF 0.05 (5,841,562 out of 28,342,758 sites), mac3 (5,841,562 out of 5,841,562 sites), minDP3 (5,841,562 out of 5,841,562 sites), minQ20 (5,841,562 out of 5,841,562 sites)in this order. Doubt: When I did the analysis running all the filters together, with the same thresholds, I obtained less SNPs (5462016 out of 28342758 Sites): vcftools --gzvcf SNP.vcf.gz --mac 3 --minQ 20 --minDP 3 --max-missing 0.9 --maf 0.05 --recode --recode-INFO-all --out SNP.mac3.minQ20.minDP3.miss10.maf0.05. Filter again max-missing 0.9 and I obtained 5474907 out of 5841562 Sites, still a bit more than using all the filters together. Shoud I use the most restrictive one for LD and population structure analysis?
- [x] GWAS in Panel 2: following steps: Filter by LD in plink (threshold r2<0.2) to retain only independent SNPs for population structure and kindship analyses.
- [ ] Performing PCA and the discriminant analysis of principal components (DAPC) and calculate kindship (using EMMAX) to account for population structure and kindship during GWAS analysis (in operation)
- [ ] Reference panel and panel 3 (assay). Figuring out the code to get the flanking regions using the SNP positions (in operation)
- [x] Instructions manual for the workflow (with all the checks, input/output formats, programs and flags to include in TreeGenes Galaxy) [here](https://gitlab.com/TreeGenes/CartograTree/-/blob/cartogratrees_master/galaxy_workflows/Imputation-based.meta-analysis/Imputation-based_meta-analysis_workflow_for_GWAS_and_GEA_analysis_in_CartograTree.md)
- [ ] Build a recombination map for *P. trichocarpa*
- [x] Wild relatives added to this [page](https://docs.google.com/spreadsheets/d/1ZMcREgdogpOPbW_aWdD2TL7QkdszAdqSjwNBH-DRXzA/edit#gid=1963599508) (CartograPlant/Wildrelatives) 


**15 June 2020**
- [x] Running GATK Haplotype caller in some samples  It is running using 12 cores, 2 hours per sample more or less (running the last one): gatk --java-options  "-Xmx16g -XX:ParallelGCThreads=1" HaplotypeCallerSpark --spark-master local[12] -R Ptrichocarpa_533_v4.0.fa -I 213.sort.mkdup.bam -O 213.sort.mkdup.g.vcf --emit-ref-confidence GVCF. (finished but I found an error in one of them: ERROR Utils: Uncaught exception in thread Executor task launch worker for task 46
java.lang.NullPointerException
	at org.apache.spark.scheduler.Task$$anonfun$run$1.apply$mcV$sp(Task.scala:142)
	at org.apache.spark.util.Utils$.tryLogNonFatalError(Utils.scala:1340)
	at org.apache.spark.scheduler.Task.run(Task.scala:140)
	at org.apache.spark.executor.Executor$TaskRunner$$anonfun$10.apply(Executor.scala:408)

- [x] GWAS in Panel2 (in operation using this [tutorial](http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/vcf-activity/)). Quality filtering is being performed sequencially (1. Missing data, 2. Frequency, 3. Quality). Missing data 10% threshold already run (21,694,149 out of a possible 28,342,758 sites), MAF 0.05 (5,841,562 out of 28,342,758 sites), mac3 (5,841,562 out of 5,841,562 sites), minDP3 (5,841,562 out of 5,841,562 sites), minQ20 (5,841,562 out of 5,841,562 sites)in this order. Doubt: When I did the analysis running all the filters together, with the same thresholds, I obtained less SNPs (5462016 out of 28342758 Sites): vcftools --gzvcf SNP.vcf.gz --mac 3 --minQ 20 --minDP 3 --max-missing 0.9 --maf 0.05 --recode --recode-INFO-all --out SNP.mac3.minQ20.minDP3.miss10.maf0.05. Filter again max-missing 0.9 and I obtained 5474907 out of 5841562 Sites, still a bit more than using all the filters together.
- [x] GWAS in Panel 2: following steps: Filter by LD in plink (threshold r2<0.2) to retain only independent SNPs for population structure and kindship analyses, perform the discriminant analysis of principal components (DAPC) and calculate kindship (using EMMAX) to account for population structure and kindship during GWAS analysis (in operation)
- [x] Reference panel and panel 3 (assay). Following IMPUTEv2 best practices: Getting SNP flanking regions for Reference panel SNPs to remap them against P.trichocarpa v4.0. 1. Reference panel in .hmp format, I transformed it to .vcf format using Tassel: ./tassel5-standalone/run_pipeline.pl -Xmx5g -fork1 -h Poplar_Genotype_01.hmp.txt -export -exportType VCF -runfork1 2. I have 1 hmp file per chromosome, I want to concantenate them in a single file, thus, I compressed them to vcf.gz (bcftools view) and index them (bcftools index). I tried to concatenate them using bcftools concate but it gave me an error (Checking the headers and starting positions of 6 files
[E::hts_hopen] Failed to open file Poplar_Genotype_01.in.vcf.gz
[E::hts_open_format] Failed to open file "Poplar_Genotype_01.in.vcf.gz" : Exec format error
Failed to open: Poplar_Genotype_01.in.vcf.gz)
I tried transforming one of the vcf files to plink format using "vcftools --vcf input_data.vcf --plink --out out_in_plink" and getting the flanking regions using bedtools "slopBed -i Poplar_Genotype_01.bed -g Ptrichocarpa_533_v4.0.fa.fai -b 10 > chr01flankregions10bases.bed", but it also gave me an error: It looks as though you have less than 3 columns at line 1 in file Poplar_Genotype_01.bed.  Are you sure your files are tab-delimited?
- [ ] Build a recombination map for *P. trichocarpa*
- [x] Wild relatives added to this [page](https://docs.google.com/spreadsheets/d/1ZMcREgdogpOPbW_aWdD2TL7QkdszAdqSjwNBH-DRXzA/edit#gid=1963599508) (CartograPlant/Wildrelatives) 

**11 June 2020**
- [x] Running GATK Haplotype caller in some samples  It is running using 12 cores, 2 hours per sample more or less (running the last one): gatk --java-options  "-Xmx16g -XX:ParallelGCThreads=1" HaplotypeCallerSpark --spark-master local[12] -R Ptrichocarpa_533_v4.0.fa -I 213.sort.mkdup.bam -O 213.sort.mkdup.g.vcf --emit-ref-confidence GVCF. 
- [x] GWAS in Panel2 (in operation using this [tutorial](http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/vcf-activity/)). Quality filtering is being performed sequencially (1. Missing data, 2. Frequency, 3. Quality). Missing data 10% threshold already run (21,694,149 out of a possible 28,342,758 sites), now running MAF 0.05, following steps: mac3, minDP3, minQ20 in this order.
- [ ] GWAS in Panel 2: following steps: Filter by LD (threshold r2<0.2) to retain only independent SNPs for population structure and kindship analyses, perform the discriminant analysis of principal components (DAPC) and calculate kindship (using EMMAX) to account for population structure and kindship during GWAS analysis
- [x] Reference panel and panel 3 (assay). Following IMPUTEv2 best practices: Getting SNP flanking regions for Reference panel SNPs to remap them against P.trichocarpa v4.0 (in operation) (thinking about different options [this](https://www.researchgate.net/post/How_get_SNP_flanking_sequence), [this](https://www.biostars.org/p/68344/), [this](https://www.biostars.org/p/391066/) or [this](https://www.itb.cnr.it/web/bioinformatics/g-snpm) [reference article](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4015528/)
- [x] Workflow uploaded with the raw sequences steps (including checks)
- [ ] Build a recombination map for *P. trichocarpa*
- [x] Wild relatives added to this [page](https://docs.google.com/spreadsheets/d/1ZMcREgdogpOPbW_aWdD2TL7QkdszAdqSjwNBH-DRXzA/edit#gid=1963599508) (CartograPlant/Wildrelatives) 

**8 June 2020**
- [x] Running GATK Haplotype caller in some samples  It is running using 12 cores, 2 hours per sample more or less (still running): gatk --java-options  "-Xmx16g -XX:ParallelGCThreads=1" HaplotypeCallerSpark --spark-master local[12] -R Ptrichocarpa_533_v4.0.fa -I 213.sort.mkdup.bam -O 213.sort.mkdup.g.vcf --emit-ref-confidence GVCF. 
- [x] GWAS in Panel2 (in operation using this [tutorial](http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/vcf-activity/)). Now, quality filtering: MAF 0.05, missing data 10%, and I have included: mac 3 (remove SNPs with minor allele count less than 3), minQ 20 (minimun quality score of 20), minDP 3 (recode genotypes with less than 3 reads). After filtering, all individuals were kept and 5,462,016 out of a possible 28,342,758 Sites (maybe should I change the threshold?) In addition, 1.	Perform an histogram of % of missing data (X axes) per Number of occurrences (Y axes)  to select the individuals you want to remove in the following step. A good threshold could be to remove individuals with more than 25% of missing data. But the histogram will give you an idea. If the majority of the individuals show less than 0.5 missing data, go for this cutoff. 
- [x] GWAS in Panel 2: following steps: Filter by LD (threshold r2<0.2) (running now) to retain only independent SNPs for population structure and kindship analyses, perform the discriminant analysis of principal components (DAPC) and calculate kindship (using EMMAX) to account for population structure and kindship during GWAS analysis
- [x] Reference panel and panel 3 (assay): Transform them into IMPUTE v2 format (reading the IMPUTEv2 best practices)
- [x] Workflow uploaded with the raw sequences steps (including checks)
- [ ] Look for a recombination map for *P. trichocarpa*
- [ ] Look for georeferenced SNP-based articles, wild relatives for apple

**1 June 2020**
- [x] Running GATK Haplotype caller in some samples (running time issues). Trying with multithread flag using this code: gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=12" HaplotypeCaller -R Ptrichocarpa_533_v4.0.fa -I 210.sort.mkdup.bam -ERC GVCF --output 210.sort.mkdup.g.vcf (super slow). Jill's code: gatk --java-options  "-Xmx16g -XX:ParallelGCThreads=1" HaplotypeCallerSpark --spark-master local[2] -R Ptrichocarpa_533_v4.0.fa -I 213.sort.mkdup.bam -O 213.sort.mkdup.g.vcf --emit-ref-confidence GVCF (faster, but 5 hours per sample). I am trying using 12 cores: gatk --java-options  "-Xmx16g -XX:ParallelGCThreads=1" HaplotypeCallerSpark --spark-master local[12] -R Ptrichocarpa_533_v4.0.fa -I 213.sort.mkdup.bam -O 213.sort.mkdup.g.vcf --emit-ref-confidence GVCF. If it is still slow, I will try performing the SNP calling on individual chromosomes and scaffolds to reduce the run time.
- [x] GWAS in Panel2 (in operation using this [tutorial](http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/vcf-activity/)). Now, quality filtering: MAF 0.05, missing data 10% or 5%
- [x] Reference panel received in hmp format for 19 chromosomes separately and 1 scaffold (learning about hmp format). Objective: Transform it into .gens format for IMPUTE v2
- [x] Workflow uploaded with the raw sequences steps (including checks)
- [ ] Look for a recombination map for *P. trichocarpa*
- [ ] Look for georeferenced SNP-based articles, wild relatives for apple


**26 May 2020**
- [x] Re-run bwa with the -R flag, samtools (.sam to .bam) and picard sortSam and MarkDuplicates
- [ ] Realign indels using GATK indel realigner and base recalibration (GATK) (GATK best practices). Afterwards, SNP calling using GATK
- [x] Best linear unbiased predictors (BLUPs) (phenotypes) calculated for Panels 2 and 3 (I already had the BLUPs for Panel 1) (doubts about the model used in lmer4, how I handle NAs in both panels and the absence of multiple common garden locations in panel 3: **Panel 2** hmodel = lmer(HEIGHT ~ (1|LINE) + (1|LOC) + (1|YEAR) + (1|REP%in%LOC:YEAR) + (1|LINE:LOC) + (1|LINE:YEAR),na.action = na.omit)), **Panel 3** (one only location for common garden): hmodel = lmer(HEIGHT ~ (1|LINE) + (1|YEAR) + (1|REP:YEAR) + (1|LINE:YEAR),na.action = na.omit)
- [x] GWAS in Panel2 (in operation, learning how to perform a PCA with SNPs from a .vcf file using this [tutorial](http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/vcf-activity/))
- [ ] Waiting for the reference panel upload to TreeGenes (Guerra et al. 2019)
- [ ] Look for a recombination map for *P. trichocarpa*

**20 May 2020**
- [x] Panel 1: Reading the GATK tutorial, doubt, should I have used the "-R flag" to add the sample ID in the bwa code to align fastq.gz samples against the reference genome? This is the code I used "bwa mem Ptrichocarpa_533_v4.0.fa SRR1819206_1.fastq.gz SRR1819206_2.fastq.gz>206.sam" . The manual says that it does not have to be done at this stage (picard has a tool for adding read groups to alignment files). 
- [x] Panel 1: sort the reads in the alignment files by their positions in the reference genome (Picard SortSam function) (still running)
- [ ] Panel 1: Remove duplicate reads (Picard MarkDuplicates function)
- [x] Reference panel: SNP datafile will be uploaded through this week
- [x] Estimating Best linear unbiased predictors (BLUPs) for panel 2 (in operation)

**18 May 2020**
- [x] Panel 1: Still mapping against the latest reference genome (v4). Afterwards: samtools to convert .sam to .bam files and duplicate reads removal (in the article they used Picard and GATK to do so, I have previously used samtools to perform this duplicate reads removal, is there any substancial difference? 
- [x] Reference panel: Problem solved, I talked with the authors (actually there was only 39 individuals updated to SRA, they are a subset from the 461 clones included in the association population. These were sequenced in the first stage for constructing the exome libraries (and SNP discovery). Then, these libraries/information were utilized for genotyping the entire population. The entire population SNP database is going to be uploaded to TreeGenes)
- [x] Estimating Best linear unbiased predictors (BLUPs) for each trait of interest (heigh, bud set and bud flush) from panel 2 using lmer4 in R (learning the [code](https://pbgworks.org/sites/pbgworks.org/files/TBRT2011Script.txt), [.csv file phenotype quality tomato](https://pbgworks.org/sites/pbgworks.org/files/TBRTQuality.csv)).
- [x] Compilation of a list of literature related to wild crop relatives: Looking for SNP-based studies using geolocated wild crop relatives of apple (Malus domestica), one of the first crop genomes sequenced: M. floribunda, M. baccata jackii, M. micromalus (interest in Apple scab (Venturia inaequalis) resistance). Thinking also about cacao and banana (reference genomes as well, but not sure about the amount of genomic resources). Reference article: https://doi.org/10.3389/fpls.2017.00460 
- [ ] Recombination map latests poplar genetic map [2018](https://link.springer.com/article/10.1007/s00425-018-2958-y#Sec15]) (889 and 1650 SNPs formed the female and male genetic maps, respectively.) and [2019](https://link.springer.com/article/10.1007/s00468-018-01809-y#Abs1), 913 female 252 male respectively or could I use a reference genome as recombination map, what is the difference? Or would it be better to build one? There should be a [high-density genetic map for P. trichocarpa](https://phytozome.jgi.doe.gov/pz/portal.html#!info?alias=Org_Ptrichocarpa) out there but I do not manage to find it
- [ ] Visual Imputation-based meta-analysis workflow (in operation)

**13 May 2020**
- [x] Data quality filtering of panel 1 (fastqc in operation and afterwards, multiqc to summarize the results) and try to perform GWAS with SNPs in panel 2 (.vcf file)
- [x] Reference panel, problem with the low number of individuals (39): I have read the [Bioproject](https://www.ncbi.nlm.nih.gov/sra?term=SRA058855&cmd=DetailsSearch) and the article results, nothing found 
- [x] Recombination map latests poplar genetic map [2018](https://link.springer.com/article/10.1007/s00425-018-2958-y#Sec15]) (889 and 1650 SNPs formed the female and male genetic maps, respectively.) and [2019](https://link.springer.com/article/10.1007/s00468-018-01809-y#Abs1), 913 female 252 male respectively or could I use a reference genome as recombination map, what is the difference? Or would it be better to build one? There should be a [high-density genetic map for P. trichocarpa](https://phytozome.jgi.doe.gov/pz/portal.html#!info?alias=Org_Ptrichocarpa) out there but I do not manage to find it
- [ ] Mapping raw reads against the reference genome using BWA (Panel 1, v4.1, and reference panel, v2.2) 
- [ ] SNP calling of Panel 1 and reference panel using GATK
- [ ] Visual Imputation-based meta-analysis workflow (in operation)
- [ ] Compilation of a list of literature related to wild crop relatives (in operation)

**11 May 2020**
- [x] Data quality filtering of panel 1 (problems with reference panel) and try to perform GWAS with SNPs in panel 2 (.vcf file)
- [x] Mapping raw reads against the reference genome (Panel 1, v4.1, and reference panel, v2.2) 
- [ ] Recombination map (in operation)
- [ ] Data filtering of Genome-wide Exome-wide panels (in operation)
- [ ] Visual Imputation-based meta-analysis workflow (in operation)
- [ ] Compilation of a list of literature related to wild crop relatives (in operation)

**6 May 2020**
- [x] Mapping raw reads against the reference genome (Panel 1, v4.1, and reference panel, v2.2) 
- [ ] Recombination map (in operation)
- [ ] Data filtering of Genome-wide Exome-wide panels (in operation)
- [ ] Visual Imputation-based meta-analysis workflow (in operation)
- [ ] Compilation of a list of literature related to wild crop relatives (in operation)

**4 May 2020**
- [x] Reference panel selection
- [x] Workflow update with things to consider during reference panel selection
- [x] Reference panel download
- [x] P.trichocarpa reference genomes download (v2.0, 2.2, 3.0, 4.1)
- [ ] Mapping raw reads against the reference genome (Panel 1 ?? and reference panel, v2.2)
- [ ] Recombination map (in operation)
- [ ] Data filtering of Genome-wide Exome-wide panels (in operation)
- [ ] Visual Imputation-based meta-analysis workflow (in operation)
- [ ] Compilation of a list of literature related to wild crop relatives (in operation)

**29 April 2020**
- [x] Phenotype table: Different criteria for measuring phenotypes through panels
- [x] Panel locations map in R with sequencing (different shapes): Some individuals overlap
- [ ] Reference panel (in operation)
- [ ] Recombination map (in operation)
- [ ] Data filtering of  WGS sequences (in operation)
- [ ] Visual Imputation-based meta-analysis workflow (in operation)
- [ ] Compilation of a list of literature related to wild crop relatives (in operation)

**22 April 2020**
- [ ] Visual Imputation-based meta-analysis workflow (in operation)
- [x] Panel selection
- [x] Panel locations map in R
- [x] Table with Panel information
- [x] Get SNPs, phenotype and environmental panel information
- [x] Be familiar with Gitlab (Guide gitimmersion.com)
- [x] Get assay design information for Panel 3
- [ ] Recombination map
- [ ] Reference panel
- [ ] Compilation of a list of literature related to wild crop relatives
