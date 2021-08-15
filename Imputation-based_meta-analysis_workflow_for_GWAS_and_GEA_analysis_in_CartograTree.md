# Imputation-based meta-analysis workflow for GWAS and GEA analyses in CartograTree

## Index
1. [Mapping raw reads](#rawreads) (start here if you have **.fastq.gz** files)
2. [Data quality filtering](#dat) (start here if you have high density SNPs -**WGS or Exome capture**)
3. [Imputation](#imp) (start here if you have low density SNPs -**assay**s-)

<a name="rawreads"></a>
### 1. Mapping raw reads

Start here if you have raw **fastq.gz** files

**1.1.** [Multiqc](https://multiqc.info/docs/ ). Use trimmomatic if it is necessary to remove adapters or overrepresented sequences.

**1.2.** Mapping against the reference genome using **bwa**

```
bwa index Ptrichocarpa_533_v4.0.fa
```
 (index reference genome)

 ```
 module load bwa/0.7.17
bwa mem -t 12 -R '@RG\tID:206\tSM:206' Ptrichocarpa_533_v4.0.fa SRR1819206_1.fastq.gz SRR1819206_2.fastq.gz>206.sam

```

Flags:

- "-t": indicates the number of CPUs
- "-R": individual label
- Afterwards: Reference genome in fasta format (Ptrichocarpa_533_v4.0.fa) and the two pair end reads in *fastq.gz* format (SRR1819206_1.fastq.gz and SRR1819206_2.fastq.gz). The output is a *.sam* file (206.sam).

**1.3.** Transform *.sam* to *.bam* (**samtools**)

```
module load samtools/1.10
samtools view -bhS 206.sam > 206.bam

```

**1.4.**	**Picard** sort sam (sort the reads in the chromosome order) (only modify input and output names, the rest remain as default)

```
module load picard/2.9.2
export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch
java -jar /isg/shared/apps/picard/picard-tools-2.9.2/picard.jar SortSam INPUT=207.bam OUTPUT=207.sort.bam SORT_ORDER=coordinate CREATE_INDEX=True

```

**1.5**. Mark duplicate reads produced during PCR (**picard**) (only modify the input, output and metricfile names, the rest remains as default)

```
module load picard/2.9.2
export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch
java -jar /isg/shared/apps/picard/picard-tools-2.9.2/picard.jar MarkDuplicates INPUT=575.sort.bam OUTPUT=575.sort.mkdup.bam METRICS_FILE=575.metrics.txt ASSUME_SORT_ORDER=coordinate CREATE_INDEX=True

```

**1.6.**	Index reference genome (**samtools**)

```
module load samtools/1.10
samtools faidx Ptrichocarpa_533_v4.0.fa

```

**1.7.** SNP calling (**GATK**)

```
module load GATK/4.1.3.0
gatk --java-options "-Xmx16g -XX:ParallelGCThreads=1" HaplotypeCallerSpark --spark-master local[12] -R Ptrichocarpa_533_v4.0.fa -I 206.sort.mkdup.bam -O 206.sort.mkdup.g.vcf --emit-ref-confidence GVCF

```

Flags

- "--spark-master local[12]": Use 12 CPUs
- "-R": Reference genome
- "-I": Input (bam files)
- "-O": output (vcf files)
- Leave the other two flags (""-Xmx16g -XX:ParallelGCThreads=1" and "--emit-ref-confidence GVCF") as default

**1.8.** Group *g.vcf* files together 

[GATK 4.1.8.1 documentation](https://gatk.broadinstitute.org/hc/en-us/sections/360009803432-4-1-8-1?page=6#articles)

[GATK 4.1.8.1 CombineGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360047217791-CombineGVCFs)
```
module load GATK/4.1.8.1
#IMPORTANT: The -Xmx value the tool is run with should be less than the total amount of physical memory available by at least a few GB, as the native TileDB library requires additional memory on top of the Java memory. Failure to leave enough memory for the native code can result in confusing error messages!
gatk --java-options "-Xmx80g -Xms80g" CombineGVCFs -R Ptrichocarpa_533_v4.0.fa --variant 206.sort.mkdup.g.vcf --variant 207.sort.mkdup.g.vcf --variant 208.sort.mkdup.g.vcf -O 206.207.208.g.vcf
```
Flags

- "-R": Reference genome in fasta format
- "--variant": Inputs (single g.vcf or g.vcf.gz files)
- "-O": output (combined g.vcf or g.vcf.gz file)
- --java-options "-Xmx80g -Xms80g": memory used (80 GB in this case)

Working directory: /labs/Wegrzyn/IreneCobo/Ptrichocarpa/Panel1-3/FASTQ

Code: GATKcombineGVCF.sh

**1.9.** Transform the resulting **g.vcf** file containing all the individuals into a **vcf** file

[GATK 4.1.8.1 documentation](https://gatk.broadinstitute.org/hc/en-us/sections/360009803432-4-1-8-1?page=6#articles)

[GATK 4.1.8.1 GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360047218551-GenotypeGVCFs)

```
module load GATK/4.1.8.1

#IMPORTANT: The -Xmx value the tool is run with should be less than the total amount of physical memory available by at least a few GB, as the native TileDB library requires additional memory on top of the Java memory. Failure to leave enough memory for the native code can result in confusing error messages!
gatk --java-options "-Xmx80g -Xms80g" GenotypeGVCFs -R Ptrichocarpa_533_v4.0.fa -V 206.207.208.g.vcf -O 206.207.208.vcf
```
Flags

- "-R": Reference genome in fasta format
- "-V": Input (combined g.vcf or g.vcf.gz file)
- "-O": output (combined .vcf or .vcf.gz file)
- --java-options "-Xmx80g -Xms80g": memory used (80 GB in this case)

Working directory: /labs/Wegrzyn/IreneCobo/Ptrichocarpa/Panel1-3/FASTQ

Code: GATKGenotypeGVCF.sh

<a name="dat"></a>
### 2. Data quality filtering

Start here if you have a *.vcf.gz* containing **high density SNPs** (WGS or exome capture). If the SNPs are in a different format, transform them into a *.vcf* file and compress it to a *.vcf.gz* file using the following code in **bcftools**. Data quality filtering for vcf files based on this [tutorial](http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/vcf-activity/):     

```
module load bcftools/1.10.2
bcftools view SNP.vcf -Oz -o SNP.vcf.gz
```

Once obtained the .vcf file after SNP calling or from your selected panels:

**2.1**.	Perform an histogram of % of missing data (X axes) per Number of occurrences (Y axes)  to select the individuals you want to remove in the following step (**it would be so good an interactive histogram allowing the user to perform the filtering directly on it**). A good threshold could be to remove individuals with more than 25% of missing data. But the histogram will give you an idea. If the majority of the individuals show less than 0.5 missing data, go for this cutoff. The individuals you want to remove can be also obtained using the following code:

```
module load vcftools/0.1.16
vcftools --gzvcf SNP.vcf.gz --missing-indv --out SNP.ed
```

This way, you will obtain a list of your individuals with their frequency of missing data (fifth column):

INDV N_DATA N_GENOTYPES_FILTERED N_MISS F_MISS

BS001 1184094 0 53664 0.0453207
BS002 1184094 0 124941 0.105516
BS003 1184094 0 84133 0.0710526

You can sort this fifth column from high to small values and select the individuals with a frequency of missing data higher than 25 %

```
cat data.ed.imiss | sort -rk5
```

**2.2.**	Perform the filtering using **VCFtools**. Include the individuals you want to remove using the following flag (--remove-ind):

```
 vcftools --gzvcf SNP.vcf.gz --remove-ind individual_name
 ```
**2.3.** Perform the rest of data quality filtering using **VCFtools**:

```
module load vcftools/0.1.16
vcftools --gzvcf SNP.vcf.gz --mac 3 --minQ 20 --minDP 3 --max-missing 0.9 --maf 0.05 --recode --recode-INFO-all --out SNP.mac3.minQ20.minDP3.miss10.maf0.05
```

Flags:

-	mac 3 (remove SNPs with minor allele count less than 3)
-	minQ 20 (minimun quality score of 20)
-	minDP 3 (recode genotypes with less than 3 reads)
-	max missing 0.9 (remove SNPs with more than 0.10 frequency missing individuals per genotype)
-	maf 0.05 (remove SNPs with an allele frequency lower than 0.05). This is the standard for high density SNP panels (WGS or exome capture). For assays, follow the formula: Number of chromosomes/(2 × Number of individuals) (Aulchenko, Y. S., Ripke, S., Isaacs, A. & van Duijn, C. M. GenABEL: an R library for genome-wide association analysis. Bioinformatics 23, 1294–1296 (2007).)

It is advisable to perform this quality filtering sequencially following the same flag order than in the previous code. This way, it is possible to check how many SNPs are retained after each filtering step

This data quality filtering can be also performed in **plink**

Removing bad SNPs and individuals: First, remove any individuals who have less than, say, 95% genotype data (--mind 0.05); and then remove SNPs that have less than, say 1% minor allele frequencies (--maf 0.01); and then remove SNPs that have less than, say, < 90% genotype call rate or >10% genotype error rate (--geno 0.1).
removing individuals with genotyping error >5% and SNPs with maf <1% and genotype missing data <90% and SNPs with pvalues < 0.05 of deviation from HWE :

```
plink --vcf SNP.vcf.gz --mind 0.10 --maf 0.05 --geno 0.05 --hwe 0.05 --recode vcf --out SNPclean
```

**2.3.**	Only for population structure and kindship analysis (we want to account for population structure and kinship during the GWAS and GEA), filter by linkage disequilibrium (LD) (we only want independent SNPs). Afterwards, we will perform a PCA (we are going to use PLINK). Based on [this tutorial](https://speciationgenomics.github.io/pca/). Steps:

-	Transform the *.vcf.gz* or *.vcf* file to plink format (.ped) using **vcftools** and the following code:

```
vcftools --gzvcf SNP. mac3.minQ20.minDP3.miss10.maf0.05.vcf.gz --plink --out SNP. mac3.minQ20.minDP3.miss10.maf0.05
```

-	Filter by LD (threshold r2<0.2) to retain only independent SNPs (**plink**)

```
plink --file SNP.mac3.minQ20.minDP3.miss10.maf5 --indep-pairwise 50 5 0.2 --out SNP.mac3.minQ20.minDP3.miss10.maf0.05.LD0.2
```

or using directly the .vcf file (not necessary to previously transform it into plink format)

```
module load plink/1.90.beta.4.4
plink --vcf SNP.mac3.minQ20.minDP3.miss10.maf0.05.recode.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 5 0.2 --recode vcf --out SNP.mac3.mi$
```

Flags:

- --vcf - specified the location of our VCF file.
- --double-id - told plink to duplicate the id of our samples (this is because plink typically expects a family and individual id - i.e. for pedigree data - this is not necessary for us.
- --allow-extra-chr - allow additional chromosomes beyond the human chromosome set. This is necessary as otherwise plink expects chromosomes 1-22 and the human X chromosome.
- --set-missing-var-ids - also necessary to set a variant ID for our SNPs. Human and model organisms often have annotated SNP names and so plink will look for these. We do not have them so instead we set ours to default to chromosome:position which can be achieved in plink by setting the option @:# - see here for more info.
- --indep-pairwise - finally we are actually on the command that performs our linkage pruning! The first argument, 50 denotes we have set a window of 50 Kb. The second argument, 10 is our window step size - meaning we move 10 bp each time we calculate linkage. Finally, we set an r2 threshold - i.e. the threshold of linkage we are willing to tolerate. Here we prune any variables that show an r2 of greater than 0.1.
- --out Produce the prefix for the output data.

-	Perform the PCA (plink --pca)

```
module load plink/1.90.beta.4.4
plink --vcf SNP.mac3.minQ20.minDP3.miss10.maf0.05.vcfLD0.2.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out PCAPanel2
```

This is very similar to our previous command. What did we do here?

- --extract - this just lets plink know we want to extract only these positions from our VCF - in other words, the analysis will only be conducted on these.
- --make-bed - this is necessary to write out some additional files for another type of population structure analysis - a model based approach with admixture.
- --pca - fairly self explanatory, this tells plink to calculate a principal components analysis.

Once the command is run, we will see a series of new files. We will break these down too:

PCA output:

- PCAPanel2.eigenval - the eigenvalues from our analysis
- PCAPanel2.eigenvec - the eigenvectors from our analysis

plink binary output

- PCAPanel2.bed - the bed file - this is a binary file necessary for admixture analysis. It is essentially the genotypes of the pruned dataset recoded as 1s and 0s.
- PCAPanel2.bim - a map file (i.e. information file) of the variants contained in the bed file.
- PCAPanel2.fam - a map file for the individuals contained in the bed file.

The PCA will be performed in R using the following code:

```
getwd()
# load tidyverse package
library(tidyverse)

# read in data
pca <- read_table2("./PCAPanel2.eigenvec", col_names = FALSE)
eigenval <- scan("./PCAPanel2.eigenval")
str(eigenval)

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
head(pca)
individuals<-pca$ind
individuals[1:20]
write.csv(individuals, "IDorder.csv", row.names = FALSE)

#Include rivers names
ID2<-read.csv("IDorder2.csv")
head(ID2)
IDriverfigure<-read.csv("IDriverfigure.csv")
head(IDriverfigure)
str(IDriverfigure)
str(ID2)
Riverfigure<-merge(ID2,IDriverfigure,by.x = "IDorder", all.x = T)
head(Riverfigure)
str(Riverfigure)
tail(Riverfigure)
head(individuals)
tail(individuals)
head(pca)

#remake pca dataframe with the rivers names
pca2 <- as.tibble(data.frame(Riverfigure,pca))
head(pca2)

# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot 
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

# plot pca
b <- ggplot(pca2, aes(PC1, PC2, col = Location)) + geom_point(size = 3)
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

```

A PCA similar to this will be obtained (in this case, Panel 2, *P. trichocarpa*)

![PCAPanel2locations](galaxy_workflows/Imputation-based.meta-analysis/PCAPanel2location.png "PCAPanel2location")

However, Discriminant Analysis of Principal Components (DAPC) has shown a better performance in clustering individuals in higher in size and complexity genetic datasets of biological populations than Bayesian Methods such as STRUCTURE, which require considerable computation time, and than other more fast and flexible Multivariate analysis such as PCA [Jombart et al. 2010](https://doi.org/10.1186/1471-2156-11-94). DAPC is a wonderful tool for exploring structure of populations based on PCA and DA without making assumptions of panmixia. Thus, this technique provides a robust alternative to Bayesian clustering methods like STRUCTURE (Pritchard et al., 2000) that should not be used for clonal or partially clonal populations. Using this [tutorial](https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html) for DAPC analysis. This is for [.vcf files](https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html)

**DAPC** is performed in R using the following code:

For the transformation of vcf file to genlight object: https://cran.r-project.org/web/packages/vcfR/vignettes/converting_data.html

For the DAPC analysis: https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf 

```
library(adegenet)
library(vcfR)

#The only input file needed is a vcf or vcf.gz file
vcf <- read.vcfR("SNP2_subset.vcf")

#transform vcf to genlight
x <- vcfR2genlight(vcf) 
x

grp <- find.clusters(x, max.n.clust=40)

#DAPC is implemented by the function "dapc", which first transforms the data using PCA, and then performs a Discriminant Analysis on the retained principal components.Two graphs will appear as a result to help decide the number of principal components and the number of discriminant functions to be retained 
dapc1 <- dapc(x, grp$grp)  
dapc1

#plot the groups in a scatterplot and the DA eigenvalues
scatter(dapc1)

#The scatterplot can be customized (dot size)
scatter(dapc1, cell = 0, pch = 18:23, cstar = 0, mstree = TRUE, lwd = 2, lty = 2)

#position of the DA eigenvalues plot
scatter(dapc1, posi.da="bottomleft", bg="white", pch=17:22)

#Change of colors, position of the PCA eigenvalues plot 
myCol <- c("darkblue","purple","green","orange","red","blue")
scatter(dapc1, posi.da="bottomleft", bg="white", legend = TRUE,
        pch=17:22, cstar=0, col=myCol, scree.pca=TRUE,
        posi.pca="bottomright")

#Membership probabilities
#Exploration
class(dapc1$posterior)
dim(dapc1$posterior)
round(head(dapc1$posterior),3)


#Change formatting for storage and be used to correct by population structure in the GWAS (assignedPopframe)
round<-round(dapc1$posterior,3)
assignedPop<-apply(round, 1, which.max)
assignedPopframe<-as.data.frame(assignedPop)
write.csv(assignedPopframe, "assignedPopPanel2subsetDAPC.csv")

```

**2.4.** Perform fastSTRUCTURE analysis

However, sometimes can be more appropriate to perform a fastSTRUCTURE analysis. A STRUCTURE-like approach assumes that markers are not linked and that populations are panmictic (Pritchard et al., 2000). To use model-free methods K-means clustering based on genetic distance or DAPC are more convenient approaches for populations that are clonal or partially clonal.

The .bed, .fam and .bim files resulting from the PCA analysis in plink (which are the outputs together with the .eigenval and .eigenvec files) are the input files needed for the fastSTRUCTURE analysis (and in addition, they contain independent SNPs r2<0.02):

```
 module load fastStructure/1.0
 structure.py -K 4 --input=Panel2/PCAPanel2 --output=Panel2/PCAPanel2out --full --seed=100
```
"K" was set to 4 because there were 4 original populations in our dataset and the PCA showes a population structure among them

Working directory: /labs/Wegrzyn/IreneCobo/Ptrichocarpa/Panel2/Referencepanel/fastSTRUCTURE/fastStructure

Inside the "fastStructure" folder, there is another folder called "Panel2" containing the three input files: PCAPanel2.bed, PCAPanel2.bim, PCAPanel2.fam



<a name="imp"></a>
### 3. Imputation

Start here if you have **low density SNPs** (assays).

-	Imputation methods work by combining a **reference panel** of individuals genotyped at a dense set of polymorphic sites (usually single-nucleotide polymorphisms, or “SNPs”) with a **study sample** collected from a genetically similar population and genotyped at a subset of these sites.
-	Imputation methods predict unobserved genotypes in the study sample by using a population genetic model to extrapolate allelic correlations measured in the reference panel.
-	The imputed genotypes expand the set of SNPs that can be tested for association, and this more comprehensive view of the genetic variation in a study can enhance true association signals and facilitate meta-analysis

To accomplish this imputation, we need a **Study panel**, a **Reference panel** and a **Recombination map**.

**3.1.** Reference panel: high density SNPs (WGS or exome capture). Before using it for imputation, perform the Data quality filtering (section 2). For the study panel(s), the data quality filtering will be performed after imputation.

**3.2.** The Reference panel and the Study panel have to be mapped against the same version of the reference genome. If this is not the case, the Study panel has to be remapped against the same version of the reference genome than the Reference panel using the following steps:

- Get the flanking regions of each SNP in the study panel using the following code:

1. Index reference genome used to obtain the SNPs in your SNP panel and your reference panel. In my case, the authors used v3 (for both the SNP panel and the reference panel)

```
samtools faidx Ptrichocarpa_210_.v3.0.fa
```
Input reference genome in fasta format, unzipped

2. Obtain the sizes of the chromosomes in this reference genome (v3 in my case). 
```
cut -f1,2 Ptrichocarpa_210_.v3.0.fa.fai > sizesv3.genome
```
-f flag: get 1 and 2 columns of the indexed reference genome (first column, chromosome name, 2nd column, chromosome lenght)


3. Get 200 pb flanking regions in .bed format
```
flankBed -i Referencepanel.bed -g sizesv3.genome -b 200 > Referencepanelflank200pb.bed
```
-i flag: .bed file with the SNP positions in your study panel (format, first column chromosome name, second and third column, SNP position, repeated)

```
chr01	52	52
chr01	117	117
chr01	133	133
chr01	2387	2387
chr01	3111	3111
```

-g flag: the output of the previous step (2. step)

-b flag: the number of pair bases you want for the flanking regions (I recommend more than 70 pb to be able to use bwa mem in the mapping steps)

The output will be a .bed file with the flanking regions positions of the SNPs in your study panel

4. Keep only the flanking region on the right

```
awk 'NR%2==0' Referencepanelflank200pb.bed > Referencepanelflank200pbpares.bed
```

5. Transform the flanking regions in .bed format to .fasta to be mapped against the newer version of the reference genome using bwa mem
```
bedtools getfasta -fi Ptrichocarpa_210_v3.0.fa -bed Referencepanelflank200pbpares.bed -fo Referencepanelflank200pbpares.fasta
```

-fi flag: Reference genome from which the SNP positions of your study panel where obtained (in fasta format, unzipped)

-bed flag: Flanking region positions in your study panel in .bed format, output from the previous step (3. step)

-fo flag: The output (flanking regions in fasta format)


- Mapping the flanking regions in fasta format against the same version of the reference genome used in the Reference panel using bwa:

1. Index the reference genome

```
module load bwa/0.7.17
bwa index Ptrichocarpa_533_v4.0.fa
```

2. Mapping

```
bwa mem Ptrichocarpa_533_v4.0.fa Panel4flank200pbpares.fasta>Panel4flank200pbpares.sam
```

3. Code to extract the positions:
Remove header from sam file (header lines start with the at sign @)

```
sed '/@/d' ./Panel4flank200pbpares.sam>Panel4flank200pbparesnoheader.sam
```
In the alignment, 3rd column is the chromosome, 4th column is the position in which the flanking region start mapping against the reference genome. 1st column is the flanking region name and 5th and 6th columns, mapping quality. Get these columns using the following code:

```
awk '{print $1,$3,$4=$4-1,$5,$6}' Panel4flank200pbparesnoheader.sam > Panel4posv4.txt
```

 "-1" was also substracted from the position column (forth column), so that the v4 positions are obtained

4. After bwa mapping, some of the scaffolds can be duplicated, some of them showing different qualities, others just doubled. Remove the doubled ones, keeping just one per couple, and in the case of those with different qualities, if they are few (192 out of 34,000 in my case), just remove them. It is also possible to filter by quality using samtools: 

5. Include the positions in the SNPs file, using "merge" function in R

6. Build a high density recombination map from the Reference panel for imputation, since a fine-scale recombination map for the region to be analyzed is required in Imputev2. This file should have three columns: physical position (in base pairs), recombination rate between current position and next position in map (in cM/Mb), and genetic map position (in cM). The file should also have a header line with an unbroken character string for each column (e.g., "position COMBINED_rate(cM/Mb) Genetic_Map(cM)"). 
