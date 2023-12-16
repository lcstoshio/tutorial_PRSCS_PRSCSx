```
Script to run PRS-CSx
Last edited in December 16th, 2023
PGC: Trainers Development Program Brazil 
Lucas Toshio Ito
```

# PRS-CSx

PRS-CSx is a Python based command line tool that integrates GWAS summary statistics and external LD reference panels from multiple populations to improve cross-population polygenic prediction. Posterior SNP effect sizes are inferred under coupled continuous shrinkage (CS) priors across populations. PRS-CSx is an extension of the Bayesian polygenic prediction method PRS-CS (https://github.com/getian107/PRScs), described in the article:

T Ge, CY Chen, Y Ni, YCA Feng, JW Smoller. Polygenic Prediction via Bayesian Regression and Continuous Shrinkage Priors. Nature Communications, 10:1776, 2019.

The development and evaluation of PRS-CSx are described in:

Y Ruan, YF Lin, YCA Feng, CY Chen, M Lam, Z Guo, Stanley Global Asia Initiatives, L He, A Sawa, AR Martin, S Qin, H Huang, T Ge. Improving polygenic prediction in ancestrally diverse populations. Nature Genetics, 54:573-580, 2022.

An application of the "meta" and "auto" version of PRS-CSx is described in:

T Ge et al. Development and validation of a trans-ancestry polygenic risk score for type 2 diabetes in diverse populations. Genome Medicine, 14:70, 2022.

## Getting Started - Already done (check https://github.com/getian107/PRScs for setup details)

- Clone this repository (executable files).

`git clone https://github.com/getian107/PRScs.git`

- Download the LD reference panels:

    LD reference panels constructed using the 1000 Genomes Project phase 3 samples (should use the function `tar -zxvf ldblk_1kg_POP.tar.gz` before using):
    
     African - /home/PGC-TDP/refs/ldblk_1kg_afr   
     American - /home/PGC-TDP/refs/ldblk_1kg_amr  
     East Asian - /home/PGC-TDP/refs/ldblk_1kg_eas  
     European - /home/PGC-TDP/refs/ldblk_1kg_eur  
     South Asian - /home/PGC-TDP/refs/ldblk_1kg_sas  
     
- PRScs requires Python packages **scipy** (https://www.scipy.org/) and **h5py** (https://www.h5py.org/) installed.
  
- Download the SNP information file:

    1000 Genomes reference: SNP info - /home/PGC-TDP/refs/snpinfo_mult_1kg_hm

- Once Python and its dependencies have been installed, and will show up the command-line options.

  PRScsx.py --help
  PRScsx.py -h

## Formatting the GWAS Summary Statistics

The summary statistics file must include the SNP identifier (in rsid format), Allele 1 (effective allele) and 2 (alternative allele), BETA/OR (effect/odds ratio), SE/P (standard error/ pvalue) and follow this set format. Also the columns **must** follow this column order and header line names of this set example, changing to OR or BETA and P or SE based on the value that you are using.

    SNP          A1   A2   OR/BETA   P/SE
    rs4970383    C    A    0.7164   0.4778
    rs4475691    C    T    1.0145   0.1245
    rs13302982   A    G    0.5232   0.2429
    ...

## Okay let's run the *Test Data*

The test data contains GWAS summary statistics for chromosome 22 and plink files.

### Creating work folder for PRS
```
mkdir prscsx
cd prscsx
```

### Copying all of the test data to your new folder (includes summary statistics from GWAS and plink files - .bed.bim.fam)
```
cp /home/PGC-TDP/test_prs/* ./
```

### Running PRS-CS
```
refs=/home/PGC-TDP/refs

PRScsx.py \
	--ref_dir=$refs \
	--bim_prefix=./cc.clean \
	--sst_file=./EUR_sumstats.txt,./EAS_sumstats.txt \
	--n_gwas=200000,100000 \
	--pop=EUR,EAS \
	--chrom=22 \
	--meta=True \
	--out_dir=./ \
	--out_name=PRSCSX_TEST
```
	
### Merging the results of each chromosome
```
cat ./PRSCSX_TEST_pst_eff_a1_b0.5_phiauto_chr* > ./PRSCSX_TEST_pst_eff_a1_b0.5_phiauto_all.txt
```

### Calculating the score in PLINK
```
# Results should be written in ./PRSCSX_TEST_Score.profile
plink --bfile ./cc.clean --score ./PRSCSX_TEST_pst_eff_a1_b0.5_phiauto_all.txt 2 4 6 sum center --out ./PRSCSX_TEST_Score
```

### Principal Components Analysis

```
plink --bfile ./cc.clean --pca 10 --out ./cc.clean_pca10
```

## Interpreting results in R

```{r}
# Packages
# install.packages(c("dplyr", "plyr", "ggplot2", "rcompanion"))
library(dplyr)
library(plyr)
library(ggplot2)
library(rcompanion)

# Importing the data
prs_score <- read.table("/home/PGC-TDP/scores/PRSCSX_TEST_Score.profile", header=T)
colnames(prs_score)[6] <- "PRS"

fam <- read.table("/home/PGC-TDP/scores/cc.clean.fam", header=F)
colnames(fam) <- c("famID", "IID", "fatID", "motID", "sex", "pheno")
fam <- mutate(pheno, pheno=factor(pheno, labels=c('Controle', 'Caso')))

eigenvec <- read.table ("./cc.clean_pca10.eigenvec", header = F)
colnames(eigenvec) <- c("IID", paste0("PC", 1:10))

# Merging PRS + Phenotype and PCs
final <- join_all(list(prs_score, fam, eigenvec), by="IID", type="inner)

# Regression with covariables
fit <- lm (PRS ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=final)
final$PRS_Residuos <- residuals(fit)

# Graph
# Calculating means for each group
mean <- ddply(final, "pheno", summarise, grp.mean=mean(PRS_Residuos))

# Counting cases and controls
table(final$pheno)

ggplot(final, aes(x=PRS_Residuos, group=pheno, fill=pheno))+
     geom_density(alpha=0.5)+
     geom_vline(data=mean, aes(xintercept=grp.mean, colour=pheno), linetype="dashed", size=0.8)+
     scale_fill_brewer(palette="Set1") +
     scale_color_brewer(palette="Set1") +
     theme_classic()+
     labs(x="Polygenic Score", y="Density")

# R2 Nagelkerke
Null <- glm(formula = pheno ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial(link="logit"), data = final)
Full <- glm(formula = pheno ~ PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial(link="logit"), data = final)
nagelkerke(Full_SD, null=Null_SD)$Pseudo.R.squared.for.model.vs.null[3]
```
