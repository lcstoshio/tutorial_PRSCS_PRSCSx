# PRS-CS

**PRS-CS** is a Python based command line tool that infers posterior SNP effect sizes under continuous shrinkage (CS) priors
using GWAS summary statistics and an external LD reference panel. Details of the method are described in the article:

T Ge, CY Chen, Y Ni, YCA Feng, JW Smoller. Polygenic Prediction via Bayesian Regression and Continuous Shrinkage Priors. *Nature Communications*, 10:1776, 2019.

## Getting Started - Already done (check https://github.com/getian107/PRScs for setup details)

- Clone this repository (executable files).

`git clone https://github.com/getian107/PRScs.git`

- Download the LD reference panels:

    LD reference panels constructed using the 1000 Genomes Project phase 3 samples:
    
     African - /home/PGC-TDP/refs/ldblk_1kg_afr.tar.gz     
     American - /home/PGC-TDP/refs/ldblk_1kg_amr.tar.gz        
     East Asian - /home/PGC-TDP/refs/ldblk_1kg_eas.tar.gz        
     European - /home/PGC-TDP/refs/ldblk_1kg_eur.tar.gz     
     South Asian - /home/PGC-TDP/refs/ldblk_1kg_sas.tar.gz
     
- PRScs requires Python packages **scipy** (https://www.scipy.org/) and **h5py** (https://www.h5py.org/) installed.
  
- Download the SNP information file:

    1000 Genomes reference: SNP info - /home/PGC-TDP/refs/snpinfo_mult_1kg_hm

- Once Python and its dependencies have been installed, and will show up the command-line options.

  PRScs.py --help
  PRScs.py -h

## Formatting the GWAS Summary Statistics

The summary statistics file must include the SNP identifier (in rsid format), Allele 1 (effective allele) and 2 (alternative allele), BETA/OR (effect/odds ratio), SE/P (standard error/ pvalue) and follow this set format. Also the columns **must** follow this column order and header line names of this set example, changing to OR or BETA and P or SE based on the value that you are using.

    SNP          A1   A2   OR/BETA   P/SE
    rs4970383    C    A    0.7164   0.4778
    rs4475691    C    T    1.0145   0.1245
    rs13302982   A    G    0.5232   0.2429
    ...

## Okay let's run the *Test Data*

The test data contains GWAS summary statistics and a bim file for 1,000 SNPs on chromosome 22.

