RVS_R_package
=============

Robust Variance Statistic Score R package

The robust variance score (RVS) is a likelihood-based method/test to compute p-values to analyze and interpret genetic association data. It substitutes genotype calls by their expected values given observed sequence data.

RVS eliminates read depth bias in the estimation of minor allele frequency. It controls Type I error and has comparable power to the 'gold standard' analysis with the true underlying genotypes for both common and rare variants.

The R current package contains all the core functions but does not process or filter the vcf file yet.