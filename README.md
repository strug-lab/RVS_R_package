RVS_R_package
=============

Robust Variance Score R package

The robust variance score (RVS) is a novel likelihood-based method for genetic association with NGS data from external control groups. RVS substitutes genotype calls by their expected values given observed sequence data and implements a robust variance estimate for the score statistic.

RVS eliminates read depth bias in the estimation of minor allele frequency. It controls Type I error and has comparable power to the 'gold standard' analysis with the true underlying genotypes for both common and rare variants.

The R current package contains core functions but does not process or filter the vcf file yet.
