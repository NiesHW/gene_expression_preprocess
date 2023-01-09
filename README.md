Gene Expression Data Pre-processing

Because the raw gene expression data consist of missing and repeated gene Entrez ID, dataset pre-processing was performed. Missing and repeated data can lead to poor survival analysis and the incorrect interpretation of predictors like the diagnosis stage [1]. Based on [2,3 the missing gene Entrez IDs were removed, and the gene expression values of the repeated gene Entrez IDs were averaged across all of the samples. Table 1 presents the de-tails of the gene expression data used in this research. 

References
1. Nur, U.; Shack, L.G.; Rachet, B.; Carpenter, J.R.; Coleman, M.P. Modelling relative survival in the presence of incomplete data: A tutorial. Int. J. Epidemiol. 2009, 39, 118–128.
2. Liu, W.; Wang, W.; Tian, G.; Xie, W.; Lei, L.; Liu, J.; Huang, W.; Xu, L.; Li, E. Topologically inferring pathway activity for precise survival outcome prediction: Breast cancer as a case. Mol. Biosyst. 2017, 13, 537–548.
3. Mohammed, A.; Biegert, G.; Adamec, J.; Helikar, T. Identification of potential tissue-specific cancer biomarkers and devel-opment of cancer versus normal genomic classifiers. Oncotarget 2017, 8, 85692–85715, doi:10.18632/oncotarget.21127.
