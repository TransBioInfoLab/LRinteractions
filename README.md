## LRinteractions: An R Program for Identifying Ligand-Receptor Gene Interactions for scRNA-seq Using Random Forest Based Methods

Min Lu, Yifan Sha, Tiago C. Silva1, Antonio Colaprico, Xiaodian Sun, Yuguang Ban, Lily Wang, Brian D. Lehmann and X. Steven Chen
<br>
Maintainer: Min Lu \<m.lu6@umiami.edu\>

**Cite as**<br>
Lu M., Sha Y., Silva T.C., Colaprico A., Sun X., Ban Y., Wang L., Lehmann B.D. and Chen X.S. (2022) *LRinteractions: an R program for identifying ligand-receptor gene interactions for scRNA-se using Random Forest based methods*, available at https://github.com/TransBioInfoLab/LRinteractions


There are two steps to implement our LR hunting algorithm for detecting gene-gene interactions among DIFFERENT cell types in any scRNA-seq data with an illustration example to reproduce the results of the TNBC data analysis in Lu et. al (2021). We have two .R files match these two steps, named as “1_Multiple_Imputation.R” with dataset “LogNormalized_9types.data”, and “2_Ensemble_Results_Pvalues.R” with dataset “LRdb.RData”. The TNBC dataset that contains 7 patients of 9 cell types including "Myeloid", "Epithelial_Basal", "Endothelial", "Myoepithelial", "Epithelial_Luminal_Mature", "iCAFs", "myCAFs", "CD4+ T-cells", "T-Regs" are stored in “LogNormalized_9types.data”.

### Reference<br>
Lu M., Sha Y., Silva T.C., Colaprico A., Sun X., Ban Y., Wang L., Lehmann B.D. and Chen X.S. (2021). *LR Hunting: A Random Forest Based Cell–Cell Interaction Discovery Method for Single-Cell Gene Expression Data*. Frontiers in Genetics. 2021 Aug 20;12:708835.
