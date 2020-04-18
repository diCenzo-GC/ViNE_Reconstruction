# ViNE: A Virtual Nodule Environment

ViNE is a metabolic model of a nodulated legume, consisting of the interacting partners *Medicago truncatula* and *Sinorhizobium meliloti*. The final model contains root and shoot tissues, as well as nodule tissue subdivided into five distinct developmental zones. The nodule tissue contains representations of the metabolism of both the plant and the bacterial partners, as well as cross-talk between the organisms.

Within this directory you can find all of the scripts used to build ViNE, as well as the scripts used to build an updated metabolic model of free-living *S. meliloti*. In addition, all scripts for the analysis of ViNE in our associated manuscript are provided.

## Citations

Construction of ViNE and our analyses of ViNE, as well as the construction of the updated *S. meliloti* model, are described in the following preprint:\
diCenzo GC, Tesi M, Pfau T, Mengoni A, Fondi M (2019) [A Virtual Nodule Environment (ViNE) for modelling the inter-kingdom metabolic integration during symbiotic nitrogen fixation](https://www.biorxiv.org/content/10.1101/765271v1). bioRxiv. doi:10.1101/765271

The  *M. truncatula* model used as a starting point in this work was described in the following publication:\
Pfau T, Christian M, Masakapalli SK, Sweetlove, LJ, Poolman MG, Ebenhöh O (2018) [The intertwined metabolism during symbiotic nitrogen fixation elucidated by metabolic modelling](https://www.nature.com/articles/s41598-018-30884-x). Scientific Reports. 8:12504.

The RNA-seq data used to constrain the nodule reaction space was described in the following publication:\
Roux B, Rodde H, Jardinaud M-F, Timmers T, Sauviac L, Cottret L, Carrère S, Sallet E, Courcelle E, Moreau S, Debellé F, Capela D, de Carvalho-Niebel F, Gouzy J, Bruand C, Gamas P (2014) [An integrated analysis of plant and bacterial gene expression in symbiotic root nodules using laser‐capture microdissection coupled to RNA sequencing](https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.12442). The Plant Journal. 77:817-837.

## Required Inputs

### Dependencies

If you wish to run our scripts to build ViNE and/or to analyze ViNE, the following dependencies are required. In all cases, the versions that we use are listed.
* [MATLAB](https://www.mathworks.com/products/matlab.html) R2016b
* [libSBML](https://www.sourceforge.net/projects/sbml/files/libsbml) version 5.13.0
* [SBMLToolbox](https://www.sourceforge.net/projects/sbml/files/SBMLToolbox) version 4.1.0
* [COBRA Toolbox](https://opencobra.github.io/cobratoolbox/stable/) commit 9b10fa1
* [Tn-Core Toolbox](https://www.github.com/diCenzo-GC/Tn-Core) version 2.2
* [TIGER Toolbox](https://csbl.bitbucket.io/tiger/download.html) version 1.2-beta
* [FASTCORE](https://www.uni.lu/forschung/fstc/life_sciences_research_unit/research_areas/systems_biology/software/fastcore) version 1.0
* [iLOG CPLEX Studio](https://www.ibm.com/products/ilog-cplex-optimization-studio) version 12.7.1

### Input Files

If you wish to run our scripts to build ViNE, a few additional files need to be downloaded and placed in the correct directory.
* The previously published version of the *M. truncatula* metabolic reconstruction needs to be [downloaded](https://github.com/sysbiolux/MedicagoScripts/blob/master/Matlab/Data/MedicagoTruncatula.xml) and placed in the directory 2_ViNE_Reconstruction/4_Integrate_Models/1_Produce_Plant/Data as well as 2_ViNE_Reconstruction/4_Integrate_Models/2_Produce_Nodule
* The *chem_xref.tsv* and *reac_xref.tsv* files must be downloaded from [MetaNetX](https://www.metanetx.org/mnxdoc/mnxref.html), placed in the directory 2_ViNE_Reconstruction/4_Integrate_Models/5_Finalize_Nodule, and the files renamed *chem_xref.txt* and *reac_xref.txt*, respectively.

If you wish to run our scripts to analyze ViNE, a couple additional files need to be downloaded and placed in the correct directory.
* The *chem_xref.tsv* and *reac_xref.tsv* files must be downloaded from [MetaNetX](https://www.metanetx.org/mnxdoc/mnxref.html), placed in the directory 3_ViNE_Analyses/4_Carbon_Source/Add_Sucrose_Metabolism/Finalize_Nodule, and the files renamed *chem_xref.txt* and *reac_xref.txt*, respectively.

## Contacts

If you have questions about ViNE, please send inquiries to George diCenzo (george.dicenzo [at] queensu.ca) or Marco Fondi (marco.fondi [at] unifi.it).
