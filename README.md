# OPM
Combining phenotypic and genomic data to improve prediction of binary traits

Each folder of AENOPM, ALASSOOPM, RAWOPM, and SCADOPM contain the R code for the respective case. We have added OPMpackage folder with R package for easier use.

To install the package use the following R command

devtools::install_github("royarkaprava/OPM/OPMpackage")

Open the Usage.R file. The first two lines may also be used to install and load the package. After that you need to download the simulated data.rda or the real data from "OPM\Realdata" folder separately and edit the line load("path to the data") in Usage.R accordingly. Alternatively one may "clone" this github folder and use the Usage.R directly. The coded paths will work fine then.

The data is from 

Liang, Zhikai, Yumou Qiu, and James C. Schnable. "Genome–Phenome Wide Association in Maize and Arabidopsis Identifies a Common Molecular and Evolutionary Signature." Molecular plant 13.6 (2020): 907-922.
