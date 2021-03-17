# OPM
Software to implement the computing described in:

Combining phenotypic and genomic data to improve prediction of binary traits
By Jarquin, Roy, Clarke, and Ghoshal.

Each folder AENOPM, ALASSOOPM, SCADOPM, and RAWOPM, contains the R code for the cases adaptive elastic net with the one-pass method, adaptive LASSO with the one-pass method, smoothly clipped absolute deviation with the one-pass method, and the raw data with the one-pass method. We have also added the OPMpackage folder for easier use.

To install the package use the following R command:

install.packages(devtools)

devtools::install_github("royarkaprava/OPM/OPMpackage")

Next, open the Usage.R file. (The first three lines of this script may also be used to install and load OPMpackage.) After that, you must download the simulated from "OPM\Simulated data" or the real data from "OPM\Real data" folder separately and edit the line load("path to the data") in Usage.R accordingly. Alternatively, one may "clone" this github folder and use the
Usage.R directly. The coded paths should work fine then.

The real data is from:

Liang, Zhikai, Yumou Qiu, and James C. Schnable. &quot;Genome–Phenome Wide
Association in Maize and Arabidopsis Identifies a Common Molecular and Evolutionary
Signature.&quot; Molecular plant 13.6 (2020): 907-922.
