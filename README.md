# Propagation_Analysis_Code
Propagation Analysis Code version 1.0 © GPL-3.0, 2019, ETH Zurich, Institute of Biochemistry, Ulrich Berge 

Description
This set of functions computes out of a matrix with lineage-analysis based single cell data. The principle of this script is illustrated in Figure S7 of the manuscript:
A) An observed division type matrix as described in Figure 3 of the manuscript. This division type matrix is based on the parameter (here cell cycle duration), where cells are classified in outlier cells and not outlier cells (ODM).
B) The empirical probability to obtain A) based on a probability function estimated by permuting the cell IDs
C) Estimation of the probability distribution function of obtaining any possible occurrence value per division type
D) Power of the test visualized as a histogram
Of note:
•	Although in our study the analysed parameter is cell cycle duration, also other parameters obtained as read-outs from biological reports could be studied. 
•	Due to the nature of a randomization reproducibility is limited but minor.
•	This is a script, therefore no source code, system requirements and installation guidelines are given.
•	There are no dependencies.
•	The runtime of the script depends on the input file size and the number of requested permutations

Computation
The following MATLAB functions run on MATLAB R2016b and requires the statistics toolbox.
All MATLAB functions need to be in the same directory. 
For each function, the input and output is defined and explained at the beginning of the code.
1.	Read in your data (or the provided TestDataset.dat)
2.	Run PropagationAnalysis.m; this will, given a provided number of permutations, the estimated probabilities as described in Figure 3 of the manuscript.
a.	This function will call the function ODM.m (for its principle see Figure S5 in the manuscript)
b.	This function will call the function makeAbsolTransMatrix.m to generate the empirical division type matrix based on the classification in step 2a (above) of the computation
3.	Run Plotting.m; this will plot the data obtained in step 2. (above) of the computation to generate figure parts as in Figure 3.


Example
We provide the Source Data as a separate Zip-file with which the results of Figure 3 and 5 of the manuscript can be, due to the randomization with minor variations, reproduced.
The computation of this data set with 10,000 permutations (and 2 permutation rounds takes in total about 5-10 h, with about 1 to 3 seconds for each permutation.


Reference
The manuscript/accepted publication of Berge et al. «Asymmetric division events promote variability in cell cycle duration in animal cells and E.coli «






