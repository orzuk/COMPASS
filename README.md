COMPASS
=======
Convex Optimization for Microbial Profiling by Aggregating Short Sequence reads


(1.) General:

The COMPASS repository contains programs for reconstruction of microbial species identities 
and frequencies from massively parallel sequencing short-reads data. 

The repository includes source code, executables and example datasets, and can be used for 
microbial profiling of experimental and simulated reads, for simulating reads, and for evaluating
prifling results. 


(2.) Installation:
Download the file COMPASS.tgz from https://github.com/orzuk/COMPASS
Extract the COMPASS.tgz file to a directory of your choice, add to Matlab path and run.
Alternatively, you can clone the repository directly from github, by running the following: 

git ... [T.B.D.]


2.1 Requirements: 


2.2. Dependencies:
COMPASS requires the installation of additional packages for some of its operations. 
To ensure operations of all COMPASS modules, you will need to download the following packages: 

(i) CVX - Matlab convex optimization toolbox. 
Please download and install from http://cvxr.com/cvx/.
The package is required for mixture reconstruction module.

(ii) Mothur - Microbial ecology bioinformatics package. 
   Please download from http://www.mothur.org/ , install in a directory of your choice and add to path.  
This package is required only for the solution comparison module of COMPASS. 

2.3. Mex files:
To enhance performance, parts of COMPASS are written in C and compiled to .mex files accessible by Matlab. 
We’ve provided .mex files for linux in the directory /src/mex. 
If you need to create .mex files for your environment (e.g. Mac/Windows), you will need to run 
the following command in Matlab: 
		CompileCOMPASSMexFiles 
or contact shental@openu.ac.il for help in creating a Mac/Windows version. 

	
2.4 Demo/Usage Example: 

2.4.1. After installation, you can run the following script/function example_of_a_single_simulation.m
This script demonstrates most of COMPASS' main modules. 
It simulates a bacterial mixture, sample reads, reconstruct identities, and compare reconstructed solution to original vector. 

In addition, we preovide two more examples: 

2.4.2. Example of simulations similar to those performed to produce Figure 3 in the manuscript [1]. 
Namely varying the number of reads, read length and number of bacteria – See file: "example_simulations.m"

2.4.3.Example of analysing the Drosophila larva sample L2 from [1] using COMPASS – 
See file: " example_experimental_reads_sample_L2.m "  




3. Directories in the COMPASS package: 

The following directories are avialable after installation: 
3.1. mFiles (or src?) – Matlab files and mex files.
3.2. results – A directory used to save simulation results for the examples. Note that we have already saved several example results that appear in the introductory examples.
3.3. database – the Greengenes database used in a COMPASS format. Two databases are provided – the full 16S rRNA gene database and the database of 750bp long sequences covering variable regions V3-V6 (see manuscript)
3.4. experimentalReads –  example files for larva sample L2. Other experimental data presented in the paper are available at the MG-RAST website: http://metagenomics.anl.gov/linkin.cgi?project=5237

In addition, you will need to add the following directories to the COMPASS master directory (see Section 2.2 Dependencies)
3.5. CVX – the optimization software used, downloaded from the CVX website (http://cvxr.com/cvx/) 
3.6. mothur – the mothur software (http://www.mothur.org/) used to calculated distances between correct and reconstructed bacteria.



4. Overview of main modules

4.1. Preparing the database:
This module prepares a sequence database in .mat format from a fasta file containing 
16S rRNA sequences. 
To prepare sequence database for a region of interest, perform the following steps:
(i) download the desired database as a fasta file 
(ii) If a specific PCR-amplified region is needed, perform in-silico PCR for the desired primers, and output the region from all amplified sequences to a new fasta file (we will name it for convenience db.fasta).
(iii) Run in the command line the PrepareDB script:
PrepareDB database-name output-mat-file-name (T.B.D.)

T.B.D.:
repeatWhenLowerThanThisValue
realReadsFromFileFlag
preprocessByExistanceFlag
Provide code for splitting the database


This will result in output-mat-file-name, a .mat file with prepared database in STIRP format, which is then used in the main solver (section 3.2). 

The program comes with the 16S sequence database described in the paper (see section 4). 

4.2. Reconstruct mixture identities from reads:
This is the main utility of COMPASS. It reconstructs a microbial sample from short read data, 
and a prepared sequence database (see Section 3.1).

Main solver program: SolveMixingMatrixFromReads

Input:
Reads (fasta file)
Database (T.B.D.)
Parameters

Output:
T.B.D.

Usage: SolveMixingMatrixFromReads reads-file database-file parameters-file output-file-name



4.3. Experiment processing:

Experimental reads – Reads should be prepared in a specific way for COMPASS. 
For an example of loading a fasta file and preparing the needed input for COMPASS, 
see file – " example_read_fasta.m "

script
Input:
reads (fasta)
output:
T.B.D.
Algorithm:
(i) Local noise processing (k-mer to k-10-mer) → normalized reads (T.B.D.)
(ii) Global noise processing (median filtering) → normalized reads (T.B.D.)
(iii) Main solver
(iv) L1-solution minimization

Usage: PreprocessData

4.4. Read Simulation Module:
This module simulates short-read data obtained from a bacterial mixture. 

Reads can be simulated using the file  "createReads_package.m/SimulateReadsFromSequences.m", 
that appears as part of the introductory examples.

Input:
number of reads, error model, number of bacteria, freq. distribution, read length, bacteria ID (optional)
Output:
simulated reads (fasta? or reads for solver?)

Usage: SimulateReadsFromSequences ….


4.5. Evaluation of Solution Accuracy:
This module evaluates the accuracy of the reconstructed solution, when a ‘ground-truth’ solution 
is available (either in simulations or using other technologies, e.g. Sanger sequencing). 


The file "example_evaluating_simulation_results.m" presents an example of calculating weighted 
specificity and sensitivity for a specific simulation, as used in Figure 3 and Figure 4 in the manuscript [1].


Alternative: 
EvaluateSolutionAccuracy

Input:
Original mixture file (Two columns: ID+Freq)
Reconstructed mixture (Two columns: ID+Freq)
Similarity threshold
Frequency threshold

Output
Weighted recall and precision (or other metrics?) 

Usage: 

EvaluateSolutionAccuracy orig-mixture-file recons-mixture-file sim-metric sim-threshold freq-threshold

5. Data:
The following files are available in the ‘database’ directory: 
current_prokMSA_unaligned.fasta.gz - a fasta file with 16s database - Downloaded from greengenes database  (T.B.D.)
s16_data_uni_packed.mat - a preprocessed version of the 16s database. Can be generated by running the following: (T.B.D.)



6. Acknowledgment

COMPASS was developed by: 
Amnon Amir, Noam Shental, Amit Zeisel and Or Zuk, as part of work on the papers,


[1] “High Resolution Microbial Community Reconstruction by Integrating Short Reads from Multiple 16S rRNA Regions“ 
A. Amir, A. Zeisel, O. Zuk, M. Elgart, S. Stern, O. Shamir, P.J. Turnbaugh, Y. Soen and N. Shental (submitted).

[2] “Accurate Profiling of Microbial Communities from Massively Parallel Sequencing using Convex Optimization”, O. Zuk, A. Amir, A. Zeisel, O. Shamir and N. Shental (SPIRE13).

Please cite the above papers if using the package. 

7. Support
For any questions or comments, please contact: 
- Noam Shental: shental@openu.ac.il
- Or Zuk: orzuk@ttic.edu

