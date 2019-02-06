%Release v1_1 - 06/02/2019
% Author: Andrea Iannelli (University of Bristol)
% Supervisor: Dr. Andrés Marcos (University of Bristol)
% contact: andrea.iannelli@bristol.ac.uk (alternatively, andrea.iannelli.uni@gmail.com)
%------------------------------------------------------------------

The MATLAB files in this folder present an implementation of the algorithms
developed in the paper "Robust estimations of the region of attraction using invariant sets", accepted to 
the Journal of The Franklin Institute. The pre-print of the paper is in the folder under the name "ROA_Invariant_prePrint".

  
Files structure:
There are 2 Main files. Specifically:
-JFI_ERA.m allows to perform nominal ROA analysis for two case studies: Van der Pol oscillator and short period aircraft model
-JFI_rERA.m allows to perform robust ROA analysis for two case studies: Van der Pol oscillator and short period aircraft model

Inside, these two files are similarly structured (within the inherent differences due to the different problem solved).
With reference to the nomenclature used in the paper:

%%% NOMINAL
-JFI_ERA_Alg1.m contains the implementation of (nominal) Algorithm 1. Specifically, inside are called
**JFI_ERA_StepA1_1.m-> solving Step A1-1  
**JFI_ERA_StepA1_2.m-> solving Step A1-2
-JFI_ERA_Alg2.m contains the implementation of (nominal) Algorithm 2. Specifically, inside are called
**JFI_ERA_StepA2_1.m-> solving Step A2-1  
**JFI_ERA_StepA2_2.m-> solving Step A2-2
**JFI_ERA_StepA2_3.m-> solving Step A2-3
Note that the hybrid algorithm (Algorithm 3) can be adopted by defining the corresponding option in JFI_ERA.m

%%% ROBUST
-JFI_rERA_Alg1.m contains the implementation of the extension to the uncertainty case of 
Algorithm 1 (not specifically discussed in the paper, but driven by the same rationale of Algorithm 6). Specifically, inside are called
**JFI_rERA_StepA1_1.m-> solving (extension of) Step A1-1  
**JFI_rERA_StepA1_2.m-> solving (extension of) Step A1-2
-JFI_rERA_Alg2.m contains the implementation of the extension to the uncertainty case of 
Algorithm 2 (this is Algorithm 6 in the paper). Specifically, inside are called
**JFI_rERA_StepA2_1.m-> solving (extension of) Step A2-1, i.e. Step A6-1 
**JFI_rERA_StepA2_2.m-> solving (extension of) Step A2-2, i.e. Step A6-2
**JFI_rERA_StepA2_3.m-> solving (extension of) Step A2-3, i.e. Step A6-3
Again, the hybrid algorithm can be adopted by defining the corresponding option in JFI_rERA.m


Some results obtained with the provided routines are also provided as examples:
-VdP_ERA_JFI_Alg3_deg4.mat: Nominal Van der Pol case study, quartic level set functions, Hybrid algorithm
-VdP_rERA_JFI_Alg3_deg4.mat: Uncertain Van der Pol case study, quartic level set functions, Hybrid algorithm 
-SP_ERA_JFI_Alg3_deg2.mat: Nominal short period case study, quadratic level set functions, Hybrid algorithm
-SP_ERA_JFI_Alg3_deg4.mat: Nominal short period case study, quartic level set functions, Hybrid algorithm
-SP_rERA_JFI_Alg3_deg2.mat: Uncertain short period case study, quadratic level set functions, Hybrid algorithm


NOTE: the suite of libraries "SOSAnalysis" from 
http://www.aem.umn.edu/~AerospaceControl/
must be loaded in order to perform the analyses.
e.g.  addpath(genpath('...\Software\SOSAnalysis'))
Specifically: 
- the library "multipoly" is used for manipulation of polynomials;
- the library "sosopt" is used for SOS optimization

