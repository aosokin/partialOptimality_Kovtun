
This software implements the MATLAB mex-wrapper Kovtun's algorithm to compute 
partially optimal labelings of multilabel pairwise Potts MRFs

Anton Osokin, (firstname.lastname@gmail.com)
11.08.2013

Please cite the following papers in any resulting publication:

Karteek Alahari, Pushmeet Kohli, P. H. S. Torr
Dynamic Hybrid Algorithms for MAP Inference in Discrete MRFs 
In IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI), 2010. 

I. Kovtun, "Partial optimal labeling search for a NP-Hard subclass of (max,+) problems"
DAGM Symposium, pp. 402-409, 2003.

PACKAGE
-----------------------------

./partialOptimalityMex.cpp - the C++ code of the wrapper

./buildPartialOptimalityMex.m - function to build the wrapper

./partialOptimalityMex.m - the description of the implemented function

./example_partialOptimalityMex.m - the simple example of usage

./example_partialOptimalityMex_Tsukuba.m - example of usage on the Tsukuba stereo MRF instance

./tsukubaInstance.mat - Tsukuba instance. Potentials were precomputed by Karteek Alahari: 
http://www.di.ens.fr/~alahari/papers/alahari10.pdf

./projectEnergyPotts.m - a function to project an MRF instance on a smaller one by fixing some variables

./src - C++ code by Karteek Alahari (the code is modified by Anton Osokin)
visit http://www.di.ens.fr/~alahari/code/eff_mlabel_ver1.zip for the original code
Note, that this package include quite a lot of things but the wrapper provides interface only to the Kovtun's method

./partialOptimalityMex.mexw64 - binary file for the mex-function compiled using MATLAB R2013a + MSVC 2012

USING THE CODE
-----------------------------

0) Install MATLAB and one of the supported compilers

1) Run buildGraphCutMex.m

2) Run example_partialOptimalityMex.m and example_partialOptimalityMex_Tsukuba.m to check if the code works

The code was tested using MSVC 2012 and MATLAB 2013a

