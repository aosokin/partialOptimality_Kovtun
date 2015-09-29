%   partialOptimalityMex computes partially optimal solution for multilabel pairwise MRF with Potts pairwise potentials
%  
%   A wrapper to code by K. Alahari created for the following paper:
%   Karteek Alahari, Pushmeet Kohli, P. H. S. Torr
%   Dynamic Hybrid Algorithms for MAP Inference in Discrete MRFs 
%   In IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI), 2010. 
%  
%   The method is proposed in
%   I. Kovtun, "Partial optimal labeling search for a NP-Hard subclass of (max,+) problems"
%   DAGM Symposium, pp. 402-409, 2003.
% 
%   Usage example:
%   S = partialOptimalityMex(U, P)
% 
%   INPUT:
%       U		- unary terms (double[numLabels, numNodes])
%       P		- matrix of Potts edge coefficients (sparse double[numNodes, numNodes]);
% 
%   OUTPUT: 
%       S		- labeling that has energy E, vector numNodes * 1 of type double (indices are in [0,...,numLabels]); 0 means unlabeled
%
%   Use buildPartialOptimalityMex.m to build MEX-file, run example.m for an example
%
%   Anton Osokin (firstname.lastname@gmail.com), 26.04.2013
