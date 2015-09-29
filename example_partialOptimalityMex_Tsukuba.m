% This example runs partialOptimalityMex on a Tsukuba stereo instance of MRF
% all the potentials were precomputed by Karteek Alahari:
% http://www.di.ens.fr/~alahari/data/pami10data.tgz

% load data
load('tsukubaInstance.mat');

% run Kovtun's method
labels = partialOptimalityMex(tsukubaInstance.unary, tsukubaInstance.pairwisePotts);
fprintf('Fraction of unlabeled nodes: %f\n', sum(labels == 0) / length(labels));

% project the energy (fix labeled nodes)
[ projectedUnary, projectedPairwise, projectionInfo ] = projectEnergyPotts( tsukubaInstance.unary, tsukubaInstance.pairwisePotts, labels );

% show the results: black corrsponds to unlabeled
imshow(reshape(labels, [tsukubaInstance.nRows, tsukubaInstance.nCols]), []);
