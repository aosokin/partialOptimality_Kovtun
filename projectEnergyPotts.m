function [ projectedUnary, projectedPairwise, projectionInfo ] = projectEnergyPotts( dataCost, neighbors, labels )
% projectEnergyPotts projects the energy onto the partial labelling, i.e. fixes some labels
%
% [energy, dataTerm, smoothnessTerm] = projectEnergyPotts(dataCost, neighbors, labels)
%
% INPUT
%   dataCost   - unary potentials ( double[ numLabels x numNodes ])
%   neighbors  - paiwise Potts potentials ( sparse double[ numNodes x numNodes ]). The function uses only upper triangle of this matrix.
%   labels - the partial labeling ( double[ numNodes x 1 ], 1-based indexing for labels, "0" means unlabeled)
%
% OUTPUT
%   projectedUnary  - unary potentials of the new energy  ( double[ numLabels x numNodesNew] )
%   projectedPairwise - pairwise Potts potentials of the new energy ( sparse double[ numNodesNew x numNodesNew ] )
%   projectionInfo - structure with information, allowing to restore the original labeling
%      fields:
%           initialLabels - labels of the initials energy
%           initialNodeId - ids of the new nodes in the initial energy
%           energyConstant - constant terms from the initial energy (after fixing some labels)
%
% Anton Osokin (firstname.lastname@gmail.com),  14.05.2013

if ~isnumeric(dataCost) || ~ismatrix(dataCost)
    error('projectEnergyPotts:badDataCost', 'dataCost should be a matrix numLabels x numNodes');
end
dataCost = double(dataCost);
numLabels = size(dataCost, 1);
numNodes = size(dataCost, 2);

if ~isnumeric(neighbors) || ~ismatrix(neighbors) || ~issparse(neighbors) || size(neighbors, 1) ~= numNodes || size(neighbors, 2) ~= numNodes
    error('projectEnergyPotts:badNeighbors', 'neighbors should be a sparse matrix numNodes x numNodes');
end

if ~isnumeric(labels) || ~iscolumn(labels) || length(labels) ~= numNodes || any(~ismember( unique(labels), [0, 1 : numLabels] ))
    error('projectEnergyPotts:badLabels', 'labels should be a column vector of length numNodes, 1-indexing for label IDs, 0 - unlabeled');
end

% find new node IDs
labeledMask = (labels ~= 0);
nodeIdNew = find( ~labeledMask );

%% projecting unary potentials
projectedUnary = dataCost;

% update unary potentials using their labeled neighbors
[neighRow, neighCol, neighWeight] = find( neighbors );

% leave only upper triangle to avoid edge duplication
toDeleteMask = neighCol <= neighRow;
neighRow( toDeleteMask ) = [];
neighCol( toDeleteMask ) = [];
neighWeight( toDeleteMask ) = [];

% update new unary potentials using old pairwise
id = ones(numLabels, numLabels) - eye( numLabels );

% first direction
id1 = [neighRow; neighCol];
id2 = [neighCol; neighRow];
edgeWeight = [neighWeight; neighWeight];

toUpdate = labels(id1) == 0 & labels(id2) ~= 0;
updateId = id1(toUpdate);
updateVector = bsxfun(@times, id( :, labels( id2(toUpdate) ) ), edgeWeight( toUpdate )' );
for iLabel = 1 : numLabels
    update = accumarray(updateId, updateVector(iLabel, :)', [numNodes, 1])';
    projectedUnary(iLabel, :) = projectedUnary(iLabel, :) + update;
end

% do the actual projecting
projectedUnary = projectedUnary(:, ~labeledMask);

%% projecting pairwise potentials
projectedPairwise = neighbors( ~labeledMask, ~labeledMask);

% set projection info
projectionInfo = struct;
projectionInfo.initialLabels = labels;
projectionInfo.initialNodeId = nodeIdNew;

unaryConstant = sum( dataCost( labels(labeledMask) + numLabels * (find(labeledMask) - 1) ) );
edgesLabeled = labels(neighRow) ~= 0 & labels(neighCol) ~= 0;
pairwiseConstant = sum( (labels(neighRow(edgesLabeled)) ~= labels(neighCol(edgesLabeled))) .* neighWeight(edgesLabeled) );
projectionInfo.energyConstant = unaryConstant + pairwiseConstant;

end

