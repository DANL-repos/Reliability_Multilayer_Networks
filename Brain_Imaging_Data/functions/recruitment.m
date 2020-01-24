function R = recruitment(MA,systemByNode)
%RECRUITMENT      Recruitment coefficient
%
%   R = RECRUITMENT(MA,systemByNode) calculates the recruitment coefficient
%   for each node of the network. The recruitment coefficient of a region
%   corresponds to the average probability that this region is in the same
%   network community as other regions from its own system.
%
%   Inputs:     MA,     Module Allegiance matrix, where element (i,j) 
%                       represents the probability that nodes i and j
%                       belong to the same community
%               systemByNode,	vector or cell array containing the system
%                       assignment for each node
%
%   Outputs:    I,              recruitment coefficient for each node
%   _______________________________________________
%   Marcelo G Mattar (08/21/2014) 


% Initialize output
R = zeros(length(systemByNode),1);

% Make sure the diagonal of the module allegiance is all nan
MA(logical(eye(size(MA)))) = nan;

% Calculate the recruitment for each node
if iscell(systemByNode)
    for i=1:length(systemByNode)
        thisSystem = systemByNode{i};
        R(i) = nanmean(MA(i,strcmp(systemByNode,thisSystem)));
    end
else
    for i=1:length(systemByNode)
        thisSystem = systemByNode(i);
        R(i) = nanmean(MA(i,systemByNode==thisSystem));
    end
end