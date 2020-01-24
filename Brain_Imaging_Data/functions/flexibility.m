function F = flexibility(S, nettype)
%FLEXIBILITY    Flexibility coefficient
%
%   F = FLEXIBILITY(S, NETTYPE) calculates the flexibility coefficient of
%   S. The flexibility of each node corresponds to the number of times that
%   it changes module allegiance, normalized by the total possible number
%   of changes. In temporal networks, we consider changes possible only
%   between adjacent time points. In multislice/categorical networks,
%   module allegiance changes are possible between any pairs of slices.
%
%   Inputs:     S,      pxn matrix of community assignments where p is the
%                       number of slices/layers and n the number of nodes
%
%               nettype,   string specifying the type of the network:
%                   'temp'  temporal network (default)
%                   'cat'   categorical network
%
%   Outputs:    F,      Flexibility coefficient
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: PROMISCUITY
%   _______________________________________________
%   Marcelo G Mattar (08/21/2014) 

%% CHECK INPUTS
if (nargin < 2)
	nettype = 'temp';
end

if ~(strcmp(nettype,'temp') || strcmp(nettype,'cat'))
    error('Expected input ''nettype'' to match ''temp'' or ''cat''');
end

if length(size(S)) > 2
    error('S must be a pxn matrix');
end

[numSlices, numNodes] = size(S);


%% CALCULATE FLEXIBILITY
% Pre-allocate totalChanges
totalChanges = zeros(numNodes,1);

if strcmp(nettype,'temp')
    possibleChanges = numSlices-1; % only consider adjacent slices
    for t=2:numSlices
        totalChanges = totalChanges + (S(t,:) ~= S(t-1,:))';
    end
elseif strcmp(nettype,'cat')
    possibleChanges = numSlices*(numSlices-1); % consider all pairs of slices
    for s=1:numSlices
        otherSlices = 1:numSlices;
        otherSlices = otherSlices(otherSlices~=s); % all slices but the current one
        totalChanges = totalChanges + sum(repmat(S(s,:),(numSlices-1),1) ~= S(otherSlices,:),1)';
    end
end

% Calculate output
F = totalChanges/possibleChanges;
