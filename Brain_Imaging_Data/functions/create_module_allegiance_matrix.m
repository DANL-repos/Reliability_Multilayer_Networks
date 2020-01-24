function Algnc = create_module_allegiance_matrix(S,figureFlag)
% PURPOSE:	This function generates the allegiance matrix from a community
%           assignment vector. Function will collapse data across time and
%           repetitions if applicable.
%
% INPUT:
%           S:	Community assignment matrix (function accepts two formats)
%               Format 1: CELL array
%               * Cells | repetitions
%               * Rows  | community assignment
%               * Cols  | time windows
%
%               Format 2: MAT (2-D array)
%               * Rows | community assignment
%               * Cols | time windows or repetitions
%
% OUTPUT:
%           Alij: allegiance matrix
%
%--------------------------------------------------------------------------
%
%    Author: Qawi Telesford
%      Date: 2018-01-13
%   Version: 2.0
% Update(s):
%            - Added error checking
%            - Changed function inputs to allow CELL or MAT input
%            - Displaying figure now optional
%            - Initial release only worked on static networks, can now
%              generate allegiance matrix for multilayer data
%            - node2mat is now subfunction
%
%    Date: 2015-07-03
% Version: 1.0 (Initial Release)
%
% Revised by Ryan Lim on 2019-04-18 to optimize the computation speed
%
%--------------------------------------------------------------------------
%% Error Checking
if(nargin<1 || isempty(S))
    error('Please input community structure MAT or CELL array');
end

if(nargin<2 || isempty(figureFlag))
    figureFlag = 0;
end

if(isscalar(S) || ischar(S) || isstring(S) || isstruct(S))
    % Check if data structure is scalar, character, string or struct array
    error('Please input CELL or MATT array');
elseif(~iscell(S))
    % Check is data structure is cell, then if data is vector or data array
    cellFlag = 0;
    if(~ismatrix(S))
        error('Please input CELL or MATT array');
    end
elseif(iscell(S))
    % If data structure is cell, check if all entries vector or data array
    cellFlag = 1;
    for rr = 1:length(S(:))
        if(isscalar(S{rr}) || ischar(S{rr}) || isstring(S{rr}) || isstruct(S{rr}))
            error('CELL entry contains invalid data type (all entries must vector or data array).');
        elseif(~ismatrix(S{rr}))
            error('CELL entry contains invalid data type (all entries must vector or data array).');
        end
    end
end

%% Set paramaters for analysis
% Based on error checking, data is either cell or data array. If data is
% cell array, function determines number of nodes, number of windows and
% number of repetitions
% NOTE: Organization of data is user dependent. Although data dimensions
% are labeled nodes, windows, and repetitions (N�T�R), these dimensions can
% also represent subect-level data (i.e., nodes, repetitions, subjects).
if(cellFlag == 1)
    repetitions = length(S(:)); % number of repetitions
    N = size(S{1},1); % number of nodes
    T = size(S{1},2); % number of windows
    
    for rr = 1:repetitions
        checkN = size(S{1},1);
        checkT = size(S{1},2);
        if(checkN ~= N || checkT ~= T)
            error('Matrices in cell array not all the same size');
        end
    end
    clear('checkN','checkT');
else
    repetitions = 1;
    N = size(S,1); % number of nodes
    T = size(S,2); % number of windows
end
%% Contstruct allegiance matrix
Algnc = zeros(N);
for rr = 1:repetitions
    for tt = 1:T
        % Generate vector of community assignments
        if(cellFlag == 1)
            commCol = S{rr}(:,tt);
        else
            commCol = S(:,tt);
        end
        
        % Find unique community assignments
        commNum = unique(commCol);
        
        for ii = 1:length(commNum)
            findComm = find(commCol==commNum(ii));
        
           
            nodePairs = zeros(length(findComm)*length(findComm),2);
            last_index = 1;
            for xx = 1:length(findComm)
                for yy = 1:length(findComm)
                    if(findComm(xx) ~= findComm(yy) && findComm(xx) < findComm(yy))
                        nodePairs(last_index, 1) = findComm(xx);
                        nodePairs(last_index, 2) = findComm(yy);
                        last_index = last_index + 1;
                        %nodePairs = [nodePairs; findComm(xx) findComm(yy)]; %#ok<AGROW>
                    end
                end
            end
            nodePairs = nodePairs(1:last_index-1,:);
            
            if(~isempty(nodePairs))
                algnc_temp = node2mat(nodePairs,N);
                Algnc = Algnc + algnc_temp;
                clear('algnc_temp');
            end
        end
    end
end
Algnc = Algnc./(repetitions*T);

%% Display allegiance matrix (flag dependent)
if(figureFlag == 1)
    figure; imagesc(Algnc,[0 1]); axis('square'); colorbar;
end

%% Clear extraneous data
clear('ii','rr','tt','xx','yy'); % loop data
clear('N','T','repetitions'); % input data parameters
clear('maxComm','findComm','commCol','commNum','nodePairs'); % community structure data
clear('cellFlag','figureFlag'); % function flags
end

%% Subfunction: node2mat
% This function generates an adjaceny matrix from list of node pairs
function Aij = node2mat(node_list,N)
    Aij_nodes = false(N);

    for ii = 1:size(node_list,1)
        Aij_nodes(node_list(ii,1),node_list(ii,2)) = 1;
        Aij_nodes(node_list(ii,2),node_list(ii,1)) = 1;
    end
    
    Aij = Aij_nodes;
end
