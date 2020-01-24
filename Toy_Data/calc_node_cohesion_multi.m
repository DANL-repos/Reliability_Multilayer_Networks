function [Cij,commFlexible,nodeInfo_multi,commInfo_multi,CijInfo] = calc_node_cohesion_multi(community_structure)
% PURPOSE: Calculates node cohesion from the community structure cell array
%          taken from multiple simulations
%
% INPUT:
%	community_structure: Cell array from multiple simulations of
%                        multilayer community structure
%
% OUTPUT:
%
%	* Cij: Mean cohesion matrix of all simulations; similar to adjacency
%          matrix but edge weights represent the number of times two nodes
%          change communities mutually.
%
%	* commFlexible: Space-time diagram of mean flexible community changes
%
%	* nodeInfo_multi: STRUCT variable containing dynamic measure statistics
%       + global: Mean, standard deviation and variance for entire network
%       + node: Mean, standard deviation and variance for individual nodes
%       + all: Dynamic measure output from all simulations
%           - cohesion: Node cohesiveness
%           - disjoint: disjointedness
%           - flexibility: flexibility
%           - strength_cohesion: cohesion strength
%
%	* commInfo_multi: STRUCT variable containing space-time diagrams
%       + commIndex: Array showing unique community shifts.
%                Column 1 - Indicates initial community of node
%                Column 2 - Indicates new community of node
%                Column 3 - Value for community shift in commChange
%       + commChanges:  Space-time array of community changes
%       + commCohesion: Space-time array of mutual community changes
%       + commDisjoint: Space-time array of independent community changes
%       + commFlexible: Space-time array of flexible community changes
%
%   * CijInfo: STRUCT variable containing cohesion matrix statistocs
%       + all: Cohesion matrix from each simulation
%       + mean: Mean of cohesion matrix edges, equal to main output Cij
%       + std: Standard devation of cohesion matrix edges
%       + var: Variance of cohesions matrix edges
%
%--------------------------------------------------------------------------
%
%  Author: Qawi Telesford
%    Date: 2019-04-17
% Version: 2.1
%
% History:  
%           2.1 (2019-04-17) | QKT
%           * Streamlined output to simplify as STRUCT variables
%
%           2.0 (2018-08-03) | QKT
%           * Minor update to naming convention for CDFS
%           * Updated function call output for 'calc_node_cohesion'
%           * Added commFlexibleMean to output (average flexibile changes)
%           * Added Cij_all (to see Cij for individual runs)
%           * Added calculation of node variance
%
%           1.2 (2015-06-14) | QKT
%           * Updated to include cohesiveness, disjointedness, flexibility,
%             and strength calculation (CDFS)
%
%           1.1 (2015-01-16) | QKT
%           * Minor bug fixes
%
%           1.0 (2015-01-16) | QKT
%           * Initial release
%   
%--------------------------------------------------------------------------
%% Error Checking
if(nargin < 1 || isempty(community_structure))
    error('Missing input, please enter cell array from multpile simulations of community structure.');
end

if(~iscell(community_structure))
    error('Input does not appear to be cell array. Please enter cell array from multpile simulations of community structure. If input is from single simulation, please use ''calc_node_cohesion'' function instead');
end

%% Initialize Variables
N = size(community_structure{1},1);
dT = size(community_structure{1},2) - 1;

repetitions = size(community_structure,1);
Cij_sum = zeros(N);

cohesion_node_all          = zeros(N,repetitions);
disjoint_node_all          = zeros(N,repetitions);
flexibility_node_all       = zeros(N,repetitions);
strength_cohesion_node_all = zeros(N,repetitions);

commInfo_multi.commIndex    = cell(repetitions,1);
commInfo_multi.commChanges  = cell(repetitions,1);
commInfo_multi.commFlexible = cell(repetitions,1);
commInfo_multi.commCohesion = cell(repetitions,1);
commInfo_multi.commDisjoint = cell(repetitions,1);

commFlexible_sum = zeros(N,dT);
CijInfo.all = zeros(N,N,repetitions);

%% Calculate CDFS
for xx = 1:repetitions
    S = community_structure{xx};
    
    [cij,nodeInfo,commInfo] = calc_node_cohesion(S);
    
    Cij_sum = Cij_sum + cij;
    CijInfo.all(:,:,xx) = cij;
    
    cohesion_node_all(:,xx)          = nodeInfo.cohesion_node;
    disjoint_node_all(:,xx)          = nodeInfo.disjoint_node;
    flexibility_node_all(:,xx)	     = nodeInfo.flexibility_node;
    strength_cohesion_node_all(:,xx) = nodeInfo.strength_cohesion_node;
    
    % Dynamic Community Changes Info
    commInfo_multi.commIndex{xx}    = commInfo.commIndex;
    commInfo_multi.commChanges{xx}  = commInfo.commChanges;
    commInfo_multi.commFlexible{xx} = commInfo.commFlexible;
    commInfo_multi.commCohesion{xx} = commInfo.commCohesion;
    commInfo_multi.commDisjoint{xx} = commInfo.commDisjoint;
    
    commFlexible_sum = commFlexible_sum + commInfo.commFlexible;
end

Cij = Cij_sum./repetitions;
commFlexible = commFlexible_sum./repetitions;

cohesion_node_mean          = mean(cohesion_node_all,2);
disjoint_node_mean          = mean(disjoint_node_all,2);
flexibility_node_mean       = mean(flexibility_node_all,2);
strength_cohesion_node_mean = mean(strength_cohesion_node_all,2);

% Global statistics: mean, standard deviation, variance
nodeInfo_multi.global.mean.cohesion          = mean(cohesion_node_mean);
nodeInfo_multi.global.mean.disjoint          = mean(disjoint_node_mean);
nodeInfo_multi.global.mean.flexibility       = mean(flexibility_node_mean);
nodeInfo_multi.global.mean.strength_cohesion = mean(strength_cohesion_node_mean);

nodeInfo_multi.global.std.cohesion          = std(cohesion_node_mean);
nodeInfo_multi.global.std.disjoint          = std(disjoint_node_mean);
nodeInfo_multi.global.std.flexibility       = std(flexibility_node_mean);
nodeInfo_multi.global.std.strength_cohesion = std(strength_cohesion_node_mean);

nodeInfo_multi.global.var.cohesion          = var(cohesion_node_mean);
nodeInfo_multi.global.var.disjoint          = var(disjoint_node_mean);
nodeInfo_multi.global.var.flexibility       = var(flexibility_node_mean);
nodeInfo_multi.global.var.strength_cohesion = var(strength_cohesion_node_mean);

% Node statistics (across simulations): mean, standard deviation, variance
nodeInfo_multi.node.mean.cohesion          = mean(cohesion_node_all,2);
nodeInfo_multi.node.mean.disjoint          = mean(disjoint_node_all,2);
nodeInfo_multi.node.mean.flexibility       = mean(flexibility_node_all,2);
nodeInfo_multi.node.mean.strength_cohesion = mean(strength_cohesion_node_all,2);

nodeInfo_multi.node.std.cohesion          = std(cohesion_node_all,0,2);
nodeInfo_multi.node.std.disjoint          = std(disjoint_node_all,0,2);
nodeInfo_multi.node.std.flexibility       = std(flexibility_node_all,0,2);
nodeInfo_multi.node.std.strength_cohesion = std(strength_cohesion_node_all,0,2);

nodeInfo_multi.node.var.cohesion          = var(cohesion_node_all,0,2);
nodeInfo_multi.node.var.disjoint          = var(disjoint_node_all,0,2);
nodeInfo_multi.node.var.flexibility       = var(flexibility_node_all,0,2);
nodeInfo_multi.node.var.strength_cohesion = var(strength_cohesion_node_all,0,2);

% Simulation data: Output from all simulations
nodeInfo_multi.all.cohesion          = cohesion_node_all;
nodeInfo_multi.all.disjoint          = disjoint_node_all;
nodeInfo_multi.all.flexibility       = flexibility_node_all;
nodeInfo_multi.all.strength_cohesion = strength_cohesion_node_all;

% Cohesion Matrix Info
CijInfo.mean = mean(CijInfo.all,3);
CijInfo.std  = std(CijInfo.all,0,3);
CijInfo.var  = var(CijInfo.all,0,3);