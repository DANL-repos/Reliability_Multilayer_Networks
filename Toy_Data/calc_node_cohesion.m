function [Cij,nodeInfo,commInfo] = calc_node_cohesion(S,varargin)
% PURPOSE:	Calculates node cohesion of a node in a multilayer network
%           Node cohesion is a measure of how often a node changes from one
%           community to another with other nodes. Highly cohesive nodes
%           change communities with other nodes often.
%
% INPUT:
%	S: Multilayer community structure with community assignment, rows
%      represent each node while columns denote layer
%
%   Optional Inputs:
%       'figure': Flag for showing output figures; three output
%           0: No figures (default)
%           1: Show figures
%           	Figure 1. Community structure of input S
%           	Figure 2. Subplot image showing community changes
%               	a. All unique community changes
%               	b. Flexibile community changes
%               	c. Cohesive community changes
%               	d. Disjoint community changes
%            	Figure 3. Cohesion Matrix
%
%       'colormap': Colormap used to display figures (Default: 'jet')          
%
% OUTPUT:
%
%	* Cij: Cohesion matrix, similar to adjacency matrix but edge weights
%          represent the number of times two nodes change communities
%          mutually. Measure is normalized by the number of times two nodes
%          could change communities mutually.
%
%	* nodeInfo: STRUCT variable containing dynamic measures of nodes:
%
%         + cohesion_node: Node cohesiveness for a given node
%            
%                             # Times Node Changes Communities Mutually
%             Cohesiveness = -------------------------------------------
%                             # Possible Times Nodes Change Communities
%
%         + disjoint_node: Node disjointedness for a given node
%            
%                            # Times Node Changes Communities Independently
%           Disjointedness = ----------------------------------------------
%                              # Possible Times Nodes Change Communities
%
%         + flexibility_node: Node flexibility for a given node
%            
%                                # Times Node Changes Communities
%              Flexibility = -------------------------------------------
%                             # Possible Times Nodes Change Communities
%
%         + strength_cohesion_node: Cohesion strength for a given node (sum
%                                   of entries in cohesion matrix)
%
%   * commInfo: STRUCT variable containing space-time diagrams of community
%               changes. In array, rows represent nodes, columns represent
%               layers/windows. STRUCT also includes reference index, that
%               describes all unique community changes.
%
%         + commIndex: Array showing unique community shifts.
%                      Column 1 - Indicates initial community of node
%                      Column 2 - Indicates new community of node
%                      Column 3 - Value for community shift in commChanges
%
%         + commChanges:  Array contains all community changes
%
%         + commCohesion: Array contains all mutual community changes
%
%         + commDisjoint: Array contains all independent community changes
%
%         + commFlexible: Array contains all flexible community changes
%
%--------------------------------------------------------------------------
%
%  Author: Qawi Telesford
%    Date: 2019-04-17
% Version: 2.1
%
%	History:	
%               2.1 (2019-04-17) | QKT
%               * Streamlined output to simplify as STRUCT variables
%
%               2.0 (2018-08-03) | QKT
%               * Overhauled community changes portion of code. Previous
%                 version created table for all possible community changes,
%                 which is unwieldy for network with large number of
%                 communities. Changed so that only unique community
%                 changes are added to commIndex
%               * Added community changes subfunction
%               * Added flexible community changes to output
%               * Altered cohesion strength matrix to include number of
%                 times two nodes change together with respect to the
%                 number of layers
%
%               1.2 (2015-06-13) | QKT
%               * Removed cohesion normalization
%               * Updated cohesion calculation, now calculates only when
%                 node changes with other nodes
%               * Added node disjointedness, a measure of how often a node
%                 changes communities by itself
%               * Added community change arrays and index
%
%               1.1 (2015-01-16) | QKT
%               * Updated cohesion normalization
%
%               1.0 (2014-10-22) | QKT
%               * Initial release
%
%--------------------------------------------------------------------------
%% Options Parsing
% Set defaults for figures (optional)
defaultFigure = 0;
defaultColormap = 'jet';

% Set input parser names
p = inputParser;
addRequired(p,'community_structure',@(x) true);

% Default input parameters
addParameter(p,'figure',defaultFigure,@(x) isnumeric(x) && isscalar(x))
addParameter(p,'colormap',defaultColormap,@ischar);

parse(p,S,varargin{:});

S = p.Results.community_structure;
figureFlag = p.Results.figure;
colormapFig = p.Results.colormap;

%% Error Checking
if(size(S,2)<2)
    error('Input matrix invalid or contains one layer, please input multilayer community structure matrix.');
end

% Check if matrix contains values of community assignment
S_check = sum(sum(S - double(int16(S))));

if(S_check ~= 0)
    error('Input matrix contains non-integer values, please input multilayer community structure matrix with integer values.');
end

% if(nargin < 2 || isempty(options))
%     options.figureFlag	= 0;
%     options.colormap	= 'jet';
% end


%% Initialize Cohesion Matrices and Vectors
N	= size(S,1);                % network size
dT  = size(S,2) - 1;
Cij	= zeros(N);                 % cohesion matrix

%% Community Structure Changes and Cohesion/Disjoint Arrays
[commChanges,commIndex] = calc_community_change(S);

commCohesion = commChanges;
commDisjoint = commChanges;

for tt = 1:dT
    for idx = 1:size(commIndex,1)
        commCount = sum(commChanges(:,tt)==commIndex(idx,3));
        if(commCount > 0)
            if(commCount == 1)
                commCohesion(commChanges(:,tt)==commIndex(idx,3),tt) = 0;
            else
                commDisjoint(commChanges(:,tt)==commIndex(idx,3),tt) = 0;
            end
        end
    end
end

%% Calculate node flexibility, cohesiveness, and disjointedness
commFlexible   = logical(commChanges);  % Indicates if node changed communities
cohesionChange = logical(commCohesion); % Indicates if node changed communities mutually
disjointChange = logical(commDisjoint); % Indicates if node changed communities independently

nodeChangeAll = sum(commFlexible,2);       % Number of times a node changed communities
nodeChangeCoh = sum(cohesionChange,2);
nodeChangeDis = sum(disjointChange,2);

flexibility_node = nodeChangeAll./dT; % node flexibility
cohesion_node    = nodeChangeCoh./dT; % node cohesiveness
disjoint_node    = nodeChangeDis./dT; % node disjointedness

%% Cohesion array

% Finds all unique node changes in each layer
communityShift = unique(commCohesion);
communityShift(communityShift==0) = []; % Remove nodes that do not change communities
cohesion = cell(size(communityShift,1),size(commChanges,2)); % Cell matrix indicating which nodes made unique community shift together

% Finds all unique node changes in each layer
if(isempty(communityShift))
    % disp('Multilayer network contains no changes in community assignment across layers.');
    for xx = 1:size(cohesion,1)
        for yy = 1:size(cohesion,2)
            cohesion{xx,yy} = 0;
        end
    end
else
    for xx = 1:size(cohesion,1)
        communityShiftXX = communityShift(xx);
        for yy = 1:size(cohesion,2)
            if(communityShiftXX == 0)
                cohesion{xx,yy} = 0; % nodes did not change communities
            else
                findCohesion = find(commChanges(:,yy)==communityShiftXX);
                
                if(isempty(findCohesion) || length(findCohesion)==1)
                    cohesion{xx,yy} = 0; % layer without unique community shift or single node community changes set to 0
                else
                    cohesion{xx,yy} = findCohesion;
                end
            end
        end
    end
end

%% Cohesion Matrix
cij = zeros(N);

for ii = 1:size(cohesion,1)
    for jj = 1:size(cohesion,2)
        if(cohesion{ii,jj} ~= 0)
            [p,q] = meshgrid(cohesion{ii,jj}, cohesion{ii,jj});
            node_pairs = [p(:) q(:)];
            
            cij_temp = node2mat(node_pairs,N);
            cij = cij + cij_temp;
        end
    end
end

Cij = Cij + cij;
Cij(logical(eye(size(Cij)))) = 0;
Cij = Cij./dT;

%% Calculate Change Frequnecy
strength_cohesion_node = sum(Cij,2);

%% Figures
if(figureFlag == 1)
    figure(1); subplot(1,1,1); imagesc(S);                                    axis('square'); xlabel('Layer'); ylabel('Nodes'); colormap(colormapFig); colorbar; title('Multilayer Community Structure');
    
    figure(2); subplot(2,2,1); imagesc(commChanges, [1 max(commIndex(:,3))]); axis('square'); xlabel('Layer'); ylabel('Nodes'); colormap(colormapFig); colorbar; title('Community Changes');
    figure(2); subplot(2,2,2); imagesc(commFlexible);                         axis('square'); xlabel('Layer'); ylabel('Nodes'); colormap(colormapFig); colorbar; title('Flexible Community Changes');
    figure(2); subplot(2,2,3); imagesc(commCohesion,[1 max(commIndex(:,3))]); axis('square'); xlabel('Layer'); ylabel('Nodes'); colormap(colormapFig); colorbar; title('Cohesive Community Changes');
    figure(2); subplot(2,2,4); imagesc(commDisjoint,[1 max(commIndex(:,3))]); axis('square'); xlabel('Layer'); ylabel('Nodes'); colormap(colormapFig); colorbar; title('Disjoint Community Changes');
    
    figure(3); subplot(1,1,1); imagesc(Cij); axis('square','off'); colormap(colormapFig); colorbar; title('Cohesion Matrix');
end

%% Output STRUCT Variables
% Dynamic Node Info
nodeInfo.cohesion_node = cohesion_node;
nodeInfo.disjoint_node = disjoint_node;
nodeInfo.flexibility_node = flexibility_node;
nodeInfo.strength_cohesion_node = strength_cohesion_node;

% Dynamic Community Changes Info
commInfo.commIndex = commIndex;
commInfo.commChanges = commChanges;
commInfo.commCohesion = commCohesion;
commInfo.commDisjoint = commDisjoint;
commInfo.commFlexible = commFlexible;

end
%% Subfunction: Edge List to Adjacency Matrix Conversion
function Aij = node2mat(node_list,N)
% This subfunction creates an adjacency matrix from an input edge list

Aij_nodes = false(N);

for ii = 1:length(node_list)
    Aij_nodes(node_list(ii,1),node_list(ii,2)) = 1;
    Aij_nodes(node_list(ii,2),node_list(ii,1)) = 1;
end
Aij = Aij_nodes;
end

%% Subfunction: Community Change Arrays and Index
function [commChanges,commIndex] = calc_community_change(S)
N = size(S,1); % Number of nodes
T = size(S,2); % Time points/windows
dT = T - 1; % Number of possible changes
commChanges_ii = zeros(N,dT);
commIndex_ii = [];
for tt = 1:dT
    % Take individual columns in community network
    S_T = [S(:,tt), S(:,tt+1)];
    
    % Find where nodes changes occur and select only those rows
    S_shift = abs(S_T(:,2) - S_T(:,1));
    nodeChanges = S_shift~=0;
    S_dT = S_T(nodeChanges,:);
    
    % Identify unique community shifts
    commShiftUnique = unique(S_dT,'rows');
    
    % Generate list of community changes
    if(isempty(commIndex_ii))
        commShiftIndex = (1:size(commShiftUnique,1))';
        commIndex_ii = [commShiftUnique, commShiftIndex];
    else
        for xx = 1:size(commShiftUnique,1)
            % Check if unique community shift exists, in commShiftList, if
            % it does not exist, add it to the list, else ignore
            [~,~,commShiftIntersect] = intersect(commShiftUnique(xx,:),commIndex_ii(:,1:2),'rows');
            if(isempty(commShiftIntersect))
                commShiftIndex = max(commIndex_ii(:,3))+1;
                commIndex_ii = [commIndex_ii; commShiftUnique(xx,:), commShiftIndex]; %#ok<AGROW>
            end
        end
    end
    
    % Fill in community changes
    for yy = 1:size(commIndex_ii,1)
        nodeShiftIndex = ismember(S_T,commIndex_ii(yy,1:2),'rows');
        commChanges_ii(nodeShiftIndex,tt) = commIndex_ii(yy,3);
    end
end

% Sort and reindex community changes
commIndex_sort = sortrows(commIndex_ii);
commIndex = [commIndex_sort(:,1:2),commIndex_ii(:,3)];

commChanges = zeros(size(commChanges_ii));
for cc = 1:(size(commIndex_ii,1))
    commChanges(commChanges_ii==cc) = commIndex_sort(cc,3);
end
end