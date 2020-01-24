function [S,Q] = genlouvain_multi(Aij,repetitions,gamma,omega)
% PURPOSE:	This function takes an input multilayer network and runs the
%           generalized Louvain algorithm on the network with user-defined
%           number of repetitions.
%
% INPUT:
%         Aij:	Multilayer network (N×1 or 1×N cell array)
%
% repetitions:	Number of simulations run on network
%               (Default: 100)
%
%       gamma: 	Resolution paramete for intraslice connections
%               (Default: 1)
%
%       omega:  Resolution parameter for interslice connections
%               (Default: 1)
%
%	outputDir:  Output directory for adjacency matrix MAT file
%               (Default: Current directory)
%
% OUTPUT:
%           S:	Cell array of community structure for each simulation
%
%           Q:  Modularity value for each simulation
%
%--------------------------------------------------------------------------
%
%    Author: Qawi Telesford
%      Date: 2018-08-06 | QKT
%   Version: 1.2
% Update(s):
%            * Update genlouvain function call parameters, should address
%              problems with data output from previous versions where no
%              changes were seen in first and last layer in algorithm
%
%    Author: Qawi Telesford
%      Date: 2018-01-12 | QKT
%   Version: 1.1
% Update(s):
%            * Updated to run new version of genlouvain algorithm
%              (older version used genlouvainmex)
%
%    Author: Qawi Telesford
%      Date: 2015-03-19 | QKT
%   Version: 1.0 (Initial Release)
%
%--------------------------------------------------------------------------
%% Error Checking
if(nargin < 1 || isempty(Aij))
    error('Missing input network, please input Aij.');
end

if(nargin < 2 || isempty(repetitions))
	repetitions = 100;
end

if(nargin < 3 || isempty(gamma))
    gamma = 1;
end

if(nargin < 4 || isempty(omega))
    omega = 1;
end

%% Initialize variables
S = cell(repetitions,1);
Q = zeros(repetitions,1);

%% Run generalized Louvain simulations
for ii = 1:size(S,1)
    [S{ii},Q(ii)] = run_genLouvain(Aij,gamma,omega);
end

%% Generalized Louvain subfunction
function [S,Q] = run_genLouvain(AijIn,gamma,omega)

N=length(AijIn{1});
T=length(AijIn);
B=spalloc(N*T,N*T,N*N*T+2*N*T);
twomu=0;

for s=1:T
    k=sum(AijIn{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=(1:N)+(s-1)*N;
    B(indx,indx)=AijIn{s}-gamma*(k'*k)/twom;
end

twomu=twomu+2*omega*N*(T-1);
B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
[S,Q] = genlouvain(B,[],0,1,'moverandw');
Q = Q/twomu;
Q = full(Q);
S = reshape(S,N,T);