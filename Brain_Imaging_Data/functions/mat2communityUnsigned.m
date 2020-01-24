function [SRep,Q,B] = mat2communityUnsigned(A, gamma, omega,repT,varargin)
% Input:
% A is a mxn cell array, 
%            m is the number of sessions
%            n is the number of slices within each session
%            each slice is a connectivity matrix
% gamma: resolution parameter, omega: coupling parameter
% repT is the repeated times of community detection
% Output:
% SRep is repT x 1 cell array of community detection results
% Q is repT x 1 array of associated modularity values of SRep
SRep = cell(repT,1);
Q = zeros(repT,1);
N=length(A{1,1});
[nS, T] = size(A);
B=spalloc(N*nS*T,N*nS*T,N*N*nS*T+2*N*nS*T);
twomu=0;
for s = 1:nS
    for t=1:T
        k=sum(A{s,t});
        twom=sum(k);
        twomu=twomu+twom;
        indx=[1:N]+((s-1)*T+t-1)*N;
        B(indx,indx)=A{s,t}-gamma*k'*k/twom;
    end
end
twomu=twomu+2*omega*N*(T-1);
diagTerm = ones(N*nS*T,2);
% define a scaling of cross-section omega, if it is zero, different
% sessions are then isoloted; if it is one, different sessions are then
% concatenated without penalty
if nargin == 4
    crossSectionOmegaScale = 0.5;
end
if nargin == 5
    crossSectionOmegaScale = varargin{1};
end
for i = 1:nS-1
    diagTerm(N*i*T-N+1:N*i*T,1) = crossSectionOmegaScale;
end
diagB = spdiags(diagTerm,[-N,N],N*nS*T,N*nS*T);
diagB = min(diagB ,diagB');
B = B + omega*diagB;
% for the slices between two sessions, decrease its interlayer connection
% into half omega

% repeatedly run genlouvain to obtain multiple results as the genlouvain
% has randomness
for i = 1:repT
    % B,limit,verbose,randord,randmove,S0
    [S,QTemp] = genlouvain(B,[],0,1,'moverandw');
    Q(i) = QTemp/twomu;
    SRep{i} = reshape(S,N,T,nS);
end
