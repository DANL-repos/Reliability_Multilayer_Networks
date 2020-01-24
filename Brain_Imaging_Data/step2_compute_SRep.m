% Perform multilayer community detection

clear
clc

addpath(genpath('./functions'))
addpath(genpath('./GenLouvain2.1.2'))

% Define input and output directory
inputDir='./example_data/';
outputDirSRep=[inputDir, 'SRep'];
outputDirQ=[inputDir, 'Q'];

if ~exist(outputDirSRep, 'dir')
    mkdir(outputDirSRep)
end

if ~exist(outputDirQ, 'dir')
    mkdir(outputDirQ)
end

taskList={'rest', 'inscapes', 'movie', 'flanker'};
sub='M00472509';

% Parameters used for community detection
gammaList=1.05;
omegaList=2.5;
% repeat the algorithm 100 times
repT=100;

for i=1:length(taskList)
    task=char(taskList{i})
    
    for j=1:length(gammaList)
        gamma=gammaList(j);
        
        for k=1:length(omegaList)
            omega=omegaList(k);
            
            disp(['Working on gamma ', num2str(gamma), ' omega ', num2str(omega)])
            
            inputFile=[inputDir, 'Aij/Aij_', sub, '_ssc_2_', task, '.mat'];
            
            if exist(inputFile, 'file')
                disp(['Working on ', inputFile])
                [filepath,name,ext]=fileparts(inputFile);
                tmp=load(inputFile);
                Aij=tmp.x;
                
                % compute modular membership and modularity
                SRepfile=[outputDirSRep, '/SRep_', name, '_gamma', num2str(gamma), '_omega', num2str(omega), '.mat'];
                Qfile=[outputDirQ, '/Q_', name,'_gamma', num2str(gamma), '_omega', num2str(omega),'.mat'];
                
                if ~exist(SRepfile, 'file')
                    disp('Compute modular structure')
                    tic
                    [SRep, Q]=mat2communityUnsigned(Aij, gamma, omega, repT);
                    
                    save(SRepfile, 'SRep')
                    save(Qfile, 'Q')
                    toc
                else
                    disp([SRepfile, ' exist.'])
                end
                
            else
                disp([inputFile, ' does not exist.'])
            end
            
        end
        
    end
    
end



