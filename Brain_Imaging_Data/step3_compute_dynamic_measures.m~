% Compute flexibility, allegiance matrix, integration, and recruitment

clear
clc

addpath(genpath('./functions'))

taskList={'rest', 'inscapes', 'movie', 'flanker'};
sub='M00472509';

% define input and output directory
inputDir='./example_data/SRep/';
outputDirFlexibility='./example_data/flexibility/';
outputDirAlgnc='./example_data/Algnc/';
outputDirI='./example_data/integration/';
outputDirR='./example_data/recruitment/';

if ~exist(outputDirFlexibility, 'dir')
    mkdir(outputDirFlexibility)
end

if ~exist(outputDirAlgnc, 'dir')
    mkdir(outputDirAlgnc)
end

if ~exist(outputDirI, 'dir')
    mkdir(outputDirI)
end
if ~exist(outputDirR, 'dir')
    mkdir(outputDirR)
end

% Parameters used for community detection
gammaList=1.05;
omegaList=2.5;
% repeat the algorithm 100 times
repT=100;

% Yeo 7 networks were used to define integration and recruitment
tmp=load('./example_data/YeoNetworkIndex.mat');
YeoIndex=tmp.networkIndex;
numNetwork=max(YeoIndex);
nROI=length(YeoIndex);


for i=1:length(taskList)
    task=char(taskList{i})
    
    for j=1:length(gammaList)
        gamma=gammaList(j);
        
        for k=1:length(omegaList)
            omega=omegaList(k);
            
            disp(['Working on gamma ', num2str(gamma), ' omega ', num2str(omega)])
            
            file=[inputDir, 'SRep_Aij_', sub, '_ssc_2_', task, '_gamma', num2str(gamma), '_omega', num2str(omega), '.mat'];
            
            
            if exist(file, 'file')
                disp(['Working on ', file])
                [filepath,name,ext]=fileparts(file);
                tmp=load(file);
                SRep=tmp.SRep;
                
                % compute flexibility
                disp('Compute flexibility.')
                outputNameFlexibility=[outputDirFlexibility, '/flexibility_', name, '.txt'];
                
                F=zeros(nROI, length(SRep));
                for p=1:repT
                    S=squeeze(SRep{p});
                    S=S';
                    F(:, p)=flexibility(S);
                end
                avgF=mean(F, 2);
                
                disp('Save flexibility measure.')
                save(outputNameFlexibility, '-ascii','-tabs', 'avgF')
                
                % compute Allegiance
                disp('Compute Allegiance.')
                outputNameAlgnc=[outputDirAlgnc, '/Algnc_', name, '.mat'];
                SRepNew=[];
                for p=1:length(SRep)
                    S=SRep{p};
                    SRepNew{p,1}=squeeze(S);
                end
                figureFlag=0;
                Algnc = create_module_allegiance_matrix(SRepNew,figureFlag);
                disp('Save Allegiance Matrix.')
                save(outputNameAlgnc, 'Algnc')
                
                % compute integration and recruitment
                disp('Compute integration and recruitment.')
                outputNameR=[outputDirR, '/recruitment_Algnc_', name, '.txt'];
                outputNameI=[outputDirI, '/integration_Algnc_', name, '.txt'];
                R = recruitment(Algnc,YeoIndex);
                I = integration(Algnc,YeoIndex);
                disp('Save integration and recruitment.')
                save(outputNameI, '-ascii', '-tabs','I')
                save(outputNameR, '-ascii', '-tabs', 'R')
                
            else
                disp([file, ' does not exist.'])
            end
        end
    end
end



