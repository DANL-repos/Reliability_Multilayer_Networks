% Convert time series to multilayer networks
clear
clc

addpath(genpath('./functions'))
addpath(genpath('./grinsted-wavelet-coherence-gus'))

% Define input and output directory
inputDir='./example_data/time_series/';
outputDir='./example_data/Aij/';
if ~exist(outputDir, 'dir')
    mkdir(outputDir)
end

fileList = qbin_create_fileList(inputDir,'file','csv');

% Define variables for multilayer network construction
TR=1.450;
passBand=[0.01 0.1];
windowSize=68; % in TR, window size ~100s
windowOverlap=0;
analysisType='wtc';

% multilayer network construction
parfor i=1:length(fileList)
    file=fileList{i};
    disp(['Working on ', file])
    [filepath,name,ext]=fileparts(file);
    TS=load(file)';
    
    AijFile=[outputDir, 'Aij_', name];
    
    if ~exist([AijFile, '.mat'], 'file')
        tic
        Aij=run_timeSeries2matFastCohi(TS, TR, passBand, windowSize, windowOverlap, analysisType, outputDir);
        toc
       
        parsave(AijFile, Aij)
    else
        disp([AijFile, '.mat file already exist.'])
    end
end
