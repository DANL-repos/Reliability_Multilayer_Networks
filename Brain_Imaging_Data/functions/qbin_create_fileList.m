function [filePath,fileName,fileDate,fileExtn] = qbin_create_fileList(inputDir,options,fileType,suppressFlag)
%% PURPOSE:	Creates list of files or subdirectories from input directory
%
% INPUT:
%           inputDir: Input directory specified by user
%                     (Default: current working directory)
%
%           options: Option to specify list of files or subdirectories
%                    'file' creates list of files
%                    'dir' creates list of subdirectories
%                    (Default: 'file')
%
%           fileType: File extension to include in list of files, input as
%                     '.ext' or 'ext'
%                     (Default: All files)
%
% OUTPUT:
%           filePath: List with full path of files or subdirectories
%
%           fileName: List with names of files or subdirectories
%
%           fileDate: List of system timestamps for files or subdirectories
%
%           fileExtn: List of file extensions, list empty if subdirectory
%
% EXAMPLES:
%           INPUT:
%               inputDir = 'C:\'; options = 'file';
%               filePath = qbin_create_fileList(inputDir,options);
%
%           OUTPUT:
%               Full path list of all files in directory: 'C:\'
%
%           INPUT:
%               inputDir = 'C:\'; options = 'dir';
%               filePath = qbin_create_fileList(inputDir,options);
%
%           OUTPUT:
%               Full path list of all subdirectories in directory: 'C:\'
%
%           INPUT:
%               inputDir = '[]'; options = 'file'; fileType = '.txt';
%               filePath = qbin_create_fileList(inputDir,options,fileType);
%                                           OR
%               filePath = qbin_create_fileList([],options,fileType);
%
%           OUTPUT:
%               Full path list of all files with file extension '.txt.' in
%               current working directory
%
% CHANGE/VERSION LOG
% Author: Qawi Telesford
%
% Date: 2017-12-08
%   Version: 2.0
%	Update(s):
%               - User can input specific file extension for list of files
%               - Added error checking for input variables
%               - Modified to make user input and function output modular,
%                 original code combined user input and list creation, now
%                 function checks user input, sets any defaults if
%                 necessary, then creates list of files/subdirectories
%
% Date: 2018-01-25
%	Version: 1.2
%   Update(s):
%               - Added text suppression for empty directories.
%
% Date: 2014-02-04
%	Version: 1.1
%   Update(s):
%               - Documentation update
%
% Date: 2013-06-24
%	Version: 1.0 (Initial Release)
%% Error checking
% Check user input for input directory, returns an error if not a valid
% directory. If no directory entered, function uses current directory
if(nargin < 1 || isempty(inputDir))
    % error('Please enter input directory.');
    disp('User did not enter input directory, using current working directory.');
    inputDir = pwd;
else
    if(~isdir(inputDir))
        error('User input not a directory, please enter input directory.');
    end
end

% Check user input for "options"
if(nargin < 2 || isempty(options))
    options = 'file';
else
    % If string in "options" does not match 'file' or 'dir', 'file' is used
    if(~strcmpi(options,'dir') && ~strcmpi(options,'file'))
        options = 'file';
        disp('User input for options not valid, using default: ''file''');
    end
end

% Check user input "fileType"
% If no file extension is given, function includes all files in directory
if(nargin < 3 || isempty(fileType))
    fileType = [];
else
    if(strcmp(options,'file'))
        % Adds '.' to file extension if not included in "fileType"
        if(~strcmp(fileType(1),'.'))
            fileType = ['.',fileType];
        end
    elseif(strcmp(options,'dir'))
        % If 'dir' specified in "options", function ignores "fileType"
        disp('User requested directory list and specific file type, function will only output directory list.');
        fileType = [];
    end
end

% Text suppression flag
if(nargin < 4 || isempty(suppressFlag))
   suppressFlag = 0;
end

%% Initialize file lists
filePath = {}; % Full path of file or subdirectory
fileName = {}; % Name of file of subdirectory
fileDate = []; % Timestamp of file
fileExtn = {}; % Extension of filename, empty if subdirectory

%% Create list of files or directories
dirList = dir(inputDir); % STRUCT list of files/subdirectories in inputDir

if(strcmp(options,'file'))
    for ii = 1:size(dirList,1)
        if(dirList(ii,1).isdir~=1)
            [~,~,fileExt] = fileparts(dirList(ii,1).name);
            
            % Checks if file extensions matches input "fileType"
            % includeFlag: specifies if file should be included in list
            if(~isempty(fileType))
                if(strcmpi(fileType,fileExt))
                    includeFlag = 1;
                else
                    includeFlag = 0;
                end
            else
                includeFlag = 1;
            end
            
            if(includeFlag == 1)
                filePath = [filePath; fullfile(inputDir,dirList(ii,1).name)]; %#ok<*AGROW>
                fileName = [fileName; dirList(ii,1).name];
                fileDate = [fileDate; dirList(ii,1).datenum];
                fileExtn = [fileExtn; fileExt];
            end
        end
    end
    
    % Display if no files or matching file extension found in directory
    if(suppressFlag == 1)
        if(isempty(filePath))
            if(~isempty(fileType))
                disp(['No files matching file extenstion (',fileType,') found in input directory (',inputDir,')']);
            else
                disp(['No files found in input directory (',inputDir,')']);
            end
        end
    end
elseif(strcmp(options,'dir'))
    for ii = 1:size(dirList,1)
        if(dirList(ii,1).isdir==1 && ~strcmp(dirList(ii,1).name,'.') && ~strcmp(dirList(ii,1).name,'..'))
            filePath = [filePath; fullfile(inputDir,dirList(ii,1).name)];
            fileName = [fileName; dirList(ii,1).name];
            fileDate = [fileDate; dirList(ii,1).datenum];
        end
    end
    
    % Display if no subdirectories found in directory
    if(suppressFlag == 1)
        if(isempty(filePath))
            disp(['No subdirectories found in input directory (',inputDir,')']);
        end
    end
end