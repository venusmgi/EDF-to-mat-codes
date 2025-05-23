function Process_EDF_To_Mat (desiredSamplingRate,desiredChannelOrder,removableChannels, channelsToBeReplaced, ...
    newChannelNames,channInfoName,removeChannels, renameChannels, headerFormat)


% Previously was called Convert_EDFtomat_Desired_Fs_channOrder_outputExcel_replaceChans
% then was renamed to Convert_and_Standerdize_EDF_to_MAT
% This function processes EDF files in two modes:
% 1. ParentPath mode: Processes all EDF files in a parent directory structure
%    (diagnosis and follow-up folders for each patient)
% 2. EDF mode: Processes selected EDF file(s) directly
%
% The function performs the following operations:
% - Reads EDF files
% - Reorders channels according to specified order
% - Removes specified channels
% - Replaces specified channels with new names
% - Resamples data to desired sampling rate
% - Saves processed data in MAT format
%
% Inputs:
%   desiredSamplingRate - Target sampling rate for the output (Hz)
%   desiredChannelOrder - Cell array of channel names in desired order
%   removableChannels - Cell array of channels to be excluded
%   channelsToBeReplaced - Cell array of channel names to be replaced
%   newChannelNames - Cell array of new names for channels to be replaced
%   channInfoName - Name for the channel information output file
%   if you chose 'EDF', you can direcly select the EDF files you want to 
%   convert from any disred path. If you choose'ParentPath', the edf files
%   should be organized as below:
%   -Parent folder with a name of your choosing(e.g. 10.CHOC)
%       -patient folder with a name of your choosing(e.g. 10-0001 (2017))
%           -a folder named 'diagnosis'
%              edfs in are(not) in 'diagnosis' folder 
%           -a folder named'follow up'
%              edfs in are(not) in 'diagnosis' folder
%                  
%   headerFormat - 'EEGlab' for EEGlab format, empty for default format
%
%   removeChannels - 'remove' if there is any spesific channel you want to 
%   remove, input the channel names that you want to remove in 
%   "removableChannels" If there is no spesific channel, and you input 
%   'remove' it will automatically remove the channels that are pre-defined 
%   in "Remove_Channel_Names.m" function.
%
%   renameChannels - 'rename' if there is any spesific channel you want to 
%   rename and input the channel names that you want to rename in 
%   "channelsToBeReplaced" and their new names in "newChannelNames". 
%   If there is no spesific channel that you wat to rename and you input
%  'rename' it will automatically rename the channels that are pre-defined 
%   in "Rename_Channel_Labels.m " function.
%
% Outputs:
%   Saves processed data as .mat files in the same directory as input files
%   If header_format is 'EEGlab', saves in EEGlab format
%   Otherwise, saves in standard format with reordered_record and reordered_hdr
%   It also a second mat file that records the changes that has been done to
%   channel names
%
% Example:
%   EDF_Processor(200, {'Fp1', 'Fp2', 'F3'}, {'ECG'}, {'Fp1'}, {'Fp1_new'}, 'channel_info.mat', 'ParentPath', '')
%   EDF_Processor(200, {'Fp1', 'Fp2', 'F3'}, {'ECG'}, {'Fp1'}, {'Fp1_new'}, 'channel_info.mat', 'EDF', 'EEGlab')



% Add current directory to MATLAB path
functionPath = pwd;
addpath(functionPath)

% Set default header format if not provided
if nargin < 9
    headerFormat = '';  % Default value (regular header format)
end


%  ask user to choose the batch processing mode of EDF
 selectionMode = questdlg('Select a processing Mode:', 'Processing Mode', ...
        'ParentPath', 'EDF(s)','ParentPath');


% Process based on selected mode
switch lower(selectionMode)
    case 'parentpath'
        % Get parent direcotry from the user
        ParentPath = uigetdir('', 'Select Parent Directory');
        if ParentPath == 0
            fprintf('Operation canceled by the user.\n')
            return;
        end

        % Get list of patient folders
        patientFolders = Find_Folders(ParentPath);

        % Process each patient folder
        for i = 1:length(patientFolders)
            fprintf('Processing patient folder: %s\n', patientFolders{i});

            % Process diagnosis folder
            Process_Folder(ParentPath, patientFolders{i}, 'diagnosis', desiredSamplingRate, ...
                desiredChannelOrder, removableChannels, channelsToBeReplaced, ...
                newChannelNames, channInfoName,removeChannels, renameChannels,headerFormat);

            % Process follow-up folder
            Process_Folder(ParentPath, patientFolders{i}, 'follow up', desiredSamplingRate, ...
                desiredChannelOrder, removableChannels, channelsToBeReplaced, ...
                newChannelNames, channInfoName,removeChannels, renameChannels, headerFormat);
        end

        fprintf('All files processed successfully.\n');


    case 'edf(s)'

        % Get EDF file(s) from user
        [fileNames, filePath] = uigetfile('*.edf', 'Select EDF file(s)', 'MultiSelect', 'on');

        if isequal(fileNames,0)
            fprintf('Operation canceled by the user.\n')
            return;
        end

        % Convert to cell array if single file selected
        if ~iscell(fileNames)
            fileNames = {fileNames};
        end

        % Process each selected file
        for i = 1:length(fileNames)
            fprintf('Processing file %d%d: %s',i,length(fileNames),fileNames{i})

            EDF_To_Mat_With_Channel_Standardizing (filePath,fileNames{i},desiredSamplingRate,desiredChannelOrder, ...
                removableChannels,channelsToBeReplaced, newChannelNames,filePath,channInfoName,removeChannels,renameChannels,headerFormat)

        end


end

end



%% subfunctions



function [folderNames] = Find_Folders(ParentPath)
% FIND_FOLDERS Finds all non-hidden folders in the specified parent path
%
% Inputs:
%   ParentPath - Path to the parent directory
%
% Outputs:
%   folderNames - Cell array of folder names found in the parent directory

j = 1;
folderNames = {};
allFilesAndFolders = dir(ParentPath);

% Process each item in the directory
for i = 1: length(allFilesAndFolders)


    % Get the folder name
    myFolderName = allFilesAndFolders(i).name;

    % Skip non-folders and hidden folders
    if allFilesAndFolders(i).isdir &  ~((strcmpi(myFolderName, '..')) || (strcmpi(myFolderName, '.')))

        % Get the folder name
        folderNames{j} = allFilesAndFolders(i).name;
        j = j+1;


    end
end
end


function Process_Folder(ParentPath, patientFolder, folderType, desiredSamplingRate, ...
    desiredChannelOrder, removableChannels, channelsToBeReplaced, ...
    newChannelNames, channInfoName,removeChannels, renameChannels, headerFormat)
    % This helper function processes all EDF files in a specified folder (diagnosis or follow-up)
    % and converts them to MAT format with the specified parameters.
    %
    % Inputs:
    %   ParentPath - Parent directory path
    %   patient_folder - Name of the patient folder
    %   folder_type - Type of folder ('diagnosis' or 'follow up')
    %   desired_sampling_rate - Target sampling rate
    %   desired_channel_order - Cell array of channel names in desired order
    %   removable_channels - Cell array of channels to be excluded
    %   channel_tobe_replaced - Cell array of channel names to be replaced
    %   new_channel_names - Cell array of new names for channels to be replaced
    %   header_format - 'EEGlab' for EEGlab format, empty for default format
    %   chann_info_name - Name for the channel information output file
    
    % Construct folder path
    folder_path = fullfile(ParentPath, patientFolder, folderType);
    if ~exist(folder_path, 'dir')
        fprintf('Warning: %s folder not found for patient %s\n', folderType, patientFolder);
        return;
    end
    
    % Change to folder and get EDF files
    cd(folder_path);
    edf_files = dir('*.edf');
    if isempty(edf_files)
        fprintf('No EDF files found in %s folder for patient %s\n', folderType, patientFolder);
        return;
    end
    
    % Process each EDF file
    for j = 1:length(edf_files)
        edf_name = edf_files(j).name;
        fprintf('Processing file: %s\n', edf_name);
    
        % Check if output files already exist
        matFileName1 = strrep(edf_name, '.edf', '_reordered.mat');
        matFileName2 = strrep(edf_name, '.edf', '_reordered_resampled.mat');
        matFilePath1 = fullfile(folder_path, matFileName1);
        matFilePath2 = fullfile(folder_path, matFileName2);
    
    
        if exist(matFilePath1, 'file') || exist(matFilePath2, 'file')
            fprintf('Skipping %s - output files already exist\n', edf_name);
            continue;
            % Convert the file
        end
    
        EDF_To_Mat_With_Channel_Standardizing(folder_path, edf_name, desiredSamplingRate, ...
            desiredChannelOrder, removableChannels, channelsToBeReplaced, ...
            newChannelNames, ParentPath, channInfoName, removeChannels, renameChannels, headerFormat);
    
    
    end
end



