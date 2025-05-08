function Convert_and_Standerdize_EDF_to_MAT (desiredSamplingRate,desiredChannelOrder,removableChannels, channelsToBeReplaced, ...
    newChannelNames,channInfoName,selectionMode, headerFormat)
% CONVERT_EDF2MAT Converts a single EDF file to MAT format with specified parameters
%
% This function converts an EDF file to MAT format, reordering channels and
% resampling the data according to specified parameters.
%
% Inputs:

%   desiredSamplingRate - Target sampling rate for the output (e.g 200 Hz)
%   desiredChannelOrder - Cell array of channel names in desired order
%   removableChannels - Cell array of channels to be excluded
%   channelsToBeReplaced - Cell array of channel names to be replaced
%   newChannelNames - Cell array of new names for channels to be replaced
%   channInfoName - Name for the channel information output file
%   selectionMode - To convert in edfs in folders or directly selecting EDF
%   headerFormat - 'EEGlab' for EEGlab format, empty for default format
%
% Outputs:
%   Saves a .mat file of re-ordered EEG with either EEGlab format or standard format
%   save a second mat file that records the changes that has been done to
%   channel names

functionPath = pwd;
addpath(functionPath)

if nargin < 8
    headerFormat = '';  % Default value (regular header format)
end

if nargin <7 
    selectionMode = questdlg('Select a processing Mode:', 'Processing Mode', ...
        'ParentPath', 'EDF','ParentPath');
end

switch lower(selectionMode)
    case 'parentpath'
    ParentPath = uigetdir('', 'Select Parent Directory');
    if ParentPath == 0
        fprintf('Operation canceled by the user.\n')
        return;
    end

    patientFolders = Find_Folders(ParentPath);
    % Process each patient folder
    for i = 1:length(patientFolders)
        fprintf('Processing patient folder: %s\n', patientFolders{i});
        
        % Process diagnosis folder
        Process_Folder(ParentPath, patientFolders{i}, 'diagnosis', desiredSamplingRate, ...
            desiredChannelOrder, removableChannels, channelsToBeReplaced, ...
            newChannelNames, headerFormat, channInfoName);
        
        % Process follow-up folder
        Process_Folder(ParentPath, patientFolders{i}, 'follow up', desiredSamplingRate, ...
            desiredChannelOrder, removableChannels, channelsToBeReplaced, ...
            newChannelNames, headerFormat, channInfoName);
    end 
    
    fprintf('All files processed successfully.\n');
       

case 'edf'
    [fileNames, filePath] = uigetfile('*.edf', 'Select EDF file(s)', 'MultiSelect', 'on');
    
    if isequal(fileNames,0)
        fprintf('Operation canceled by the user.\n')
        return;
    end

    % Convert to cell array if single file selected
    if ~iscell(fileNames)
        fileNames = {fileNames};
    end
    for i = 1:length(fileNames)
        fprintf('Processing file %d%d: %s',i,length(fileNames),fileNames{i})

        % ProcessFolder(filePath,)
        Convert_EDF2Mat (filePath,fileNames{i},desiredSamplingRate,desiredChannelOrder, ...
            removableChannels,channelsToBeReplaced, newChannelNames,headerFormat,filePath,channInfoName)

    end



end



end

function Convert_EDF2Mat (filePath,fileName,desiredSamplingRate,desiredChannelOrder,removableChannels,channelsToBeReplaced, ...
    newChannelNames,headerFormat,PathToSaveChannelInfo,channInfoName)
% CONVERT_EDF2MAT Converts an EDF file to MAT format with specified parameters
%
% This function reads an EDF file, processes it according to the specified parameters,
% and saves the results in MAT format.
%
% Inputs:
%   filePath - Path to the folder containing the EDF file
%   fileName - Name of the EDF file
%   desiredSamplingRate - Target sampling rate
%   desiredChannelOrder - Cell array of channel names in desired order
%   removableChannels - Cell array of channels to be excluded
%   channelsToBeReplaced - Cell array of channel names to be replaced
%   newChannelNames - Cell array of new names for channels to be replaced
%   headerFormat - 'EEGlab' for EEGlab format, empty for default format
%   PathToSaveChannelInfo - Directory to save the .mat file that saves the
%   info of changes to the EEG channel names
%   channInfoName - Name for the channel information output file



addpath(filePath)
isEEGlabFormat = strcmpi(headerFormat, 'EEGlab');
numChannels = length(desiredChannelOrder);

%convert edf to .mat
[hdr, eegRecord] = edfreadUntilDone(string(fileName));

%saving header info
EEGOriginalChannelOrder = hdr.label;
% finding the original sampling rate as it is in the header
originalFs = int64(hdr.frequency(1));  %TODO: figure out a better way to find the original fs


% change the channel names and order based on the input
[reordered_record,sortedChannelsIndices] =Standardize_EEG_Channel_Order (desiredChannelOrder,removableChannels,channelsToBeReplaced, ...
    newChannelNames,eegRecord, EEGOriginalChannelOrder,PathToSaveChannelInfo,channInfoName,fileName);


% Resample if needed
if originalFs ~= desiredSamplingRate
    reordered_record = Resample_Data(reordered_record,  desiredSamplingRate, originalFs, numChannels);
    hdr.frequency(1:numChannels) = desiredSamplingRate;
end

% Determine output filename based on resampling
fname = Get_Output_Filename(fileName, originalFs, desiredSamplingRate);


% Save the processed data
fprintf('Saving processed data...\n');

if isEEGlabFormat
%in case we wanted it to be in EEGLab format
% Saving the EEG record of the channel order we want
    reordered_EEG = Make_EEG_Header(hdr,desiredChannelOrder);
    reordered_EEG.setname = fname;
    reordered_EEG.data = reordered_record;
    save(fname, 'reordered_EEG', '-v7.3');
else
    reordered_hdr = New_Header(hdr, sortedChannelsIndices, desiredChannelOrder);
    save(fname,"reordered_record","reordered_hdr", '-v7.3');
end


fprintf('Successfully saved: %s\n', fname);


end


%% subfunctions


function [reordered_hdr]= New_Header(header,sortedIndices,desired_channel_order)
reordered_hdr = header;
reordered_hdr.label = desired_channel_order;
reordered_hdr.transducer = header.transducer (sortedIndices);
reordered_hdr.units = header.units(sortedIndices);
reordered_hdr.physicalMin = header.physicalMin(sortedIndices);
reordered_hdr.physicalMax = header.physicalMax (sortedIndices);
reordered_hdr.digitalMin = header.digitalMin (sortedIndices);
reordered_hdr.digitalMax = header.digitalMax (sortedIndices);
reordered_hdr.prefilter = header.prefilter (sortedIndices);
reordered_hdr.samples = header.samples (sortedIndices);
reordered_hdr.frequency = header.frequency (sortedIndices);
reordered_hdr.ns = length(sortedIndices) ;
% reordered_hdr. = header. (sortedIndices);
end

function [EEG]= Make_EEG_Header(header,desiredChannelOrder)
load EEG_struct.mat;
EEG.chanlocs.labels = char(desiredChannelOrder);
EEG.srate = unique(header.frequency);
EEG.nbchan = length(desiredChannelOrder);

end


function [folderNames] = Find_Folders(ParentPath)
j = 1;
folderNames = {};
allFilesAndFolders = dir(ParentPath);
for i = 1: length(allFilesAndFolders)

    % skipping non-folders and folders starting with '.' (hidden folders on Unix-like systems)
    foldername = allFilesAndFolders(i).name;

    if allFilesAndFolders(i).isdir &  ~((strcmpi(foldername, '..')) || (strcmpi(foldername, '.')))

        % Get the folder name
        folderNames{j} = allFilesAndFolders(i).name;
        j = j+1;


    end
end
end


function Process_Folder(ParentPath, patientFolder, folderType, desiredSamplingRate, ...
    desiredChannelOrder, removableChannels, channelsToBeReplaced, ...
    newChannelNames, headerFormat, channInfoName)
% PROCESSFOLDER Processes EDF files in a specific folder
%
% This helper function processes all EDF files in a specified folder (diagnosis or follow-up)
% and converts them to MAT format with the specified parameters.

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
     % current_path = pwd;

    % Check if output files already exist
    matFileName1 = strrep(edf_name, '.edf', '_reordered.mat');
    matFileName2 = strrep(edf_name, '.edf', '_reordered_resampled.mat');
    matFilePath1 = fullfile(folder_path, matFileName1);
    matFilePath2 = fullfile(folder_path, matFileName2);

     % matFilePath1 = fullfile(current_path, matFileName1);
     % matFilePath2 = fullfile(current_path, matFileName2);
    
    if exist(matFilePath1, 'file') || exist(matFilePath2, 'file')
        fprintf('Skipping %s - output files already exist\n', edf_name);
        continue;
           % Convert the file
    end
    
        Convert_EDF2Mat(folder_path, edf_name, desiredSamplingRate, ...
            desiredChannelOrder, removableChannels, channelsToBeReplaced, ...
            newChannelNames, headerFormat, ParentPath, channInfoName);
   
    
end
end



function fname = Get_Output_Filename(file_name, original_Fs, desired_sampling_rate)
    % GETOUTPUTFILENAME Determines the output filename based on resampling
    if original_Fs ~= desired_sampling_rate
        fname = strcat(file_name(1:end-4), '_reordered_resampled.mat');
    else
        fname = strcat(file_name(1:end-4), '_reordered.mat');
    end
end

function resampledEEG = Resample_Data(EEG_signal,desired_sampling_rate,original_Fs,num_channels)
    % RESAMPLEDATA Resamples the EEG data to the desired sampling rate
    %
    % Inputs:
    %   EEG_signal - Original EEG data matrix (channels x time points)
    %   original_Fs - Original sampling rate
    %   desired_sampling_rate - Target sampling rate
    %   num_channels - Number of channels
    %
    % Outputs:
    %   resampledEEG - Resampled EEG data matrix with length adjusted for new sampling rate
    
    %% Resample first channel to get the exact output length
    temp_resampled = resample(EEG_signal(1,:), desired_sampling_rate, original_Fs);
    new_length = length(temp_resampled);
        
    % Pre-allocate output array with the exact length from resample
    resampledEEG = zeros(num_channels, new_length);
    
    % Store the first channel's resampled data
    resampledEEG(1,:) = temp_resampled;
    
    % Resample remaining channels
    for i = 2:num_channels
        resampledEEG(i,:) = resample(EEG_signal(i,:), desired_sampling_rate, original_Fs);
    end
end

