% header_format
% if you want the output .mat file to have the EEG_lab structure, yo should
% define it as 'EEG_lab', if not, it will have the default format defined
% in "edfreadUntilDone.m"

function Convert_EDFtomat_Desired_Fs_channOrder_Parent_Path_ExcOut_repCh (desired_sampling_rate,desired_channel_order,removable_channels, channel_tobe_replaced, new_channel_names,chann_info_name,header_format)
% CONVERT_EDF2MAT Converts a single EDF file to MAT format with specified parameters
%
% This function converts an EDF file to MAT format, reordering channels and
% resampling the data according to specified parameters.
%
% Inputs:
%   file_path - Path to the EDF file
%   file_name - Name of the EDF file
%   desired_sampling_rate - Target sampling rate for the output
%   desired_channel_order - Cell array of channel names in desired order
%   removable_channels - Cell array of channels to be excluded
%   channel_tobe_replaced - Cell array of channel names to be replaced
%   new_channel_names - Cell array of new names for channels to be replaced
%   header_format - 'EEGlab' for EEGlab format, empty for default format
%   ParentPath - Parent directory path
%   chann_info_name - Name for the channel information output file
%   edf_name - Name of the EDF file being processed
%
% Outputs:
%   Saves a .mat file with either EEGlab format or standard format

functionPath = pwd;
addpath(functionPath)

if nargin < 7
    header_format = '';  % Default value (regular header format)
end

ParentPath = uigetdir;
patient_folders = Find_folders(ParentPath);

% Process each patient folder
for i = 1:length(patient_folders)
    fprintf('Processing patient folder: %s\n', patient_folders{i});
    
    % Process diagnosis folder
    ProcessFolder(ParentPath, patient_folders{i}, 'diagnosis', desired_sampling_rate, ...
        desired_channel_order, removable_channels, channel_tobe_replaced, ...
        new_channel_names, header_format, chann_info_name);
    
    % Process follow-up folder
    ProcessFolder(ParentPath, patient_folders{i}, 'follow up', desired_sampling_rate, ...
        desired_channel_order, removable_channels, channel_tobe_replaced, ...
        new_channel_names, header_format, chann_info_name);
end 

fprintf('All files processed successfully.\n');
   
end





function Convert_EDF2Mat (file_path,file_name,desired_sampling_rate,desired_channel_order,removable_channels,channel_tobe_replaced, new_channel_names,header_format,ParentPath,chann_info_name,edf_name)
% CONVERT_EDF2MAT Converts an EDF file to MAT format with specified parameters
%
% This function reads an EDF file, processes it according to the specified parameters,
% and saves the results in MAT format.
%
% Inputs:
%   folder_path - Path to the folder containing the EDF file
%   edf_name - Name of the EDF file
%   desired_sampling_rate - Target sampling rate
%   desired_channel_order - Cell array of channel names in desired order
%   removable_channels - Cell array of channels to be excluded
%   channel_tobe_replaced - Cell array of channel names to be replaced
%   new_channel_names - Cell array of new names for channels to be replaced
%   header_format - 'EEGlab' for EEGlab format, empty for default format
%   ParentPath - Parent directory path
%   chann_info_name - Name for the channel information output file
%   edf_name - Name of the EDF file (for logging)


addpath(file_path)
isEEGlabFormat = strcmpi(header_format, 'EEGlab');
num_channels = length(desired_channel_order);
[hdr, record] = edfreadUntilDone(string(file_name));
EEG_original_channel_order = hdr.label;
%     True_Fs = hdr.samples(1)/ hdr.duration;
original_Fs = int64(hdr.frequency(1));  %TODO: figure out a better way to find the original fs

EEG_record = record;
clear record

%% removing the channels that have similarity in names with our main 19 channels


[reordered_record,sorted_channels_Indices] =Get_desired_channel_order_excelOutput_replacing_input (desired_channel_order,removable_channels,channel_tobe_replaced, ...
    new_channel_names,EEG_record, EEG_original_channel_order,ParentPath,chann_info_name,edf_name);



% Determine output filename based on resampling
fname = Get_Output_Filename(file_name, original_Fs, desired_sampling_rate);
reordered_EEG.setname = fname;


% Resample if needed
if original_Fs ~= desired_sampling_rate
    reordered_record = Resample_Data(reordered_record,  desired_sampling_rate, original_Fs, num_channels);
    hdr.frequency(1:num_channels) = desired_sampling_rate;
end

%in case we wanted it to be in EEGLab format
% Saving the EEG record of the channel order we want
reordered_EEG = make_EEG_header(hdr,desired_channel_order);
reordered_EEG.data = reordered_record;



% Save the processed data
fprintf('Saving processed data...\n');
Save_Processed_Data(fname, isEEGlabFormat, reordered_EEG, hdr, ...
    sorted_channels_Indices, desired_channel_order, reordered_record);

fprintf('Successfully saved: %s\n', fname);


clear EEG_original_channel_order EEG_record reordered_record hdr reordered_hdr reordered_EEG fname file_name file_name

end
%%subfunctions
% EEG.chanlocs.labels

function [reordered_hdr]= new_header(header,sortedIndices,desired_channel_order)
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

function [EEG]= make_EEG_header(header,desired_channel_order)
load EEG_struct.mat;
EEG.chanlocs.labels = char(desired_channel_order);
EEG.srate = unique(header.frequency);
EEG.nbchan = length(desired_channel_order);

end



function [folderNames] = Find_folders(ParentPath)
j = 1;
folderNames = {}
All_files_and_folders = dir(ParentPath);
for i = 1: length(All_files_and_folders)

    % skipping non-folders and folders starting with '.' (hidden folders on Unix-like systems)
    foldername = All_files_and_folders(i).name;

    if All_files_and_folders(i).isdir &  ~((strcmpi(foldername, '..')) || (strcmpi(foldername, '.')))

        % Get the folder name
        folderNames{j} = All_files_and_folders(i).name;
        j = j+1;


    end
end
end


function ProcessFolder(ParentPath, patient_folder, folder_type, desired_sampling_rate, ...
    desired_channel_order, removable_channels, channel_tobe_replaced, ...
    new_channel_names, header_format, chann_info_name)
% PROCESSFOLDER Processes EDF files in a specific folder
%
% This helper function processes all EDF files in a specified folder (diagnosis or follow-up)
% and converts them to MAT format with the specified parameters.

% Construct folder path
folder_path = fullfile(ParentPath, patient_folder, folder_type);
if ~exist(folder_path, 'dir')
    fprintf('Warning: %s folder not found for patient %s\n', folder_type, patient_folder);
    return;
end

% Change to folder and get EDF files
cd(folder_path);
edf_files = dir('*.edf');
if isempty(edf_files)
    fprintf('No EDF files found in %s folder for patient %s\n', folder_type, patient_folder);
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
    try
        Convert_EDF2Mat(folder_path, edf_name, desired_sampling_rate, ...
            desired_channel_order, removable_channels, channel_tobe_replaced, ...
            new_channel_names, header_format, ParentPath, chann_info_name, edf_name);
    catch ME
        fprintf('Error processing %s: %s\n', edf_name, ME.message);
    end
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

function resmapledEEG = Resample_Data(EEG_signal,desired_sampling_rate,original_Fs,num_channels)
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

function Save_Processed_Data(fname,isEEGlabFormat,reordered_EEG,hdr,...
                              sorted_channels_Indices,desired_channel_order,reordered_record)
    % Saving the EEG record of the channel order we want
    if isEEGlabFormat
        reordered_EEG = make_EEG_header_EEGlab(hdr, desired_channel_order);
    
        save(fname, 'reordered_EEG', '-v7.3');
    else
        reordered_hdr = new_header(hdr, sorted_channels_Indices, desired_channel_order);
        save(fname,"reordered_record","reordered_hdr", '-v7.3');
    end


end