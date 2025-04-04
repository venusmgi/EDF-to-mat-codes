% Convert_EDFtomat_Desired_Fs_channOrder_outputExcel_replaceChans
% This function converts EDF files to .mat files with a desired sampling rate and channel order.
% It also allows for replacing specific channel names and exporting channel information to Excel.
%
% Inputs:
%   desiredSamplingRate - The target sampling rate for resampling the EEG data.
%   desiredChannelOrder - A cell array specifying the desired order of EEG channels.
%   channsToBeRemoved    - Channels to be removed from the EEG data.
%   channsToBeReplaced - Channels with incorrect names to be replaced.
%   channsNewName     - Correct names for the channels to be replaced.
%   channsInfofileName       - Name for the channel information file.
%   headerFormat         - Format of the header ('EEGlab' for EEGlab-compatible format).

% By Venus 10.23.2024
function Convert_EDFtomat_Desired_Fs_channOrder_outputExcel_replaceChans (desiredSamplingRate,desiredChannelOrder,channsToBeRemoved,channsToBeReplaced, channsNewName,channsInfofileName, headerFormat)
% Check if the headerFormat is provided and not empty
if nargin <7
    headerFormat = '';  % Default value (regular header format)
end
isEEGlabFormat = strcmpi(headerFormat, 'EEGlab'); % Check if EEGlab format is requested
nChan= length(desiredChannelOrder); % Number of desired output channels

% Prompt user to select EDF files
[fileName, path] = uigetfile('*.edf',...
    'Select EDF File(s)', ...
    'MultiSelect', 'on');
addpath(path)

% Handle single or multiple file selection
    if ischar(fileName)
        fileName = {fileName}; % Convert to cell array for consistency
    end

% Process each selected EDF file
for i = 1:length(fileName)
        currentFileName = fileName{i};
    % Define output .mat file names
    matFileName1 = strrep(currentFileName, '.edf', '_reordered.mat');
    matFileName2 = strrep(currentFileName, '.edf', '_reordered_resampled.mat');
    matFilePath1 = fullfile(path, matFileName1);
    matFilePath2 = fullfile(path, matFileName2);


    % Check if the .mat files already exist
    if (exist(matFilePath1)==2)|(exist(matFilePath2)==2)

        display(strcat(matFileName1,' already exists'));
   
    else
        % Read the EDF file
        [hdr, record] = edfreadUntilDone(string(currentFileName));

        % Extract original channel order and sampling frequency
        OrigChanOrder = hdr.label;
        origFs = unique(hdr.frequency);




        %% standerdizing of channels process

        % Process the EEG record
        EEG_record = record;
        clear record
        % Reorder channels and replace specified channels
        [reordered_record,sortedChannelIndices] = Standardize_EEG_Channel_Order(desiredChannelOrder, ...
                                                                                channsToBeRemoved,channsToBeReplaced, ...
                                                                                channsNewName,EEG_record, OrigChanOrder, ...
                                                                                path,channsInfofileName,currentFileName);
        % [reordered_record,sortedChannelIndices] =Get_desiredChannelOrder_output_excel (desiredChannelOrder,channsToBeRemoved,EEG_record, OrigChanOrder,path,channsInfofileName,currentFileName);
        % Create a new EEG structure with the desired channel order


        % Check if resampling is needed and in that case change the EEG
        % file name
        if origFs ~= desiredSamplingRate
            % Resample the data if needed
            reorderedRecordOld = reordered_record;
            clear reordered_record
            for i = 1:nChan
                reordered_record(i,:) = resample(reorderedRecordOld(i,:) ,desiredSamplingRate,origFs);
            end

            newFileName = strcat(currentFileName(1:end-4),'_reordered_resampled.mat');

        else
            newFileName = strcat(currentFileName(1:end-4),'_reordered.mat');
        end

        


        % Update header with new sampling rate
        hdr.frequency(1:nChan) = desiredSamplingRate;


        % Save the reordered EEG data
        if isEEGlabFormat
            reorderedEEG = Make_EEGLab_header(hdr,desiredChannelOrder);
            reorderedEEG.data = reordered_record;
            reorderedEEG.setname = newFileName;
            reorderedEEG = make_EEG_header_EEGlab(hdr, desiredChannelOrder);
            save(newFileName, 'reorderedEEG', '-v7.3');
        else
            reordered_hdr = Make_New_Header(hdr, sortedChannelIndices, desiredChannelOrder);
            save(newFileName,"reordered_record","reordered_hdr", '-v7.3');
        end
        display(strcat('done saving ',newFileName))

        clear OrigChanOrder EEG_record reordered_record hdr reorderedRecordOld  reordered_hdr reorderedEEG 
        clear sortedChannelIndices newFileName currentFileName
    end
end
end
%% subfunctions
% EEG.chanlocs.labels

function [reordered_hdr]= Make_New_Header(header,sortedIndices,desiredChannelOrder)
    reordered_hdr = header;
    reordered_hdr.label = desiredChannelOrder;
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

% Create an EEG header compatible with EEGlab
function [EEG]= Make_EEGLab_header(header,desiredChannelOrder)
    load EEG_struct.mat;
    EEG.chanlocs.labels = char(desiredChannelOrder);
    EEG.srate = unique(header.frequency);
    EEG.nbchan = length(desiredChannelOrder);

end