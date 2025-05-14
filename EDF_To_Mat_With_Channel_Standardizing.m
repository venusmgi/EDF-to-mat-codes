function EDF_To_Mat_With_Channel_Standardizing (filePath,fileName,desiredSamplingRate,desiredChannelOrder,removableChannels,channelsToBeReplaced, ...
    newChannelNames,pathToSaveChannelInfo,channelReportFileName,removeChannels, renameChannels,headerFormat)
% Previously called CONVERT_EDF2MAT
% This function calls other functions named "edfreadUntilDone.m",
% "Remove_Channel_Names.m", "Rename_Channel_Labels.m",
% "Resolve_Channel_Duplicate.m", and "Save_Channel_Changes_Info.m".
%
% This function reads an EDF file, processes it according to the specified parameters,
% and saves the results in MAT format.
%
% Inputs:
%   file_path - Path to the folder containing the EDF file
%   file_name - Name of the EDF file
%   desired_sampling_rate - Target sampling rate
%   desired_channel_order - Cell array of channel names in desired order
%   removable_channels - Cell array of channels to be excluded
%   channel_tobe_replaced - Cell array of channel names to be replaced
%   new_channel_names - Cell array of new names for channels to be replaced
%   header_format - 'EEGlab' for EEGlab format, empty for default format
%   chann_info_name - Name for the channel information output file
%   removeChannels - 'remove' if there is any spesific channel you want to 
%   remove, input the channel names that you want to remove in 
%   "removableChannels" If there is no spesific channel, and you input 
%   'remove' it will automatically remove the channels that are pre-defined 
%   in "Remove_Channel_Names.m" function.
%   renameChannels - 'rename' if there is any spesific channel you want to 
%   rename and input the channel names that you want to rename in 
%   "channelsToBeReplaced" and their new names in "newChannelNames". 
%   If there is no spesific channel that you wat to rename and you input
%  'rename' it will automatically rename the channels that are pre-defined 
%   in "Rename_Channel_Labels.m " function.



% Add file path to MATLAB path
addpath(filePath)
% Determine if EEGlab format is requested
isEEGlabFormat = strcmpi(headerFormat, 'EEGlab');
% Get number of channels and read EDF file
numChannels = length(desiredChannelOrder);

%% Step1: Convert .edf to .mat
%convert edf to .mat
[hdr, signalEEG] = edfreadUntilDone(string(fileName));
%saving header info
originalChannelOrder = hdr.label;
% Create a copy of original channel labels for reference
originalChannelLabels = originalChannelOrder;
% finding the original sampling rate as it is in the header
originalFs = int64(hdr.frequency(1));  %TODO: figure out a better way to find the original fs

%% Step2: Removing undesired channele, meaning that we will replace them with "unusable"
if strcmpi(removeChannels, 'remove')
    originalChannelLabels = Remove_Channel_Names (originalChannelLabels,removableChannels);
end
%% Step3 : Renaming the unconvential channel names with the mapped channels or inputted channle names or both

% Replace incorrectly named channels with their correct names
if strcmpi(renameChannels, 'rename')
    originalChannelLabels = Rename_Channel_Labels(originalChannelLabels,channelsToBeReplaced,newChannelNames);

    %Counting duplicate channels (e.g. EKG channles will be replaced by  with EKG1, EKG2)
    originalChannelLabels = Resolve_Channel_Duplicate(originalChannelLabels);
end
%% Step 5: Find Indices and Save Channel Information

% Find what indices in that originalChannelOrder are associated with the
% desiredChannelOrder
[includedChannelLabels,includedChannelIndices]=Find_Desired_Channels_Order_And_Indices(originalChannelLabels,desiredChannelOrder,removableChannels);

% saving the changes that has been made to channel labels, from renaming to
% re-ordering
Save_Channel_Changes_Info(originalChannelOrder,desiredChannelOrder,fileName,includedChannelIndices,includedChannelLabels,pathToSaveChannelInfo,channelReportFileName)

%% Step 6: Reorder EEG Record
% Keep only the channels that are associated with desiredChannelOrder
% regardless of the order
recordDesiredChanns = signalEEG(includedChannelIndices,:);

%find in what order we should pick our indecies so we can reorder the
% eeg signal based on the desiredChannelOrder. in other words, sort them.
[~, sortedChannelsIndices] = ismember (desiredChannelOrder,includedChannelLabels);

%if there is any channels that have not been included, find them an
%include them
if any(sortedChannelsIndices ==0)
    sortedChannelsIndices =  Find_Missed_Channels(includedChannelLabels,desiredChannelOrder);
end

% Reorder the reordered EEG record
reordered_record = recordDesiredChanns(sortedChannelsIndices,:);



%% step 7: save the reordered .mat file with renamed channels

% Determine output filename based on resampling
fname = Get_Output_Filename(fileName, originalFs, desiredSamplingRate);

% Resample if needed
if originalFs ~= desiredSamplingRate
    reordered_record = Resample_Data(reordered_record,  desiredSamplingRate, originalFs, numChannels);
    hdr.frequency(1:numChannels) = desiredSamplingRate;
end



% Save the processed data
fprintf('Saving processed data...\n');

%based on the selected header format,
if isEEGlabFormat
    %in case we wanted it to be in EEGLab format
    % Saving the EEG record of the channel order we want
    reordered_EEG = Make_EEGLab_Header(hdr,desiredChannelOrder);
    reordered_EEG.setname = fname;
    reordered_EEG.data = reordered_record;
    save(fname, 'reordered_EEG', '-v7.3');
else
    reordered_hdr = Make_Standard_Header(hdr, sortedChannelsIndices, desiredChannelOrder);
    save(fname,"reordered_record","reordered_hdr", '-v7.3');
end


fprintf('Successfully saved: %s\n', fname);

end


%% Helper Functions

% subfunctions for standerdizing the EEG channel lables


function [includedChannelOrder,includedChannelIndices] = Find_Desired_Channels_Order_And_Indices(channelLabels,desiredOrder,removableChannels)

% FIND_DESIRED_CHANNELS_ORDER_AND_INDICES Identifies indices of desired EEG channels.
%
% This function evaluates a list of EEG channel labels to identify those that fit a specified
% desired order and are not listed as removable. It returns the order and indices of these channels.
%
% Inputs:
%   channelLabels - A cell array containing original channel labels (e.g., {'Fp1', 'Fp2', 'O1', ...}).
%   desiredOrder - A cell array specifying the desired order of channels to retain.
%   removableChannels - A cell array of channel names that should be
%   excluded from consideration. (optional)
%
% Outputs:
%   includedChannelOrder - A cell array of channel labels that match the desired order and are retained.
%   includedChannelIndices - A vector of indices indicating the positions of the included channels within the original list.

% Initialize logical vector for included channels.



% Initialize logical vector for included channels
includedChannelIndices = false(1, length(channelLabels));

% checks if the removableChannels has been inputted
if nargin < 3
    % Evaluate each channel label against desired order and removable list
    for i = 1:length(channelLabels)
        containsRequiredChannels = contains(lower(channelLabels{i}) , lower(desiredOrder) );
        doesNotContainRemovableChannels = ~( contains (channelLabels{i} ,removableChannels));
        % Mark channel as included if it is required and not removable
        includedChannelIndices(:,i) =  containsRequiredChannels && doesNotContainRemovableChannels;
    end
else
    % Evaluate each channel label against desired order and mark it as
    % included if it matches
    for i = 1:length(channelLabels)
        includedChannelIndices(:,i) = contains(lower(channelLabels{i}) , lower(desiredOrder) );
    end
end


% Check if the number of included channels matches the desired order length
if sum(includedChannelIndices) ~= length(desiredOrder)
    % Find missing and extra channels
    missingChannels = setdiff(desiredOrder, channelLabels(includedChannelIndices));
    extraChannels = setdiff(channelLabels(includedChannelIndices), desiredOrder);

    % Build detailed error message
    errorMsg = sprintf('Channel mismatch detected:\n');
    if ~isempty(missingChannels)
        errorMsg = [errorMsg sprintf('Missing channels: %s\n', strjoin(missingChannels, ', '))];
    end
    if ~isempty(extraChannels)
        errorMsg = [errorMsg sprintf('Extra channels: %s\n', strjoin(extraChannels, ', '))];
    end
    errorMsg = [errorMsg sprintf('\nPlease check:\n')];
    errorMsg = [errorMsg sprintf('1. If a channel is marked as 1 in includedChannelIndices but not wanted, add it to removableChannels\n')];
    errorMsg = [errorMsg sprintf('2. If a channel is marked as 0 but needed, add it to desiredChannelOrder\n')];
    error(errorMsg);
end

includedChannelIndices = find(includedChannelIndices); %TODO: isn't it better if i output just the 0 and 1s?
includedChannelOrder = channelLabels(includedChannelIndices);

end




function [sortedIndices] = Find_Missed_Channels (eegRemainingChannelLabels,desiredOrder)
% FIND_SORTED_INDICES Identifies indices in EEG channel labels matching a desired order.
%
% This function searches for and identifies the indices of EEG channels
% that match a specified desired order. It allows for flexibility by matching based on substring presence
% (case-insensitive) rather than direct string equality.
%
% Inputs:
%   EEG_remaining_channel_labels - A cell array of EEG channel labels still available for sorting.
%   desiredOrder - A cell array specifying the order of channels desired for analysis or presentation.
%
% Outputs:
%   sortedIndices - A vector of indices indicating the positions within EEG_remaining_channel_labels
%                   of channels that match the desired order.

% Initialize an empty array to hold sorted indices
sortedIndices = [];

% Initialize a counter for placing found indices
k = 1;

% Iterate over each desired channel name in the specified order
for i = 1: numel(desiredOrder)
    % Convert the desired order name to lowercase for case-insensitive comparison
    lowerDesiredOrder = lower(desiredOrder{i});

    % Iterate over each remaining EEG channel label
    for j = 1:numel(eegRemainingChannelLabels)
        originalChanLabel = lower(eegRemainingChannelLabels{j});

        % Check if the current EEG channel label contains the desired substring
        if contains(originalChanLabel,lowerDesiredOrder)
            sortedIndices(1,k) = j; % Record the index of the matching channel
            k = k+1; % Increment the placement counter
        end
    end
end


end



%% subfunctions for resampling and header format

function [reordered_hdr]= Make_Standard_Header(header,sortedIndices,desired_channel_order)
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

function [EEG]= Make_EEGLab_Header(header,desiredChannelOrder)
load EEG_struct.mat;
EEG.chanlocs.labels = char(desiredChannelOrder);
EEG.srate = unique(header.frequency);
EEG.nbchan = length(desiredChannelOrder);

end


function fname = Get_Output_Filename(file_name, original_Fs, desired_sampling_rate)
% GETOUTPUTFILENAME Determines the output filename based on resampling
%
% Inputs:
%   file_name - Original EDF file name
%   original_Fs - Original sampling rate
%   desired_sampling_rate - Target sampling rate
%
% Outputs:
%   fname - Output filename with appropriate suffix
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

