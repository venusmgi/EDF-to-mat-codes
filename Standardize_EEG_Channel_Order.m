
% 10.23.2024 by Venus
% header_format
% if you want the output .mat file to have the EEG_lab structure, yo should
% define it as 'EEG_lab', if not, it will have the default format defined
% in "edfreadUntilDone.m"
function [reorderd_record,desiredEEGChannels] = Standardize_EEG_Channel_Order (desiredChannelOrder,removableChannels,channelsToBeReplaced, newChannelNames, ...
    signalEEG,originalChannelOrder,parentPath,channelInfoFileName,edfFileName)

% Previoulsy was called Get_desiredChannelOrder_excelOutput_replacing_input
% STANDARDIZE_EEG_CHANNEL_ORDER Reorders and standardizes EEG channel names and data
%
% This function performs three main operations:
% 1. Replaces specified channels with new names
% 2. Reorders channels according to a desired order
% 3. Saves channel information and changes to an Excel file
%
% Inputs:
%   desiredChannelOrder - Cell array of channel names in the desired order
%   removableChannels - Cell array of channel names to be marked as 'unusable'
%   channelsToBeReplaced - Cell array of channel names to be replaced
%   newChannelNames - Cell array of new names for channels to be replaced
%   signalEEG - Numeric matrix of EEG data (channels x time points)
%   originalChannelOrder - Cell array of original channel names
%   parentPath - String path where output files will be saved
%   channelInfoFileName - String name for the output Excel file
%   edfFileName - String name of the input EDF file
%
% Outputs:
%   reorderd_record - Reordered EEG data matrix
%   desiredEEGChannels - Indices of channels in the desired order
%
% Example:
%   [reordered_data, channel_indices] = Standardize_EEG_Channel_Order(...
%       {'Fp1','Fp2','F3','F4'}, {'X1','X2'}, {'X1','X2'}, {'EKG1','EKG2'}, ...
%       eeg_data, original_channels, 'C:\Data', 'channel_info', 'subject1.edf');

% Note: The order and length of channels in channelsToBeReplaced and
% newChannelNames should be the same



% Input validation
if ~iscell(desiredChannelOrder) || isempty(desiredChannelOrder)
    error('desiredChannelOrder must be a non-empty cell array');
end
if ~iscell(removableChannels)
    error('removableChannels must be a cell array');
end
if ~iscell(channelsToBeReplaced) || ~iscell(newChannelNames)
    error('channelsToBeReplaced and newChannelNames must be cell arrays');
end
if ~isnumeric(signalEEG) || ~ismatrix(signalEEG)
    error('signalEEG must be a numeric matrix');
end
if ~iscell(originalChannelOrder) || isempty(originalChannelOrder)
    error('originalChannelOrder must be a non-empty cell array');
end
if ~ischar(parentPath) || ~ischar(channelInfoFileName) || ~ischar(edfFileName)
    error('parentPath, channelInfoFileName, and edfFileName must be character arrays');
end




% Verify that the number of channels to be replaced matches the number of new names
assert(length(channelsToBeReplaced)==length(newChannelNames),...
    'Number of channels to be replaced (%d) must match number of new names (%d)',...
    length(channelsToBeReplaced), length(newChannelNames));

% Create a copy of original channel labels for reference
originalChannelLabels = originalChannelOrder;

% Step 1: Mark removable channels as 'unusable'
% This step identifies channels that should be excluded from the final output
for i = 1:length(removableChannels)
    if any(contains(originalChannelOrder,removableChannels{i}))
        removableChannelIdx = contains(originalChannelOrder,removableChannels{i});
        originalChannelOrder{removableChannelIdx} = 'unusable';
    end
end


%% standardizing the Channel names

% Step 2: Standardize channel names
% Replace incorrectly named channels with their correct names
if ~isempty(channelsToBeReplaced) && ~isempty(newChannelNames)
    for i = 1:length(channelsToBeReplaced)
        chan_name{1} = channelsToBeReplaced{i};
        new_chan_name = newChannelNames{i};
        originalChannelLabels = Replace_Channel_Names(originalChannelLabels, chan_name, new_chan_name);
    end
end


% Define mappings for possible channel names and their desired replacements
% This mapping helps standardize channel names across different recording systems
channelMappings = {'Eye1',{'LUO','EEGPOLLUO','POLEOGL','POLLLC','POLLOC','EEGLOCRef', 'LOC','EOGL','LLC', 'LUE','LUOC','Reye','POLLLE',...
    'LEYE','LIO','EEGLEYERef','LOF','LEOG','POLLLE','POLLUE','LLE','EOGLT'};
    'Eye2',{'RLO','EEGPOLRLO','POLEOGR','POLROC','EEGROCRef','EOGR','ROC','RAE','RUE','RLOC','Leye','POLRUE','REYE','RIO',...
    'EEGREYERef','ROF','REOG','POLRUE','POLRLE','EOGRT'};
    'EKG1',{'ECGL','ECG1','ecg1','EKGL','LEKG','EEGLEKGRef','ECGLA','EEGECGLRef','EKGLT'};
    'EKG2',{'ECGR','ECG2','ecg2','EKGR','REKG','EEGREKGRef','ECGRA','EEGECGRRef','EKGRT'};
    'EKG',{'ECG','EEGEKGRef','POLEKG','EEGPOLEKG'};
    'EMG',{'CHIN1','CHIN2','NECK1','NECK2','NEC1','NEC2','Lleg1','Lleg2','Rleg1','Rleg2','chin','EEGNeckRef','POLNeck1',...
    'POLNeck2','neck1','neck2','LEMG1','REMG1','EEGCHIN1Ref','EEGCHIN2Ref','POLNECK1','POLNECK2','POLChin1','POLChin2','RLEG',...
    'LLEG','EMGR','EMGL','ABD1','ABD2','EEGABD1Ref','EEGABD2Ref','CHINLT','CHINRT','ABDBLK','ABDWHT','UCHIN','LCHIN'};
    'T3',{'T7','EEGT7Ref'};
    'T4',{'T8','EEGT8Ref'};
    'T5',{'P7','EEGP7Ref'} ;
    'T6', {'P8','EEGP8Ref'}
    };

% Apply the channel mappings to standardize channel names
for k = 1:length(channelMappings)
    desiredName = channelMappings{k, 1};
    possibleNames = channelMappings{k, 2};
    originalChannelLabels= Replace_Channel_Names (originalChannelLabels, possibleNames,  desiredName );
end


% Renaming and Counting Channels (EKGs with EKG1 and EKG2, and EMGs with EMG1 and EMG2)

channelCounters = struct('EMG', 0, 'EKG', 0, 'Fz', 0, 'O1', 0, 'O2', 0);
% Initialize a struct to hold channel name counters


% Iterate through the channel names to standardize and count them
for i = 1:length(originalChannelLabels)
    [originalChannelLabels{i},channelCounters] = standardizeChannelName(originalChannelLabels{i}, channelCounters);
end



%% Find Indices and Save Channel Information

% Find indices of desired channels and save updated channel labels

[includedChannelLabels,includedChannelIndices]=Find_Desired_Channels_Order_And_Indices(originalChannelLabels,desiredChannelOrder,removableChannels);
Save_Channel_Changes_Info(originalChannelOrder,desiredChannelOrder,edfFileName,includedChannelIndices,includedChannelLabels,parentPath,channelInfoFileName)



%% Reorder EEG Record
% Reorder the EEG data according to the desired channel order

desiredEEGChannels = signalEEG(includedChannelIndices,:);

[~, desiredEEGChannels] = ismember (desiredChannelOrder,includedChannelLabels);
if any(desiredEEGChannels ==0)
    desiredEEGChannels =  Find_Sorted_indices(includedChannelLabels,desiredChannelOrder);
end



% Save the reordered EEG record
reorderd_record = desiredEEGChannels(desiredEEGChannels,:);
end



%% subfunctions

function [includedChannelOrder,includedChannelIndices] = Find_Desired_Channels_Order_And_Indices(channelLabels,desiredOrder,removableChannels)

% FIND_DESIRED_CHANNELS_ORDER_AND_INDICES Finds indices of desired channels
%
% This function identifies channels that match the desired order and are not removable.
% It returns the order and indices of these channels.




for i = 1:length(channelLabels)
    contains_required_channels = contains(lower(channelLabels{i}) , lower(desiredOrder) );
    does_not_contain_removableChannels = ~( contains (channelLabels{i} ,removableChannels));
    includedChannelIndices(:,i) =  contains_required_channels && does_not_contain_removableChannels;
end




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
% contains(lower(channelLabels(1,i)) , lower(desiredOrder) );
includedChannelIndices=find(includedChannelIndices);

includedChannelOrder = channelLabels(includedChannelIndices);

end




function [sortedIndices] = Find_Sorted_indices (EEG_remaining_channel_labels,desiredOrder)
% FIND_SORTED_INDICES Finds sorted indices of channels based on desired order
%
% This function finds the indices of channels in the desired order, even if the direct match is not found.

sortedIndices = []; %sometimes the ismemeber does not work, so we go through this loop to find the indices based on our desired order
k = 1;
for i = 1: numel(desiredOrder)
    lower_desiredOrder = lower(desiredOrder{i});
    for j = 1:numel(EEG_remaining_channel_labels)
        original_chan_label = lower(EEG_remaining_channel_labels{j});
        if contains(original_chan_label,lower_desiredOrder)
            sortedIndices(1,k) = j;
            k = k+1; %%for checking contanins put a pause excecution here
        end
    end
end


end



function [outputChannelNames] = Replace_Channel_Names (originalChannelNames, channelsToBeReplace, replacement_channel)
% REPLACE_CHANNEL_NAMES Replaces specified channel names with new names
%
% This function replaces any occurrence of specified channels with a new name.

for j = 1:length(channelsToBeReplace)

    % checks if any of the possible other names of the channels that
    % should be replaced exists in the channel names
    if any(strcmp(originalChannelNames, channelsToBeReplace{j}))
        % replaces that channels with the defined  channel name
        if sum(strcmp(originalChannelNames, channelsToBeReplace{j})) ==1
            originalChannelNames{strcmp(originalChannelNames, channelsToBeReplace{j})} = replacement_channel;
        else

            IdxToBeReplaced =  find(strcmp(originalChannelNames, channelsToBeReplace{j}));
            for i = 1:length(IdxToBeReplaced)
                originalChannelNames{IdxToBeReplaced(i)} = replacement_channel;
            end
        end
    end

end

outputChannelNames = originalChannelNames;

end

function [newLabel,channelCounters] = standardizeChannelName(channelLabel, channelCounters)

% STANDARDIZECHANNELNAME Renames EEG channels to ensure standardized naming conventions.
%
% This function adjusts channel names to add numerical suffixes to duplicate channels,
% ensuring that each channel has a unique name. For main EEG channels like 'O1', 'O2', and 'Fz',
% the base name is preserved for the first occurrence without adding a numerical suffix.
%
% Inputs:
%   channelLabel - The original label of the EEG channel.
%   channelCounters - A struct containing counters for each channel type
%                     (e.g., 'EMG',0, 'EKG',0, 'Fz', 1, 'O1', 1, 'O2', 0),
%                       which tracks the number of occurrences of each channel type.
%
% Outputs:
%   newLabel - The standardized channel label. If the channel is a duplicate,
%              it includes a number suffix (e.g., 'EMG1', 'EMG2').
%   channelCounters - Updated struct with the incremented counters for each channel type.



% Define the channel names of interest
channelNames = {'EMG', 'EKG', 'Fz', 'O1', 'O2'};

% Loop through each channel name to check if the current label matches any known types
for i = 1:length(channelNames)
    channelName = channelNames{i};

    if contains(channelLabel,channelName)
        % Update the counter for the channel
        channelCounters.(channelName) = channelCounters.(channelName) + 1;
        % For the main EEG channels ('O1', 'O2', 'Fz'), preserve the
        % base name without adding a
        % number for the first encounter, as they are standard electrode names
        if (strcmp(channelName,'O1') || strcmp(channelName,'O2') ||strcmp(channelName,'Fz')) && channelCounters.(channelName) == 1
            newLabel = channelName; % Preserve the base name without number for the first instance
        else
            % Append a numeral suffix for duplicate or non-main channels
            newLabel = [channelName num2str(channelCounters.(channelName))];
        end
        return;
    end
end
% Return the unchanged label if no matching channel name is found
newLabel = channelLabel;



end



function paddedArray = Pad_Cell_Array(cellArray, numRows, numColumns)
% PADCELLARRAY Pads a cell array with empty cells to match desired dimensions
%
% Inputs:
%   cellArray - The cell array to pad
%   numRows - Number of rows in the padded array
%   numColumns - Number of columns in the padded array
%
% Outputs:
%   paddedArray - The padded cell array

% Create padding of empty cells
padding = cell(numRows, numColumns - size(cellArray, 2));

% Combine original array with padding
paddedArray = [cellArray, padding];
end
% Ensure the new and existing data have the same number of columns

function newChannelInfo = Append_New_ChannInfo (newChannelInfo, existingChannelInfo)

numColumnsChannInfo = size(newChannelInfo,2);
numColumnsExistingChans = size(existingChannelInfo,2);
%checks to see if the size of the new channel info is the same
%as the uploaded on, if not, padds to be the same size

if numColumnsExistingChans > numColumnsChannInfo
    % Pad channelInfo with empty cells to match existingChannelInfo
    newChannelInfo = Pad_Cell_Array(newChannelInfo, size(newChannelInfo, 1), numColumnsExistingChans);
elseif numColumnsExistingChans < numColumnsChannInfo
    % Pad existingChannelInfo with empty cells to match channelInfo
    existingChannelInfo = Pad_Cell_Array(existingChannelInfo, size(existingChannelInfo, 1), numColumnsChannInfo);
end


% Append the new data
newChannelInfo = [existingChannelInfo; newChannelInfo];
end



function Save_Channel_Changes_Info(originalChannelOrder,desiredChannelOrder,edfFileName,includedChannelIndices,includedChannelLabels,parentPath,channelInfoFileName)
% Create a cell array to store channel information
numChannels = length(originalChannelOrder);
channelInfo = cell(3,1 + numChannels); % 3 rows: Original Name, New Name, Reordered Name

% Set row headers
channelInfo(:,1) = {edfFileName;
    strcat('New',' ',edfFileName);
    strcat('Reordered',' ',edfFileName)
    };



% Fill in channel info
for i = 1:length(originalChannelOrder)
    channelInfo{ 1,i + 1} = originalChannelOrder{i};

    % Check if the channel was renamed
    if any(includedChannelIndices == i)
        channelInfo{2, i + 1} = includedChannelLabels{includedChannelIndices == i};
        channelInfo{3, i+1}= desiredChannelOrder{includedChannelIndices == i};


    else % skip this channel
        channelInfo{2, i + 1} = '';
        channelInfo{3, i+1}= '';
    end
end

% Define path to Excel sheet as a .mat file
excelSheetPath = fullfile(parentPath, strcat(channelInfoFileName, '.mat'));

% Load existing channel information if the file exists
if isfile(excelSheetPath)


    % Load existing data
    existingData = load(excelSheetPath);
    existingChannelInfo = existingData.channelInfo;

    previousEDFnames = existingChannelInfo(:,1);


    currentEDFname = channelInfo(1,1);

    % finding if this current edf has already been processed and the info
    % has been saved inthe ecel sheet. If this info is already availabe
    % in the excel sheet, remove that edf and an other edf after that,
    % because the rest of the files will be processed again

    EDFtoRemove =  find(strcmp(previousEDFnames, currentEDFname), 1, 'first');
    if ~isempty(EDFtoRemove)
        existingChannelInfo(EDFtoRemove:end,:) =[];
    end


    channelInfo = Append_New_ChannInfo (channelInfo, existingChannelInfo);


    % Write the data to a new Excel file
    save(excelSheetPath,'channelInfo');
end
end