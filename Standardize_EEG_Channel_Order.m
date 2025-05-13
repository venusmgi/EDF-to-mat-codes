
% 10.23.2024 by Venus
% header_format
% if you want the output .mat file to have the EEG_lab structure, yo should
% define it as 'EEG_lab', if not, it will have the default format defined
% in "edfreadUntilDone.m"
function [reorderd_record,desiredEEGChannelsIdx] = Standardize_EEG_Channel_Order (desiredChannelOrder,...
    removableChannels,channelsToBeReplaced, newChannelNames, ...
    signalEEG,originalChannelOrder,pathToSaveReport,channelInfoFileName,edfFileName)

% Previoulsy was called Get_desiredChannelOrder_excelOutput_replacing_input
% STANDARDIZE_EEG_CHANNEL_ORDER Reorders and standardizes EEG channel names and data
%
% This function performs three main operations:
% 1. Replaces specified channels with new names
% 2. Reorders channels according to a desired order
% 3. Saves channel information and changes to an .mat file
%
% Inputs:
%   desiredChannelOrder - Cell array of channel names in the desired order
%   removableChannels - Cell array of channel names to be marked as 'unusable'
%   channelsToBeReplaced - Cell array of channel names to be replaced
%   newChannelNames - Cell array of new names for channels to be replaced
%   signalEEG - Numeric matrix of EEG data (channels x time points)
%   originalChannelOrder - Cell array of original channel names
%   pathToSaveReport - String path where output files will be saved
%   channelInfoFileName - String name for the output .mat file
%   (suggestion, use different names each time you run this function, and
%   then concatante the results later)
%   edfFileName - String name of the input EDF file
%
% Outputs:
%   reorderd_record - Reordered EEG data matrix
%   desiredEEGChannels - Indices of channels in the desired order
%
% Example:
%   [reordered_data, channel_indices] = Standardize_EEG_Channel_Order(...
%       {'Fp1','Fp2','F3','F4'}, {'CO2','SpO2'}, {'EEGX3' ,'EEGX4'}, {'EKG1','EKG2'}, ...
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
if ~ischar(pathToSaveReport) || ~ischar(channelInfoFileName) || ~ischar(edfFileName)
    error('pathToSaveReport, channelInfoFileName, and edfFileName must be character arrays');
end



% Verify that the number of channels to be replaced matches the number of new names
assert(length(channelsToBeReplaced)==length(newChannelNames),...
    'Number of channels to be replaced (%d) must match number of new names (%d)',...
    length(channelsToBeReplaced), length(newChannelNames));

% Create a copy of original channel labels for reference
originalChannelLabels = originalChannelOrder;

%% Step1 Removing undesired channele, meaning that we will replace them with "unusable"
originalChannelLabels = Remove_Channel_Names (originalChannelLabels,removableChannels);



%% Step2 : Renaming the unconvential channel names with the mapped channels or inputted channle names or both

% Replace incorrectly named channels with their correct names
originalChannelLabels = Rename_Channel_Labels(originalChannelLabels,channelsToBeReplaced,newChannelNames);


%Counting duplicate channels (e.g. EKG channles will be replaced by  with EKG1, EKG2)
originalChannelLabels = Resolve_Channel_Duplicate(originalChannelLabels);



%% Find Indices and Save Channel Information

% Find indices of desired channels and save updated channel labels

[includedChannelLabels,includedChannelIndices]=Find_Desired_Channels_Order_And_Indices(originalChannelLabels,desiredChannelOrder,removableChannels);

%saving the changes that has been made to channel labels, from renaming to
%re-ordering
Save_Channel_Changes_Info(originalChannelOrder,desiredChannelOrder,edfFileName,includedChannelIndices,includedChannelLabels,pathToSaveReport,channelInfoFileName)



%% Reorder EEG Record
% Reorder the EEG data according to the desired channel order

 recordDesiredChanns = signalEEG(includedChannelIndices,:);


 %if there is any channels that have not been included, find them an
 %include them
[~, desiredEEGChannelsIdx] = ismember (desiredChannelOrder,includedChannelLabels);
if any(desiredEEGChannelsIdx ==0)
    desiredEEGChannelsIdx =  Find_Missed_Channels(includedChannelLabels,desiredChannelOrder);
end



% Save the reordered EEG record
reorderd_record = recordDesiredChanns(desiredEEGChannelsIdx,:);
end



%% subfunctions

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



function Save_Channel_Changes_Info(originalChannelOrder,desiredChannelOrder,edfFileName,includedChannelIndices,includedChannelLabels,pathToSaveReport,channelInfoFileName)
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

% Define path to save the report as a .mat file
reportPath = fullfile(pathToSaveReport, strcat(channelInfoFileName, '.mat'));

% Load existing channel information if the file exists
if isfile(reportPath)


    % Load existing data
    existingData = load(reportPath);
    existingChannelInfo = existingData.channelInfo;

    previousEDFnames = existingChannelInfo(:,1);


    currentEDFname = channelInfo(1,1);

    % finding if this current edf has already been processed and the info
    % has been saved inthe ecel sheet. If this info is already availabe
    % in the report mat file, remove that edf and an other edf after that,
    % because the rest of the files will be processed again

    % TODO : right now, if there is a file, with the same name as it is in
    % the input, it will load that file, and if it does not find a matching
    % name, it will keep everything, but maybe, some EDFs have been skipped
    % before, and now that we are processing them, they will be done again
    % (imagin that the report mat file has only info on the 3rd edf in a file,
    % and we start the process from begining, at fisrt, because the fisrt 
    % edf is beign processed, it cannot find the name of the 3rd one, and
    % does the process, until it reaches the third one, then it will remove
    % all the process, because it found the name of the third, and removed
    % everything after, which might be the first or second one.
    EDFtoRemove =  find(strcmp(previousEDFnames, currentEDFname), 1, 'first');
    if ~isempty(EDFtoRemove)
        existingChannelInfo(EDFtoRemove:end,:) =[];
    end
    channelInfo = Append_New_ChannInfo (channelInfo, existingChannelInfo);

end
   

    % Write the data to a new report mat file
    save(reportPath,'channelInfo');
end


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