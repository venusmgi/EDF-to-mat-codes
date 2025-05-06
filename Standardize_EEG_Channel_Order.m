
% 10.23.2024 by Venus
% header_format
% if you want the output .mat file to have the EEG_lab structure, yo should
% define it as 'EEG_lab', if not, it will have the default format defined
% in "edfreadUntilDone.m"
function [reorderd_record,desiredEEGChannels] = Standardize_EEG_Channel_Order (desiredChannelOrder,removableChannels,channelsToBeReplaced, newChannelNames, ...
    signalEEG,originalChannelOrder,parentPath,channelInfoFileName,edfFileName)

% Previoulsy was called Get_desiredChannelOrder_excelOutput_replacing_input

% This function reorders EEG channels based on a desired order, replaces specified channel names,
% and generates a report of the changes. It handles standard EEG, EKG, EMG, and ocular channels.
%
% Inputs:
%   desiredChannelOrder - Desired order of EEG channels
%   removableChannels - Channels to be excluded from the final order
%   channelsToBeReplaced - Channels with incorrect names to be replaced(
%   e.g. a channel that is named X1, but we know that it is EKG)
%   newChannelNames- Correct names for the channels to be replaced (EKG
%   for X1)
%   signalEEG - EEG data matrix (channels x time points)
%   originalChannelOrder - Original order of EEG channel names
%   parentPath - Path to save the channel information
%   channelInfoFileName - Name for the channel information file to be saved
%   
%   edfFileName - Name of the EDF file being processed
%
% Outputs:
%   reorderd_record - EEG data reordered according to the desired channel order
%   desiredEEGChannels - Indices of the sorted channels in the original data

% Note : the order and length of channels in channelsToBeReplaced and
% newChannelNamesshould should be the same


originalChannelLabels= originalChannelOrder;  % Copy original channel labels

%replacing any of the channels that are in the removableChannels, with
%'unusable'
for i = 1:length(removableChannels)
    if any(contains(originalChannelOrder,removableChannels{i}))
        removableChannelIdx = contains(originalChannelOrder,removableChannels{i});
        originalChannelOrder{removableChannelIdx} = 'unusable';
    end
end

% Assert that the number of channels to be replaced matches the new names
assert (length(channelsToBeReplaced)==length(newChannelNames),'number of channels that are being replaced and their new names should be the same')
%% standardizing the Channel names

% Replace incorrectly named channels with new names
if ~isempty(channelsToBeReplaced) && ~isempty(newChannelNames)
    for i = 1:length(channelsToBeReplaced)
        chan_name{1} = channelsToBeReplaced{i};
        new_chan_name = newChannelNames{i};
        originalChannelLabels= Replace_Channel_Names (originalChannelLabels, chan_name,  new_chan_name );

    end
end


% Define mappings for possible channel names and their desired replacements
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

%applying the channel maps
for k = 1:length(channelMappings)
        desiredName = channelMappings{k, 1};
        possibleNames = channelMappings{k, 2};
        originalChannelLabels= Replace_Channel_Names (originalChannelLabels, possibleNames,  desiredName );
end



%% %% Renaming and Counting Channels (EKGs with EKG1 and EKG2, and EMGs with EMG1 and EMG2)

channelCounters = struct('EMG', 0, 'EKG', 0, 'Fz', 0, 'O1', 0, 'O2', 0);
% Initialize a struct to hold channel name counters


% Iterate through the channel names to standardize and count them
for i = 1:length(originalChannelLabels)
    [originalChannelLabels{i},channelCounters] = standardizeChannelName(originalChannelLabels{i}, channelCounters);
end



%% Find Indices and Save Channel Information

% Find indices of desired channels and save updated channel labels
[includedChannelLabels,includedChannelIndices]=Find_Desired_Channels_Order_And_Indices(originalChannelLabels,desiredChannelOrder,removableChannels);

% Create a cell array to store channel information
channelInfo = cell(3,1 + length(originalChannelLabels)); % 3 columns: Original Name, New Name, Skipped
row1Name = edfFileName;
row2Name = strcat('New',' ',edfFileName);
row3Name = strcat('Reordered',' ',edfFileName);

% Fill in the header
channelInfo(:,1) = {row1Name, row2Name,row3Name};

% Fill in the data
for i = 1:length(originalChannelOrder)
    channelInfo{ 1,i + 1} = originalChannelOrder{i};

    % Check if the channel was renamed
    if any(includedChannelIndices == i)
        channelInfo{2, i + 1} = includedChannelLabels{includedChannelIndices == i};
        channelInfo{3, i+1}= desiredChannelOrder{includedChannelIndices == i};


    else
        channelInfo{2, i + 1} = '';
        channelInfo{3, i+1}= '';
    end
end

% Save the channel information to a .mat file
excelSheetPath = strcat(parentPath,'\',channelInfoFileName,'.mat');
% Check if the file already exists
if exist(excelSheetPath) == 2
    % Load existing data
    existingData = load(excelSheetPath);


    previousEDFnames = existingData.channelInfo(:,1);
    currentEDFname = channelInfo(1,1);

    % finding if this current edf has already been processed and the ibnfo
    % has been saved inthe ecel sheet. If this info is already availabe
    % in the excel sheet, remove that edf and an other edf after that, 
    % because the rest of the files will be processed again

    EDFtoRemove =  find(strcmp(previousEDFnames, currentEDFname), 1, 'first');
    if ~isempty(EDFtoRemove)
        existingData.channelInfo(EDFtoRemove:end,:) =[];
    end

    % Ensure the new and existing data have the same number of columns
    numColumns1 = size(channelInfo,2);
    numColumns2 = size(existingData.channelInfo,2);
    %checks to see if the size of the new channel info is the same
    %as the uploaded on, if not, padds to be the same size
    if numColumns2 > numColumns1
        % Pad cellArray2 with empty cells
        padding = cell(size(channelInfo, 1), numColumns2 - numColumns1);
        channelInfo = [channelInfo, padding];
    elseif numColumns2 < numColumns1
        padding = cell(size(existingData.channelInfo, 1), numColumns1 - numColumns2);
        existingData.channelInfo = [existingData.channelInfo, padding];
    end

    % Append the new data
    channelInfo = [existingData.channelInfo; channelInfo];


    % Write the updated data to the Excel file
    save(excelSheetPath,'channelInfo')
else
    % Write the data to a new Excel file
    save(excelSheetPath,'channelInfo');
end


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



% %this part replaces A1 and A2, with M1 and M2
%     if (any(contains(channelLabels,'A1').* ~strcmp(channelLabels,'POL $A1')) & any(strcmp(desiredOrder,'M1')))
%         channelLabels{  find(contains(channelLabels,'A1').* ~strcmp(channelLabels,'POL $A1'))  } ='M1';
%     elseif (any(contains(channelLabels,'M1').* ~strcmp(channelLabels,'POLM1')) & any(strcmp(desiredOrder,'A1')))
%         channelLabels{  find(contains(channelLabels,'M1').* ~strcmp(channelLabels,'POLM1'))  } ='A1';
%     end
%
%     if (any(contains(channelLabels,'A2').* ~strcmp(channelLabels,'POL$A2')) & any(strcmp(desiredOrder,'M2')))
%         channelLabels{  find(contains(channelLabels,'A2').* ~strcmp(channelLabels,'POL$A2'))  } ='M2';
%     elseif (any(contains(channelLabels,'M2').* ~strcmp(channelLabels,'POLM2')) & any(strcmp(desiredOrder,'A2')))
%         channelLabels{  find(contains(channelLabels,'M2').* ~strcmp(channelLabels,'POLM2'))  } ='A2';
%     end

for i = 1:length(channelLabels)
    contains_required_channels = contains(lower(channelLabels(1,i)) , lower(desiredOrder) );
    does_not_contain_removableChannels = ~( contains (channelLabels(1,i) ,removableChannels));
    includedChannelIndices(:,i) =  contains_required_channels && does_not_contain_removableChannels;
end

% for i = 1:length(channelLabels)
%         contains_required_channels1 = contains(lower(channelLabels{i}) , lower(desiredOrder) );
%         does_not_contain_removableChannels1 = ~( contains (channelLabels{i} ,removableChannels));
%         includedChannelIndices1(:,i) =  contains_required_channels && does_not_contain_removableChannels;
% end


if sum(includedChannelIndices) ~= length(desiredOrder)
    error(strcat('Number of channels is smaller or higher than number of desired channel. check channel labels.',...
        'check variable "includedChannelIndices" here and see what channels are missing from "chan_label" or are being added aditionally.',...
        'if "includedChannelIndices" is 1, it means the corresponding channel in "chan_label"is being included, if you do not want it put it in input "removableChannels".',...
        'if in "includedChannelIndices" the corresponding index to "chan_label" is 0, and you want that channel, add it to input "desiredChannelOrder".',...
        'if it is 1 and you do not want it add the channel name in "chan_label" to "removableChannels"' ))
    return
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
             % For the main EEG channels ('O1', 'O2', 'Fz'), preserve the base name without a
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
