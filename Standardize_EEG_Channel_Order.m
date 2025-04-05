
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



% Assert that the number of channels to be replaced matches the new names
assert (length(channelsToBeReplaced)==length(newChannelNames),'number of channels that are being replaced and their new names should be the same')
%% standardizing the Channel names
originalChannelLabels= originalChannelOrder;  % Copy original channel labels
% Replace incorrectly named channels with new names
if ~isempty(channelsToBeReplaced) && ~isempty(newChannelNames)
    for i = 1:length(channelsToBeReplaced)
        chan_name{1} = channelsToBeReplaced{i};
        new_chan_name = newChannelNames{i};
        originalChannelLabels= Replace_Channel_Names (originalChannelLabels, chan_name,  new_chan_name );

    end
end


% Define possible names for ocular, EKG, and EMG channels and standardize them
OcularChannelsPossibleNames1 = {'LUO','EEGPOLLUO','POLEOGL','POLLLC','POLLOC','EEGLOCRef', 'LOC','EOGL','LLC', 'LUE','LUOC','Reye','POLLLE',...
    'LEYE','LIO','EEGLEYERef','LOF','LEOG','POLLLE','POLLUE','LLE','EOGLT'}; %--> you did not keep thi
desiredOcularChannelName1 = 'Eye1';
originalChannelLabels= Replace_Channel_Names (originalChannelLabels, OcularChannelsPossibleNames1,  desiredOcularChannelName1 );


OcularChannelsPossibleNames2 = {'RLO','EEGPOLRLO','POLEOGR','POLROC','EEGROCRef','EOGR','ROC','RAE','RUE','RLOC','Leye','POLRUE','REYE','RIO',...
    'EEGREYERef','ROF','REOG','POLRUE','POLRLE','EOGRT'};
desiredOcularChannelName2= 'Eye2';
originalChannelLabels= Replace_Channel_Names (originalChannelLabels, OcularChannelsPossibleNames2,  desiredOcularChannelName2 );


ekgChannelsPossibleNames1 = {'ECGL','ECG1','ecg1','EKGL','LEKG','EEGLEKGRef','ECGLA','EEGECGLRef','EKGLT'};  %you can add 'ECGV2' for BCH too
desiredEKGChannelName1 = 'EKG1';
originalChannelLabels= Replace_Channel_Names (originalChannelLabels, ekgChannelsPossibleNames1,  desiredEKGChannelName1 );
ekgChannelsPossibleNames2 = {'ECGR','ECG2','ecg2','EKGR','REKG','EEGREKGRef','ECGRA','EEGECGRRef','EKGRT'};  %you can add 'ECGV2' for BCH too
desiredEKGChannelName2 = 'EKG2';
originalChannelLabels= Replace_Channel_Names (originalChannelLabels, ekgChannelsPossibleNames2,  desiredEKGChannelName2 );
ekgChannelsPossibleNames3 = {'ECG','EEGEKGRef','POLEKG','EEGPOLEKG'};
desiredEKGChannelName3 = 'EKG';
originalChannelLabels= Replace_Channel_Names (originalChannelLabels, ekgChannelsPossibleNames3,  desiredEKGChannelName3);

emgChannelsPossibleNames = {'CHIN1','CHIN2','NECK1','NECK2','NEC1','NEC2','Lleg1','Lleg2','Rleg1','Rleg2','chin','EEGNeckRef','POLNeck1',...
    'POLNeck2','neck1','neck2','LEMG1','REMG1','EEGCHIN1Ref','EEGCHIN2Ref','POLNECK1','POLNECK2','POLChin1','POLChin2','RLEG','LLEG','EMGR','EMGL',...
    'ABD1','ABD2','EEGABD1Ref','EEGABD2Ref','CHINLT','CHINRT','ABDBLK','ABDWHT','UCHIN','LCHIN'};

desiredEMGChannelName = 'EMG';
originalChannelLabels= Replace_Channel_Names (originalChannelLabels, emgChannelsPossibleNames,  desiredEMGChannelName );

% Additional standard EEG channel replacements
otherEEGChannelsPossibleNames1 = {'T7','EEGT7Ref'}; % %you can add 'CHINz' for BCH too
desiredEEGChannelName1 = 'T3';
originalChannelLabels= Replace_Channel_Names (originalChannelLabels, otherEEGChannelsPossibleNames1,  desiredEEGChannelName1 );


otherEEGChannelsPossibleNames2 = {'T8','EEGT8Ref'}; % %you can add 'CHINz' for BCH too
desiredEEGChannelName2 = 'T4';
originalChannelLabels= Replace_Channel_Names (originalChannelLabels, otherEEGChannelsPossibleNames2,  desiredEEGChannelName2 );

otherEEGChannelsPossibleNames3 = {'P7','EEGP7Ref'}; % %you can add 'CHINz' for BCH too
desiredEEGChannelName3 = 'T5';
originalChannelLabels= Replace_Channel_Names (originalChannelLabels, otherEEGChannelsPossibleNames3,  desiredEEGChannelName3 );

otherEEGChannelsPossibleNames4 = {'P8','EEGP8Ref'}; % %you can add 'CHINz' for BCH too
desiredEEGChannelName4 = 'T6';
originalChannelLabels= Replace_Channel_Names (originalChannelLabels, otherEEGChannelsPossibleNames4,  desiredEEGChannelName4 );



%% %% Renaming and Counting Channels (EKGs with EKG1 and EKG2, and EMGs with EMG1 and EMG2)

%% add a matrix that says how many
% Initialize counters
ekgCounter = 0;
emgCounter = 0;
fzCounter = 0;
o1Counter = 0;
o2Counter = 0;

% Iterate through the channel names to standardize and count them
for i = 1:length(originalChannelLabels)

    if contains(originalChannelLabels{i}, 'EMG')
        emgCounter = emgCounter + 1;
        originalChannelLabels{i} = ['EMG' num2str(emgCounter)];

    elseif  matches(originalChannelLabels{i},"EKG")
        ekgCounter = ekgCounter + 1;
        originalChannelLabels{i} = ['EKG' num2str(ekgCounter)];

    elseif contains(originalChannelLabels{i}, 'Fz')
        fzCounter = fzCounter + 1;

        if fzCounter ==1
            originalChannelLabels{i} = 'Fz';
        else
            originalChannelLabels{i} = ['Fz' num2str(fzCounter)];
        end
    elseif contains(originalChannelLabels{i}, 'O1')
        o1Counter = o1Counter + 1;

        if o1Counter ==1
            originalChannelLabels{i} = 'O1';
        else
            originalChannelLabels{i} = ['O1' num2str(fzCounter)];
        end


    elseif contains(originalChannelLabels{i}, 'O2')
        o2Counter = o2Counter + 1;


        if o2Counter ==1
            originalChannelLabels{i} = 'O2';
        else
            originalChannelLabels{i} = ['O2' num2str(fzCounter)];
        end


    end

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

