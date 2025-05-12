function [outputChannelName] = Resolve_Channel_Duplicate(allChannelLabels)

% STANDARDIZECHANNELNAME Renames EEG channels to ensure standardized naming conventions.
%
% This function adjusts channel names to add numerical suffixes to duplicate channels,
% ensuring that each channel has a unique name. For main EEG channels like 'O1', 'O2', and 'Fz',
% the base name is preserved for the first occurrence without adding a numerical suffix.
%
% Inputs:
%   allChannelLabels - A cell containing the original label of the EEG channel.
%                      (e.g. ('Fp1','Fp2','O1,'O2...})
%   channelCounters - A struct containing counters for each channel type
%                     (e.g., 'EMG',0, 'EKG',0, 'Fz', 1, 'O1', 1, 'O2', 0),
%                       which tracks the number of occurrences of each channel type.
%
% Outputs:
%   newLabel - The standardized channel label. If the channel is a duplicate,
%              it includes a number suffix (e.g., 'EMG1', 'EMG2').




% Define the channel names of interest
channelNames = {'EMG', 'EKG', 'Fz', 'O1', 'O2'};


% Initialize a struct to hold channel name counters
channelCounters = struct('EMG', 0, 'EKG', 0, 'Fz', 0, 'O1', 0, 'O2', 0);

%Picking the correct EEG channel and adding a 'correct_' to the name
allChannelLabels = Handle_Difficult_Channel (allChannelLabels);


%making sure there will not be multipe EMG1, EMG2, EKG1, and EKG2
%becuase if first we encounter a channel with EMG or EKG in the name, but it is
%not out desired channel, that one willbe named EMG1, or EKG1

%instead of this I could have started EMG and EKG numbers from 2 instead of
% 0 

channs = allChannelLabels(contains(allChannelLabels,'correct_'));
if any(contains(channs,'EMG'))
    channelCounters.('EMG') = sum(contains(allChannelLabels,'EMG'));
elseif any(contains(channs,'EKG'))
    channelCounters.('EKG') = sum(contains(allChannelLabels,'EKG'));
end



% Iterate through the channel names to standardize and count them
for j = 1:length(allChannelLabels)
    channelLabel = allChannelLabels{j};
    % Loop through each channel name to check if the current label matches any known types
    for i = 1:length(channelNames)
        channelName = channelNames{i};
    
        if contains(lower(channelLabel),lower(channelName))
            % Update the counter for the channel
            channelCounters.(channelName) = channelCounters.(channelName) + 1;
            % For the main EEG channels ('O1', 'O2', 'Fz'), preserve the
            % base name without adding a number, as they are standard electrode names
            if strcmp(channelLabel,'correct_O1') || strcmp(channelLabel,'correct_O2') ||strcmp(channelLabel,'correct_Fz')
                 newLabel = channelName; % Preserve the base name without number for the correct instance
            elseif strcmp(channelLabel,'correct_EMG1') 
                newLabel = 'EMG1';
                elseif strcmp(channelLabel,'correct_EMG2') 
                newLabel = 'EMG2';
                elseif strcmp(channelLabel,'correct_EKG1') 
                newLabel = 'EKG1';
                elseif strcmp(channelLabel,'correct_EKG2') 
                newLabel = 'EKG2';
            else
                %% We don't need extra channels then keep the correct ones
                % Append a numeral suffix for duplicate or non-main channels
                newLabel = [channelName num2str(channelCounters.(channelName))];
            end

            allChannelLabels{j} = newLabel;
            
        end
    end

    
end
% Return the unchanged label if no matching channel name is found
outputChannelName = allChannelLabels;



end

function [outPutChannelLabels] = Handle_Difficult_Channel (allChannelLabels)

difficultChannelNames = {'O1','O2','Fz','EMG1','EMG2','EKG1', 'EKG2'};
for i = 1:length(difficultChannelNames)
    channelName = difficultChannelNames{i};
    channelOptions = allChannelLabels(contains(lower(allChannelLabels),lower(channelName)));
   if ~isempty(channelOptions)
        if isscalar(channelOptions) %if we have one option, that is the correct one
    
            allChannelLabels (strcmp(allChannelLabels,channelOptions)) = {['correct_' channelName]};
        else % if we have more than one option, ask the user
    
        selection = questdlg(['Which one of these ', channelName, ...
            ' channels is in the main 19 EEG electrodes?'],...
            'Select Main EEG Electrode', channelOptions{:},channelOptions{1});
            allChannelLabels (strcmp(allChannelLabels,selection)) = {['correct_' channelName]};
        end
   end
            

    
end

outPutChannelLabels = allChannelLabels;

end