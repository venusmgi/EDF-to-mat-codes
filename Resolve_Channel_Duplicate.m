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



% Removing all the channels that had similarity in names but were not
% selected by the user
for i = 1:length(channelNames)
    channelName = channelNames{i};
    idxToInvestigate = find(contains(allChannelLabels,channelName));
    userSelectedIdx = find(contains(allChannelLabels,'correct_'));
    idxChannelToRemove = setdiff(idxToInvestigate,userSelectedIdx); %fiding channels that have our desired pattern but are not the correct channels
    if ~isempty(idxChannelToRemove)
        allChannelLabels{idxChannelToRemove} = 'useless';
    end
end


% Loop through each channel name to check if the current label matches
% any known types and count them

aa = allChannelLabels;
for j = 1:length(allChannelLabels)
    channelLabel = allChannelLabels{j};

    for k = 1:length(channelNames)
        channelName = channelNames{k};

        if contains(lower(channelLabel),lower(channelName))
            channelCounters.(channelName) = channelCounters.(channelName)+1;

            if strcmp(channelLabel,'correct_O1') || strcmp(channelLabel,'correct_O2') ||strcmp(channelLabel,'correct_Fz')
                newLabel = channelName; % Preserve the base name without number for the correct instance
            else

               newLabel = [channelName num2str(channelCounters.(channelName))];
            end
            allChannelLabels{j} = newLabel;
        end

    
    end
end


% Return the unchanged label if no matching channel name is found
outputChannelName = allChannelLabels;



end

% % function [outPutChannelLabels,channelCounters] = Handle_Difficult_Channel (allChannelLabels,channelCounters)

function [outPutChannelLabels] = Handle_Difficult_Channel (allChannelLabels)
difficultChannelNames = {'O1','O2','Fz','EMG','EKG'};
for i = 1:length(difficultChannelNames)
    channelName = difficultChannelNames{i};
    channelOptions = allChannelLabels(contains(lower(allChannelLabels),lower(channelName)));
    uniqueChannelOptions = unique(channelOptions);

    
   if ~isempty(uniqueChannelOptions)
       
        if isscalar(uniqueChannelOptions) %if we have one option, that is the correct one
            % % %adding channel counters
            % % numChannelsCounted = sum(strcmp(allChannelLabels,uniqueChannelOptions));
            % % channelCounters.(channelName) = channelCounters.(channelName) + numChannelsCounted;

            allChannelLabels (strcmp(allChannelLabels,uniqueChannelOptions)) = {['correct_' channelName]};
            

        else % if we have more than one option, ask the user
    
        selection = questdlg(['Which one of these ', channelName, ...
            ' channels is in the main 19 EEG electrodes?'],...
            'Select the Correct Name,(hint: if it is not one of the main EEG electrodes, and you have multiple of the same name, just pick one!)', channelOptions{:},channelOptions{1});
        
        % % %adding channel counters
        % % numChannelsCounted = sum(strcmp(allChannelLabels,selection));
        % % channelCounters.(channelName) = channelCounters.(channelName) + numChannelsCounted;
        allChannelLabels (strcmp(allChannelLabels,selection)) = {['correct_' channelName]};
        end
   end
            

    
end

outPutChannelLabels = allChannelLabels;

end