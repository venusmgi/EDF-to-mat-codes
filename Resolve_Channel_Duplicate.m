function [outputChannelName] = Resolve_Channel_Duplicate(allChannelLabels)

% RESOLVE_CHANNEL_DUPLICATE Standardizes EEG channel names to ensure uniqueness.
%
% This function adjusts channel names by appending numerical suffixes to resolve duplicates,
% ensuring each channel has a unique name while preserving standard names for key EEG channels
% like 'O1', 'O2', and 'Fz'. It allows user input to resolve ambiguities among similarly named channels.
%
% Inputs:
%   allChannelLabels - A cell array containing the original labels of the EEG channels
%                      (e.g., {'Fp1', 'Fp2', 'O1', 'O2', ...}).
%
% Outputs:
%   outputChannelName - A cell array containing standardized channel names. For duplicates,
%                       names include a numeric suffix (e.g., 'EMG1', 'EMG2').




% Define the key EEG channels of interest
channelNames = {'EMG', 'EKG', 'Fz', 'O1', 'O2'};

% Initialize a struct to count occurrences of each channel type
channelCounters = struct('EMG', 0, 'EKG', 0, 'Fz', 0, 'O1', 0, 'O2', 0);

% Handle user input for difficult channel selection and mark the correct choice
allChannelLabels = Handle_Difficult_Channel(allChannelLabels);

% Remove labels that are similar but not selected as correct by the user
for i = 1:length(channelNames)
    channelName = channelNames{i};
    idxToInvestigate = find(contains(allChannelLabels,channelName));  % Indices of channels matching the name
    userSelectedIdx = find(contains(allChannelLabels,'correct_')); % Indices of user-corrected channels
    idxChannelToRemove = setdiff(idxToInvestigate,userSelectedIdx); %finding channels that have our desired pattern but are not the correct channels
    if ~isempty(idxChannelToRemove)
        allChannelLabels{idxChannelToRemove} = 'useless';
    end
end


% Loop through each channel name to check if the current label matches
% any known types and count them

for j = 1:length(allChannelLabels)
    channelLabel = allChannelLabels{j};

    for k = 1:length(channelNames)
        channelName = channelNames{k};

        if contains(lower(channelLabel),lower(channelName))
            % Increment the counter for the channel type
            channelCounters.(channelName) = channelCounters.(channelName)+1;

            if strcmp(channelLabel,'correct_O1') || strcmp(channelLabel,'correct_O2') ||strcmp(channelLabel,'correct_Fz')
                newLabel = channelName; % Preserve the base name without number for O1, O2, and Fz channels
            else
                % Append numeric suffix for EMGs andEKGs
               newLabel = [channelName num2str(channelCounters.(channelName))];
            end

             % Update the standardized channel name in the array
            allChannelLabels{j} = newLabel;
        end

    
    end
end


% Output the processed channel labels and return the unchanged label if no matching channel name is found
outputChannelName = allChannelLabels;



end

function [outPutChannelLabels] = Handle_Difficult_Channel (allChannelLabels)
% HANDLE_DIFFICULT_CHANNEL Asks the user to resolve ambiguities in key channel naming.
%
% This helper function identifies channels that are difficult to standardize and
% prompts the user to select the correct name among potentially ambiguous options.
% Selected names are marked with 'correct_' for easier identification later.
%
% Inputs:
%   allChannelLabels - A cell array of channel labels where potential ambiguities exist.
%
% Outputs:
%   outPutChannelLabels - Updated channel labels with user-confirmed names marked as 'correct_'.

% Define the key channels to address user selection

difficultChannelNames = {'O1','O2','Fz','EMG','EKG'};
for i = 1:length(difficultChannelNames)
    channelName = difficultChannelNames{i};
    channelOptions = allChannelLabels(contains(lower(allChannelLabels),lower(channelName)));
    uniqueChannelOptions = unique(channelOptions); % Identify unique options

   if ~isempty(uniqueChannelOptions)
        if isscalar(uniqueChannelOptions) %if we have one option, that is the correct one
            allChannelLabels (strcmp(allChannelLabels,uniqueChannelOptions)) = {['correct_' channelName]};
            
        else % if we have more than one option, ask the user
        selection = questdlg(['Which one of these ', channelName, ...
            ' channels is in the main 19 EEG electrodes?'],...
            'Select the Correct Name,(hint: if it is not one of the main EEG electrodes, and you have multiple of the same name, just pick one!)', channelOptions{:},channelOptions{1});
      
        % Mark the user's choice with 'correct_'
        allChannelLabels (strcmp(allChannelLabels,selection)) = {['correct_' channelName]};
        end
   end
            

    
end
% Output the processed channel labels and return the unchanged label if no
% change was applied
outPutChannelLabels = allChannelLabels;

end