function [outputChannelNames] = Rename_Channel_Labels(originalChannelLabels,channelsToBeReplaced,newChannelNames)
% RENAME_CHANNEL_LABELS Standardizes and replaces EEG channel names.
%
% This function updates EEG channel names by first replacing specified names with new ones (if inputted),
% and then mapping a broad set of potential names to a standardized nomenclature.
%
% Inputs:
%   originalChannelLabels - A cell array containing the original labels of the EEG channels.
%   channelsToBeReplaced - A cell array of channel names to replace in the
%   original labels. (optional, if this is inputted the newChannelNames 
%   should also be an input)
%   newChannelNames - A cell array with new names corresponding to each
%   channel name to be replaced. (optional, if this is inputted the newChannelNames
%   should also be an input)
%
% Outputs:
%   outputChannelNames - A cell array containing updated channel names.
%  
% Note: channelsToBeReplaced and newChannelNames should be the same size 

%% Step 1: Replace specified channels with the provided new names

if nargin > 1
    % Verify that the number of channels to be replaced matches the number of new names
    assert(length(channelsToBeReplaced)==length(newChannelNames),...
        'Number of channels to be replaced (%d) must match number of new names (%d)',...
        length(channelsToBeReplaced), length(newChannelNames));

    %change the channle names
    for i = 1:length(channelsToBeReplaced)
        chanName{1} = channelsToBeReplaced{i}; %added {1} to make sure I am passing a cell
        newChanName= newChannelNames{i}; % this should be a character
        originalChannelLabels = Replace_Channel_Names(originalChannelLabels, chanName, newChanName);
    end
end


%% Step 2: Map all remaining channels to standardized names
% Define potential names and their standardized replacements for various channel types


        channelMappings = {'Eye1',{'LUO','EEGPOLLUO','POLEOGL','POLLLC','POLLOC','EEGLOCRef', 'LOC','EOGL','LLC', 'LUE','LUOC','Reye','POLLLE',...
                                   'LEYE','LIO','EEGLEYERef','LOF','LEOG','POLLLE','POLLUE','LLE','EOGLT'};
        'Eye2',{'RLO','EEGPOLRLO','POLEOGR','POLROC','EEGROCRef','EOGR','ROC','RAE','RUE','RLOC','Leye','POLRUE','REYE','RIO',...
                'EEGREYERef','ROF','REOG','POLRUE','POLRLE','EOGRT'};
        'EKG',{'EKG1','ECGL','ECG1','ecg1','EKGL','LEKG','EEGLEKGRef','ECGLA','EEGECGLRef','EKGLT',... %these are called EKG1 usually
               'EKG2','ECGR','ECG2','ecg2','EKGR','REKG','EEGREKGRef','ECGRA','EEGECGRRef','EKGRT',...%these are called EKG2 usually
               'ECG','EEGEKGRef','POLEKG','EEGPOLEKG'};
        'EMG',{'CHIN1','CHIN2','NECK1','NECK2','NEC1','NEC2','Lleg1','Lleg2','Rleg1','Rleg2','chin','EEGNeckRef','POLNeck1',...
               'POLNeck2','neck1','neck2','LEMG1','REMG1','EEGCHIN1Ref','EEGCHIN2Ref','POLNECK1','POLNECK2','POLChin1',...
               'POLChin2','RLEG','LLEG','EMGR','EMGL','ABD1','ABD2','EEGABD1Ref','EEGABD2Ref','CHINLT','CHINRT','ABDBLK',...
               'ABDWHT','UCHIN','LCHIN'};
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

    % Return the updated channel labels
    outputChannelNames = originalChannelLabels;
end







function [outputChannelNames] = Replace_Channel_Names (originalChannelNames, channelsToBeReplace, replacementChannelName)
% REPLACE_CHANNEL_NAMES Replaces specified channel names with new names.
%
% This function iterates over a list of original channel names and replaces any occurrence of
% specified names with a new name provided by the user.
%
% Inputs:
%   originalChannelNames - A cell array of the channel names to process.
%   channelsToBeReplace - A cell array of channel names that should be replaced.
%   replacement_channel - The new name to assign to each matching channel.
%   a string of characters
%
% Outputs:
%   outputChannelNames - Modified cell array of channel names with replacements applied.

for j = 1:length(channelsToBeReplace)
    channelName = channelsToBeReplace{j};

    % checks if any of the possible other names of the channels that
    % should be replaced exists in the channel names
    if any(strcmp(originalChannelNames, channelName))
        % replaces that channels with the defined  channel name
            %find the indices of the matched channel
            IdxToBeReplaced = strcmp(originalChannelNames, channelName);
            originalChannelNames{IdxToBeReplaced} = replacementChannelName;

    end

end

outputChannelNames = originalChannelNames;

end