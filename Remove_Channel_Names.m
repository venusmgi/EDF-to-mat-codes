function [outPutChannelNames] = Remove_Channel_Names  (originalChannelOrder,removableChannels)

% REMOVE_CHANNEL_NAMES Marks specific EEG channels for removal based on a
% predefined set of channels. Either inputted by the user, or the ones that
% are predefined in this function.
%
% This function checks each channel in the original EEG channel order against a list of channels
% to be removed. It marks any matching channels as 'unusable', allowing the user to effectively
% exclude these channels from further analysis or processing.
%
% Inputs:
%   originalChannelOrder - A cell array of the original EEG channel labels (e.g., {'Fp1', 'Fp2', 'O1', ...}).
%   removableChannels - A cell array of additional channel names that users
%   specifically wish to remove. (optional)
%
% Outputs:
%   outPutChannelNames - A cell array containing the EEG channel names with specified removable channels
%                        marked as 'unusable'.

% Predefined list of channel names that are commonly found to be non-essential or unwanted
% These channel names has been found from the dataset NIMBIS, you might
% need to add your own channels
channelsToBeRemoved = {'SpO2','SpO2Org','etCO2','CO2Wave','CO2WaveOrg',...
    'POLCO2Wave','STERNO1','STERNO2','O11','O21','EtCO2','EtCO2Org',...
    'POLEtCO2','POLSpO2',...
    'AF3', 'AF4','AF7','AF8','AFZ','DIF3','DIF4','DIF7','DIF8',...
    'CP3','DC3','DC4','DC01','DC02','DC03','DC04',...,
    'POLBP3','POLBP4',...
    'EEGFPZRef','EEGFpZRef','EEGFpzRef','FCZ','FPZ','Fpz','Fz2',...
    'PARA1','PARA2','POLA1','POLA2','EEGPOLM1','EEGPOLM2',...
    'LEMG1','LEMG2','LEMG3','LEMG4','POLEMG1','POLEMG2',...
    'REMG1','REMG2','REMG3','REMG4'};

% Append user-provided channels to remove to the predefined list
if ~ isempty(removableChannels)
    channelsToBeRemoved = [channelsToBeRemoved removableChannels];
end


% Mark removable channels as 'unusable'
% Iterate through the list of channels to be removed
for i = 1:length(channelsToBeRemoved)
    % Check if the current channel in the removal list is in the original channel order
    if any(contains(originalChannelOrder,channelsToBeRemoved{i}))
        % Find the index of the removable channel in the original list
        removableChannelIdx = contains(originalChannelOrder,channelsToBeRemoved{i});
        % Mark the channel as 'unusable' in the output list
        originalChannelOrder{removableChannelIdx} = 'unusable';
    end
end
% Set the output variable with the modified channel list
outPutChannelNames = originalChannelOrder;
end