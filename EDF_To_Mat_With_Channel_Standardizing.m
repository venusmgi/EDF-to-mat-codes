function EDF_To_Mat_With_Channel_Standardizing (filePath,fileName,desiredSamplingRate,desiredChannelOrder,removableChannels,channelsToBeReplaced, ...
    newChannelNames,headerFormat,PathToSaveChannelInfo,channInfoName)
% previously called CONVERT_EDF2MAT
% CONVERT_EDF2MAT Converts an EDF file to MAT format with specified parameters
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


% Add file path to MATLAB path
addpath(filePath)

% Determine if EEGlab format is requested
isEEGlabFormat = strcmpi(headerFormat, 'EEGlab');

% Get number of channels and read EDF file
numChannels = length(desiredChannelOrder);
%convert edf to .mat
[hdr, eegRecord] = edfreadUntilDone(string(fileName));
%saving header info
EEGOriginalChannelOrder = hdr.label;
% finding the original sampling rate as it is in the header
originalFs = int64(hdr.frequency(1));  %TODO: figure out a better way to find the original fs


% Process channel order and get reordered data
[reordered_record,sortedChannelsIndices] =Standardize_EEG_Channel_Order (desiredChannelOrder,removableChannels,channelsToBeReplaced, ...
    newChannelNames,eegRecord, EEGOriginalChannelOrder,PathToSaveChannelInfo,channInfoName,fileName);

% Determine output filename based on resampling
fname = Get_Output_Filename(fileName, originalFs, desiredSamplingRate);

% Resample if needed
if originalFs ~= desiredSamplingRate
    reordered_record = Resample_Data(reordered_record,  desiredSamplingRate, originalFs, numChannels);
    hdr.frequency(1:numChannels) = desiredSamplingRate;
end



% Save the processed data
fprintf('Saving processed data...\n');

if isEEGlabFormat
    %in case we wanted it to be in EEGLab format
    % Saving the EEG record of the channel order we want
    reordered_EEG = Make_EEG_Header(hdr,desiredChannelOrder);
    reordered_EEG.setname = fname;
    reordered_EEG.data = reordered_record;
    save(fname, 'reordered_EEG', '-v7.3');
else
    reordered_hdr = New_Header(hdr, sortedChannelsIndices, desiredChannelOrder);
    save(fname,"reordered_record","reordered_hdr", '-v7.3');
end


fprintf('Successfully saved: %s\n', fname);


end




%% subfunctions

function [reordered_hdr]= New_Header(header,sortedIndices,desired_channel_order)
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

function [EEG]= Make_EEG_Header(header,desiredChannelOrder)
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

