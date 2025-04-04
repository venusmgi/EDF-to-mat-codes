% header_format
% if you want the output .mat file to have the EEG_lab structure, yo should
% define it as 'EEG_lab', if not, it will have the default format defined
% in "edfreadUntilDone.m"

function Convert_EDFtomat_Desired_Fs_channOrder_Parent_Path_ExcOut_repCh (desired_sampling_rate,desired_channel_order,removable_channels, channel_tobe_replaced, new_channel_names,chann_info_name,header_format)
% Check if the header format is EEGlab-compatible
% Check if the header_format is provided and not empty

functionPath = pwd;
addpath(functionPath)

if nargin < 7
    header_format = '';  % Default value (regular header format)
end

ParentPath = uigetdir;

patient_folders = Find_folders(ParentPath);
% Construct the full path to the diagnosis and follow up folder in each patient folder
for i = 1:length(patient_folders)

    All_DX_edfs = fullfile(ParentPath, patient_folders(i),'\diagnosis\');
    cd (All_DX_edfs{1})
    DX_edfs = dir ('*.edf');
    current_path = pwd;



    for j = 1:length(DX_edfs)

        DX_name = DX_edfs(j).name;
        DX_filePath = fullfile(current_path,DX_name);
        matFileName1 = strrep(DX_name, '.edf', '_reordered.mat');
        matFileName2 = strrep(DX_name, '.edf', '_reordered_resampled.mat');
        matFilePath1 = fullfile(current_path, matFileName1);
        matFilePath2 = fullfile(current_path, matFileName2);


        %%checks to see if there is any mat files present in the folder
        if (exist(matFilePath1)==2)|(exist(matFilePath2)==2)
            display(strcat(matFileName1,' already exists'));
        else
            Convert_EDF2Mat (All_DX_edfs{1},DX_name,desired_sampling_rate,...
                desired_channel_order,removable_channels,channel_tobe_replaced, new_channel_names,...
                header_format,ParentPath,chann_info_name,DX_name);

        end
    end

    All_FU_edfs = fullfile(ParentPath, patient_folders(i),'\follow up\');
    cd (All_FU_edfs{1})
    FU_edfs = dir ('*.edf');
    current_path = pwd;



    for j = 1:length(FU_edfs)

        FU_name = FU_edfs(j).name;
        FU_filePath = fullfile(current_path,FU_name);
        matFileName1 = strrep(FU_name, '.edf', '_reordered.mat');
        matFileName2 = strrep(FU_name, '.edf', '_reordered_resampled.mat');
        matFilePath1 = fullfile(current_path, matFileName1);
        matFilePath2 = fullfile(current_path, matFileName2);


        %%checks to see if there is any mat files present in the folder
        if (exist(matFilePath1)==2) |(exist(matFilePath2)==2)
            display(strcat(matFileName1,' already exists'));
        else
            Convert_EDF2Mat (All_FU_edfs{1},FU_name,desired_sampling_rate,...
                desired_channel_order,removable_channels,channel_tobe_replaced,new_channel_names,...
                header_format,ParentPath,chann_info_name,FU_name);
        end
    end


end
end

%taking all the files that end with .edf




function Convert_EDF2Mat (file_path,file_name,desired_sampling_rate,desired_channel_order,removable_channels,channel_tobe_replaced, new_channel_names,header_format,ParentPath,chann_info_name,edf_name)

addpath(file_path)

isEEGlabFormat = strcmpi(header_format, 'EEGlab');
num_chann = length(desired_channel_order);
[hdr, record] = edfreadUntilDone(string(file_name));
EEG_original_channel_order = hdr.label;
%     True_Fs = hdr.samples(1)/ hdr.duration;
original_Fs = hdr.frequency(1);


%% removing the channels that have similarity in names with our main 19 channels
EEG_record = record;
clear record

[reordered_record,sorted_channels_Indices] =Get_desired_channel_order_excelOutput_replacing_input (desired_channel_order,removable_channels,channel_tobe_replaced, new_channel_names,EEG_record, EEG_original_channel_order,ParentPath,chann_info_name,edf_name);
% Saving the EEG record of the channel order we want
reordered_EEG = make_EEG_header(hdr,desired_channel_order);
reordered_EEG.data = reordered_record;


if original_Fs ~= desired_sampling_rate
    extention = 'reordered_resampled_EEG';
    fname = strcat(file_name(1:end-4),'_reordered_resampled.mat');

else
    extention = 'reordered_EEG';
    fname = strcat(file_name(1:end-4),'_reordered.mat');

end
reordered_EEG.setname = fname;



reordered_record_old = reordered_record;
clear reordered_record

original_Fs = int64(original_Fs);
for i = 1:num_chann
    reordered_record(i,:) = resample(reordered_record_old(i,:) ,desired_sampling_rate,original_Fs);
end

hdr.frequency(1:num_chann) = desired_sampling_rate;


% Saving the EEG record of the channel order we want
if isEEGlabFormat
    reordered_EEG = make_EEG_header_EEGlab(hdr, desired_channel_order);

    save(fname, 'reordered_EEG', '-v7.3');
else
    reordered_hdr = new_header(hdr, sorted_channels_Indices, desired_channel_order);
    save(fname,"reordered_record","reordered_hdr", '-v7.3');
end

display(strcat('done saving ',fname))

clear EEG_original_channel_order EEG_record reordered_record hdr reordered_hdr reordered_EEG fname file_name file_name

end
%%subfunctions
% EEG.chanlocs.labels

function [reordered_hdr]= new_header(header,sortedIndices,desired_channel_order)
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

function [EEG]= make_EEG_header(header,desired_channel_order)
load EEG_struct.mat;
EEG.chanlocs.labels = char(desired_channel_order);
EEG.srate = unique(header.frequency);
EEG.nbchan = length(desired_channel_order);

end



function [folderNames] = Find_folders(ParentPath)
j = 1;
folderNames = {}
All_files_and_folders = dir(ParentPath);
for i = 1: length(All_files_and_folders)

    % skipping non-folders and folders starting with '.' (hidden folders on Unix-like systems)
    foldername = All_files_and_folders(i).name;

    if All_files_and_folders(i).isdir &  ~((strcmpi(foldername, '..')) || (strcmpi(foldername, '.')))

        % Get the folder name
        folderNames{j} = All_files_and_folders(i).name;
        j = j+1;


    end
end
end