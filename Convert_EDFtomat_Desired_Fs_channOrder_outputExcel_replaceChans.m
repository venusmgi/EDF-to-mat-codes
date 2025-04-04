%by Venus 10.23.2024
function Convert_EDFtomat_Desired_Fs_channOrder_outputExcel_replaceChans (desired_sampling_rate,desired_channel_order,removable_channels,channel_tobe_replaced, new_channel_names,chann_info_name, header_format)
    % Check if the header format is EEGlab-compatible
      % Check if the header_format is provided and not empty
    if nargin <7
        header_format = '';  % Default value (regular header format)
    end
    isEEGlabFormat = strcmpi(header_format, 'EEGlab');
    num_chann = length(desired_channel_order);
    
    [filename, path] = uigetfile('*.edf',...
            'Select EDF File(s)', ...
            'MultiSelect', 'on');
    addpath(path)
    
    if class(filename) == 'char'
      total = 1;
    else
        total = length(filename);
    end

    for i = 1:total
         if total == 1
             current_file_name = filename;
            
         else
             current_file_name = filename{i};
         end
        matFileName1 = strrep(current_file_name, '.edf', '_reordered.mat');
        matFileName2 = strrep(current_file_name, '.edf', '_reordered_resampled.mat');
        matFilePath1 = fullfile(path, matFileName1);
        matFilePath2 = fullfile(path, matFileName2);


                %%checks to see if there is any mat files present in the folder
        if (exist(matFilePath1)==2)|(exist(matFilePath2)==2)
            display(strcat(matFileName1,' already exists'));
        else  

             [hdr, record] = edfreadUntilDone(string(current_file_name));
        
        
        
        
    
            EEG_original_channel_order = hdr.label;
        %     True_Fs = hdr.samples(1)/ hdr.duration;
            original_Fs = hdr.frequency(1);
        
         
        
    
            %% removing the channels that have similarity in names with our main 19 channels
            
        
        
            EEG_record = record;
            clear record
            
          [reordered_record,sorted_channels_Indices] = Get_desired_channel_order_excelOutput_replacing_input(desired_channel_order,removable_channels,channel_tobe_replaced, new_channel_names,EEG_record, EEG_original_channel_order,path,chann_info_name,current_file_name);
            % [reordered_record,sorted_channels_Indices] =Get_desired_channel_order_output_excel (desired_channel_order,removable_channels,EEG_record, EEG_original_channel_order,path,chann_info_name,current_file_name);
            % Saving the EEG record of the channel order we want
            reordered_EEG = make_EEG_header(hdr,desired_channel_order);
            reordered_EEG.data = reordered_record;
    
    
            if original_Fs ~= desired_sampling_rate 
                extention = 'reordered_resampled_EEG';
                fname = strcat(current_file_name(1:end-4),'_reordered_resampled.mat');
            
            else
                extention = 'reordered_EEG';
                fname = strcat(current_file_name(1:end-4),'_reordered.mat');
            end
            reordered_EEG.setname = fname;
    
             reordered_record_old = reordered_record;
            clear reordered_record
    
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
        
            clear EEG_original_channel_order EEG_record reordered_record hdr reordered_hdr reordered_EEG fname current_file_name
        end
        end
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