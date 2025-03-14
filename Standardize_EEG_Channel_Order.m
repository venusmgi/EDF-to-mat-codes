
% 10.23.2024 by Venus
function [reorderd_record,sorted_channels_Indices] = Standardize_EEG_Channel_Order (desired_channel_order,removable_channels,channel_tobe_replaced, new_channel_names, ...
    EEG_record,EEG_original_channel_order,ParentPath,chann_info_name,edf_name)

% Previoulsy was called Get_desired_channel_order_excelOutput_replacing_input

% This function reorders EEG channels based on a desired order, replaces specified channel names,
% and generates a report of the changes. It handles standard EEG, EKG, EMG, and ocular channels.
%
% Inputs:
%   desired_channel_order - Desired order of EEG channels
%   removable_channels - Channels to be excluded from the final order
%   channel_tobe_replaced - Channels with incorrect names to be replaced(
%   e.g. a channel that is named X1, but we know that it is EKG)
%   new_channel_names - Correct names for the channels to be replaced (EKG
%   for X1)
%   EEG_record - EEG data matrix (channels x time points)
%   EEG_original_channel_order - Original order of EEG channel names
%   ParentPath - Path to save the channel information
%   chann_info_name - Name for the channel information file
%   edf_name - Name of the EDF file being processed
%
% Outputs:
%   reorderd_record - EEG data reordered according to the desired channel order
%   sorted_channels_Indices - Indices of the sorted channels in the original data

% Note : the order and length of channels in channel_tobe_replaced and 
% new_channel_names should be the same




assert (length(channel_tobe_replaced)==length(new_channel_names),'number of channels that are being replaced and their new names should be the same')
%% standardizing the Channel names
  Orig_chn_label = EEG_original_channel_order;  
    if ~isempty(channel_tobe_replaced) && ~isempty(new_channel_names)
        %replacing channels that were named not-correctly
        for i = 1:length(channel_tobe_replaced)
            chan_name{1} = channel_tobe_replaced{i};
            new_chan_name = new_channel_names{i};
            Orig_chn_label = replace_channel_names (Orig_chn_label, chan_name,  new_chan_name );
    
        end
    end

    
    
    Possible_ocular_channel_names1 = {'LUO','EEGPOLLUO','POLEOGL','POLLLC','POLLOC','EEGLOCRef', 'LOC','EOGL','LLC', 'LUE','LUOC','Reye','POLLLE',...
        'LEYE','LIO','EEGLEYERef','LOF','LEOG','POLLLE','POLLUE','LLE','EOGLT'}; %--> you did not keep thi
    Ocular_channel_name1 = 'Eye1';
    Orig_chn_label = replace_channel_names (Orig_chn_label, Possible_ocular_channel_names1,  Ocular_channel_name1 );
   

    Possible_ocluar_channel_names2 = {'RLO','EEGPOLRLO','POLEOGR','POLROC','EEGROCRef','EOGR','ROC','RAE','RUE','RLOC','Leye','POLRUE','REYE','RIO',...
        'EEGREYERef','ROF','REOG','POLRUE','POLRLE','EOGRT'};
    Ocualr_channel_name2= 'Eye2';
    Orig_chn_label = replace_channel_names (Orig_chn_label, Possible_ocluar_channel_names2,  Ocualr_channel_name2 );
   
    
    Possible_EKG_channel_names1 = {'ECGL','ECG1','ecg1','EKGL','LEKG','EEGLEKGRef','ECGLA','EEGECGLRef','EKGLT'};  %you can add 'ECGV2' for BCH too
    EKG_channel_name1 = 'EKG1';
    Orig_chn_label = replace_channel_names (Orig_chn_label, Possible_EKG_channel_names1,  EKG_channel_name1 );
     Possible_EKG_channel_names2 = {'ECGR','ECG2','ecg2','EKGR','REKG','EEGREKGRef','ECGRA','EEGECGRRef','EKGRT'};  %you can add 'ECGV2' for BCH too
    EKG_channel_name2 = 'EKG2';
    Orig_chn_label = replace_channel_names (Orig_chn_label, Possible_EKG_channel_names2,  EKG_channel_name2 );
    Possible_EKG_channel_names3 = {'ECG','EEGEKGRef','POLEKG','EEGPOLEKG'};
    EKG_channel_name3 = 'EKG';
    Orig_chn_label = replace_channel_names (Orig_chn_label, Possible_EKG_channel_names3,  EKG_channel_name3);


    

   % %you can add 'CHINz' for BCH too
   % % also you can add 'CHINz',,'LEMG2','REMG2','LEMG3','REMG3','LEMG4','REMG4'for CCRO
    Possible_EMG_channel_names = {'CHIN1','CHIN2','NECK1','NECK2','NEC1','NEC2','Lleg1','Lleg2','Rleg1','Rleg2','chin','EEGNeckRef','POLNeck1',...
        'POLNeck2','neck1','neck2','LEMG1','REMG1','EEGCHIN1Ref','EEGCHIN2Ref','POLNECK1','POLNECK2','POLChin1','POLChin2','RLEG','LLEG','EMGR','EMGL',...
        'ABD1','ABD2','EEGABD1Ref','EEGABD2Ref','CHINLT','CHINRT','ABDBLK','ABDWHT','UCHIN','LCHIN'};
    
    
    EMG_channel_name = 'EMG';
    Orig_chn_label = replace_channel_names (Orig_chn_label, Possible_EMG_channel_names,  EMG_channel_name );

     Possible_EEG_channel_name1 = {'T7','EEGT7Ref'}; % %you can add 'CHINz' for BCH too
    EEG_channel_name1 = 'T3';
    Orig_chn_label = replace_channel_names (Orig_chn_label, Possible_EEG_channel_name1,  EEG_channel_name1 );


    Possible_EEG_channel_name2 = {'T8','EEGT8Ref'}; % %you can add 'CHINz' for BCH too
    EEG_channel_name2 = 'T4';
    Orig_chn_label = replace_channel_names (Orig_chn_label, Possible_EEG_channel_name2,  EEG_channel_name2 );

    Possible_EEG_channel_name3 = {'P7','EEGP7Ref'}; % %you can add 'CHINz' for BCH too
    EEG_channel_name3 = 'T5';
    Orig_chn_label = replace_channel_names (Orig_chn_label, Possible_EEG_channel_name3,  EEG_channel_name3 );

    Possible_EEG_channel_name4 = {'P8','EEGP8Ref'}; % %you can add 'CHINz' for BCH too
    EEG_channel_name4 = 'T6';
    Orig_chn_label = replace_channel_names (Orig_chn_label, Possible_EEG_channel_name4,  EEG_channel_name4 );



    %% renaming EKGs with EKG1 and EKG2, and EMGs with EMG1 and EMG2

    %% add a matrix that says how many 
    % Initialize counters
    ekgCounter = 0;
    emgCounter = 0;
    fzcounter = 0;
    o1counter = 0;
    o2counter = 0;

    % Iterate through the channel names
    for i = 1:length(Orig_chn_label)
        % Check if the channel is EKG


        if contains(Orig_chn_label{i}, 'EMG')
            emgCounter = emgCounter + 1;
            Orig_chn_label{i} = ['EMG' num2str(emgCounter)];

        elseif  matches(Orig_chn_label{i},"EKG")
            ekgCounter = ekgCounter + 1;
            Orig_chn_label{i} = ['EKG' num2str(ekgCounter)];

        elseif contains(Orig_chn_label{i}, 'Fz')
             fzcounter = fzcounter + 1;

             if fzcounter ==1
                Orig_chn_label{i} = 'Fz';
             else
                 Orig_chn_label{i} = ['Fz' num2str(fzcounter)];
             end
       elseif contains(Orig_chn_label{i}, 'O1')
           o1counter = o1counter + 1;

              if o1counter ==1
                Orig_chn_label{i} = 'O1';
             else
                 Orig_chn_label{i} = ['O1' num2str(fzcounter)];
              end
              

       elseif contains(Orig_chn_label{i}, 'O2')
           o2counter = o2counter + 1;


             if o2counter ==1
                Orig_chn_label{i} = 'O2';
             else
                 Orig_chn_label{i} = ['O2' num2str(fzcounter)];
             end
             

        end

    end


    
    %% first finding the indecies pof the channels we want to keep and saving the new channel label info
   
    [including_channels_lables,including_channels_indecies]=find_indecies_and_channs_of_desired(Orig_chn_label,desired_channel_order,removable_channels);
            % [~,Indecies_of_new_channel_names] = sort(Untoched_sortedIndices)

            
            
  


 excel_sheet_path = strcat(ParentPath,'\',chann_info_name,'.mat');
% Create a cell array to store the information
chanInfo = cell(3,1 + length(Orig_chn_label)); % 3 columns: Original Name, New Name, Skipped
row1Name = edf_name;
row2Name = strcat('New',' ',edf_name);
row3Name = strcat('Reordered',' ',edf_name);
% Fill in the header
chanInfo(:,1) = {row1Name, row2Name,row3Name};

% Fill in the data
for i = 1:length(EEG_original_channel_order)
    chanInfo{ 1,i + 1} = EEG_original_channel_order{i};
    
    % Check if the channel was renamed
    if any(including_channels_indecies == i)
        chanInfo{2, i + 1} = including_channels_lables{including_channels_indecies == i};
        chanInfo{3, i+1}= desired_channel_order{including_channels_indecies == i};
         
        
    else
        chanInfo{2, i + 1} = '';
        chanInfo{3, i+1}= '';
    end
end

excel_sheet_path = strcat(ParentPath,'\',chann_info_name,'.mat');
% Check if the file already exists
if exist(excel_sheet_path) == 2
    % Load existing data
    existingData = load(excel_sheet_path);

   
    previousEDFnames = existingData.chanInfo(:,1);
    currentEDFname = chanInfo(1,1);

    % finding if the previous data had aready done the current edf, to
    % remove that edf and anything after that
       EDFtoRemove =  find(strcmp(previousEDFnames, currentEDFname), 1, 'first');
       if ~isempty(EDFtoRemove)
         existingData.chanInfo(EDFtoRemove:end,:) =[];
       end
     
 
    numColumns1 = size(chanInfo,2);
    numColumns2 = size(existingData.chanInfo,2);
    %checks to see if the size of the new channel infor is the same
    %as the uploaded on, if not, padds either to be the same size
    if numColumns2 > numColumns1
        % Pad cellArray2 with empty cells
        padding = cell(size(chanInfo, 1), numColumns2 - numColumns1);
        chanInfo = [chanInfo, padding];
    elseif numColumns2 < numColumns1
        padding = cell(size(existingData.chanInfo, 1), numColumns1 - numColumns2);
        existingData.chanInfo = [existingData.chanInfo, padding];
    end

    % Append the new data
    chanInfo = [existingData.chanInfo; chanInfo];
    

    % Write the updated data to the Excel file
    save(excel_sheet_path,'chanInfo')
else
    % Write the data to a new Excel file
    save(excel_sheet_path,'chanInfo');
end





    %% finding the name of the channels in the header that correspond to our 19 desired channels and their indecies
     % removing the channels that have similarity in names with our main 19 channels and % finding the indecies of the channel order we want
     % replacing the channels that we are keeping with the new name,
            % and if we are not using it, replacaing it with empty string
    
    
        
    record_desired_channs = EEG_record(including_channels_indecies,:);
    
    [~, sorted_channels_Indices] = ismember (desired_channel_order,including_channels_lables);
    if any(sorted_channels_Indices ==0)
       sorted_channels_Indices =  find_sorted_indecies(including_channels_lables,desired_channel_order);
    end



    % Saving the EEG record of the channel order we want
    reorderd_record = record_desired_channs(sorted_channels_Indices,:);
end



%% subfunctions


function [including_channels_order,including_chann_indecies] = find_indecies_and_channs_of_desired(chn_label,desired_order,removable_channels)

% %this part replaces A1 and A2, with M1 and M2
%     if (any(contains(chn_label,'A1').* ~strcmp(chn_label,'POL $A1')) & any(strcmp(desired_order,'M1')))
%         chn_label{  find(contains(chn_label,'A1').* ~strcmp(chn_label,'POL $A1'))  } ='M1';
%     elseif (any(contains(chn_label,'M1').* ~strcmp(chn_label,'POLM1')) & any(strcmp(desired_order,'A1')))
%         chn_label{  find(contains(chn_label,'M1').* ~strcmp(chn_label,'POLM1'))  } ='A1';
%     end
% 
%     if (any(contains(chn_label,'A2').* ~strcmp(chn_label,'POL$A2')) & any(strcmp(desired_order,'M2')))
%         chn_label{  find(contains(chn_label,'A2').* ~strcmp(chn_label,'POL$A2'))  } ='M2';
%     elseif (any(contains(chn_label,'M2').* ~strcmp(chn_label,'POLM2')) & any(strcmp(desired_order,'A2')))
%         chn_label{  find(contains(chn_label,'M2').* ~strcmp(chn_label,'POLM2'))  } ='A2';
%     end
        
    for i = 1:length(chn_label)
            contains_required_channels = contains(lower(chn_label(1,i)) , lower(desired_order) );
            does_not_contain_removable_channels = ~( contains (chn_label(1,i) ,removable_channels));
            including_chann_indecies(:,i) =  contains_required_channels && does_not_contain_removable_channels;
    end

    % for i = 1:length(chn_label)
    %         contains_required_channels1 = contains(lower(chn_label{i}) , lower(desired_order) );
    %         does_not_contain_removable_channels1 = ~( contains (chn_label{i} ,removable_channels));
    %         including_chann_indecies1(:,i) =  contains_required_channels && does_not_contain_removable_channels;
    % end

            
    if sum(including_chann_indecies) ~= length(desired_order)
            error(strcat('Number of channels is smaller or higher than number of desired channel. check channel labels.',...
              'check variable "including_chann_indecies" here and see what channels are missing from "chan_label" or are being added aditionally.',...
             'if "including_chann_indecies" is 1, it means the corresponding channel in "chan_label"is being included, if you do not want it put it in input "removable_channels".',...
             'if in "including_chann_indecies" the corresponding index to "chan_label" is 0, and you want that channel, add it to input "desired_channel_order".',... 
             'if it is 1 and you do not want it add the channel name in "chan_label" to "removable_channels"' ))
            return
        end
 % contains(lower(chn_label(1,i)) , lower(desired_order) );
    including_chann_indecies=find(including_chann_indecies);

    including_channels_order = chn_label(including_chann_indecies);

end

function [sortedIndices] = find_sorted_indecies (EEG_remaining_channel_labels,desired_order)
     sortedIndices = []; %sometimes the ismemeber does not work, so we go through this loop to find the indecies based on our desired order
     k = 1;
        for i = 1: numel(desired_order)
            lower_desired_order = lower(desired_order{i});
            for j = 1:numel(EEG_remaining_channel_labels)
                original_chan_label = lower(EEG_remaining_channel_labels{j});
                if contains(original_chan_label,lower_desired_order)
                    sortedIndices(1,k) = j;
                    k = k+1; %%for checking contanins put a pause excecution here
                end
            end
        end


end



function [output_channelNames] = replace_channel_names (original_channel_names, channels_to_be_replace, replacement_channel)

    for j = 1:length(channels_to_be_replace)

        %checks is any of the possible other names of the channels that
        %should be replaced
        %exists in the channel names
        if any(strcmp(original_channel_names, channels_to_be_replace{j}))
            % replaces that channels with the defined  channel name
            if sum(strcmp(original_channel_names, channels_to_be_replace{j})) ==1
                original_channel_names{strcmp(original_channel_names, channels_to_be_replace{j})} = replacement_channel;
            else
                
               Idx_tobe_replaced =  find(strcmp(original_channel_names, channels_to_be_replace{j}));
               for i = 1:length(Idx_tobe_replaced)
                   original_channel_names{Idx_tobe_replaced(i)} = replacement_channel;
               end
            end
        end
        
    end

    output_channelNames = original_channel_names;

end

