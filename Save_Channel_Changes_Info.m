

function Save_Channel_Changes_Info(originalChannelOrder,desiredChannelOrder,edfFileName,includedChannelIndices,includedChannelLabels,pathToSaveReport,channelInfoFileName)
% Create a cell array to store channel information
numChannels = length(originalChannelOrder);
channelInfo = cell(3,1 + numChannels); % 3 rows: Original Name, New Name, Reordered Name

% Set row headers
channelInfo(:,1) = {edfFileName;
    strcat('New',' ',edfFileName);
    strcat('Reordered',' ',edfFileName)
    };



% Fill in channel info
for i = 1:length(originalChannelOrder)
    channelInfo{ 1,i + 1} = originalChannelOrder{i};

    % Check if the channel was renamed
    if any(includedChannelIndices == i)
        channelInfo{2, i + 1} = includedChannelLabels{includedChannelIndices == i};
        channelInfo{3, i+1}= desiredChannelOrder{includedChannelIndices == i};


    else % skip this channel
        channelInfo{2, i + 1} = '';
        channelInfo{3, i+1}= '';
    end
end

% Define path to save the report as a .mat file
reportPath = fullfile(pathToSaveReport, strcat(channelInfoFileName, '.mat'));

% Load existing channel information if the file exists
if isfile(reportPath)


    % Load existing data
    existingData = load(reportPath);
    existingChannelInfo = existingData.channelInfo;

    previousEDFnames = existingChannelInfo(:,1);


    currentEDFname = channelInfo(1,1);

    % finding if this current edf has already been processed and the info
    % has been saved inthe ecel sheet. If this info is already availabe
    % in the report mat file, remove that edf and an other edf after that,
    % because the rest of the files will be processed again

    % TODO : right now, if there is a file, with the same name as it is in
    % the input, it will load that file, and if it does not find a matching
    % name, it will keep everything, but maybe, some EDFs have been skipped
    % before, and now that we are processing them, they will be done again
    % (imagin that the report mat file has only info on the 3rd edf in a file,
    % and we start the process from begining, at fisrt, because the fisrt 
    % edf is beign processed, it cannot find the name of the 3rd one, and
    % does the process, until it reaches the third one, then it will remove
    % all the process, because it found the name of the third, and removed
    % everything after, which might be the first or second one.
    EDFtoRemove =  find(strcmp(previousEDFnames, currentEDFname), 1, 'first');
    if ~isempty(EDFtoRemove)
        existingChannelInfo(EDFtoRemove:end,:) =[];
    end
    channelInfo = Append_New_ChannInfo (channelInfo, existingChannelInfo);

end
   

    % Write the data to a new report mat file
    save(reportPath,'channelInfo');
end


function newChannelInfo = Append_New_ChannInfo (newChannelInfo, existingChannelInfo)

numColumnsChannInfo = size(newChannelInfo,2);
numColumnsExistingChans = size(existingChannelInfo,2);
%checks to see if the size of the new channel info is the same
%as the uploaded on, if not, padds to be the same size

if numColumnsExistingChans > numColumnsChannInfo
    % Pad channelInfo with empty cells to match existingChannelInfo
    newChannelInfo = Pad_Cell_Array(newChannelInfo, size(newChannelInfo, 1), numColumnsExistingChans);
elseif numColumnsExistingChans < numColumnsChannInfo
    % Pad existingChannelInfo with empty cells to match channelInfo
    existingChannelInfo = Pad_Cell_Array(existingChannelInfo, size(existingChannelInfo, 1), numColumnsChannInfo);
end


% Append the new data
newChannelInfo = [existingChannelInfo; newChannelInfo];
end



function paddedArray = Pad_Cell_Array(cellArray, numRows, numColumns)
% PADCELLARRAY Pads a cell array with empty cells to match desired dimensions
%
% Inputs:
%   cellArray - The cell array to pad
%   numRows - Number of rows in the padded array
%   numColumns - Number of columns in the padded array
%
% Outputs:
%   paddedArray - The padded cell array

% Create padding of empty cells
padding = cell(numRows, numColumns - size(cellArray, 2));

% Combine original array with padding
paddedArray = [cellArray, padding];
end
% Ensure the new and existing data have the same number of columns