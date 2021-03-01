%
%
% Project: WORD SEEKING 
% 
% Script to perform EEG preprocessing with EEGlab and ERPlab
% 

%% Open EEGlab

eeglab; 

%% Directories

mainfolder = cd; % Current directory 
rawdata_dir = [cd '/RawData']; % Pre-existing folder with raw EEG data

%% Load Experiment Parameters (.txt file) 


%% Preprocessing

loadfailure = []; % For keeping track of failure to load the dataset


for actualsbj = sbjs % Load datasets

% Determine the name of the file 
if str2double(num2str(actualsbj)) < 10 
sbjcode = sprintf('0%s', num2str(actualsbj)); 
elseif str2double(num2str(actualsbj)) >= 10 
sbjcode = sprintf('%s', num2str(actualsbj));  
end 

actualsbj_id = sprintf('%s%s.%s', FileName{1,1}, sbjcode, FileFormat);
actualsbj_dir = fullfile(rawdata_dir, actualsbj_id);

%-------------------------------------------------------------------------
% Load the Raw Data (determine the format in the Parameters script)
%-------------------------------------------------------------------------
fprintf('Processing file number %d corresponding to participant %s\n\n\n', actualsbj, actualsbj_id);

try 
switch FileFormat
case 'vhdr'
EEG = pop_loadbv(RawDataDir, actualsbj_id); % Load the Raw file // RawFilesList(:,ActualParticipant)
case 'cnt'
EEG = pop_loadcnt(actualsbj_dir, 'dataformat', 'auto', 'memmapfile', ''); % Load the Raw file // RawFilesList(:,ActualParticipant)
end
catch CorruptFile
sprintf('err', 'Failed to load file %d corresponding to participant %s\n\n\n', actualsbj, actualsbj_id)
CorruptFiles = CorruptFiles + 1;
loadfailure(1,CorruptFiles) = str2double(actualsbj_id);
continue
end  


%-------------------------------------------------------------------------
% Re-reference
%-------------------------------------------------------------------------

if numel(RereferenceChannels) > 0
EEG = pop_reref(EEG, RereferenceChannels, 'keepref','on');
end

%-------------------------------------------------------------------------
% Notch filter
%-------------------------------------------------------------------------

if numel (NotchFilter) > 0
        EEG  = pop_basicfilter( EEG, 1:EEG.nbchan,...
                                'Boundary', BoundariesLabel,...
                                'Cutoff', NotchFilter,...
                                'Design', 'notch',...
                                'Filter', 'PMnotch', 'Order',  180);
end
    

if BipolarChannels == 1
    %-------------------------------------------------------------------------
    % Create Bipolar EOG Channels
    %-------------------------------------------------------------------------
    % Vertical (monitors blinking)
    EEG = pop_eegchanoperator( EEG, {'ch24 = ch20 - ch21 label VEOG'} , 'ErrorMsg', 'popup' );
    % Horizontal (monitors saccadic eye movements)
    EEG = pop_eegchanoperator( EEG, {'ch25 = ch22 - ch23 label HEOG'} , 'ErrorMsg', 'popup' );

    BipolarEOGChannels = zeros(1,2);
    for ActualChannel = 1:numel(EEG.chanlocs)
    if strcmp('VEOG', EEG.chanlocs(1,ActualChannel).labels) == 1
    BipolarEOGChannels(1,1) = ActualChannel;
    elseif strcmp('HEOG', EEG.chanlocs(1,ActualChannel).labels) == 1
    BipolarEOGChannels(1,2) = ActualChannel;
    end
    end
end

%-------------------------------------------------------------------------
% Create the EventList of each participant and save it (.txt)
%-------------------------------------------------------------------------
EEG  = pop_creabasiceventlist( EEG , 'AlphanumericCleaning', 'on',...
                               'BoundaryNumeric', { -99 },...
                               'BoundaryString', { 'boundary' },...
                               'Eventlist', sprintf('%sP%s_%s_EventList.txt', EventListFilesDir, actualsbj_id, ExperimentName)); 

%-------------------------------------------------------------------------
% Create the BinLister for each participant and each ERP 
%-------------------------------------------------------------------------
% The BinLister.txt files must be created beforehand and stored in the BinLister directory. 
% There must be one BinLister for each ERP component. 

EEG_raw = EEG;

for ActualERP = min(IndexActualERP):max(IndexActualERP)
   
EEG  = pop_binlister( EEG_raw, 'BDF', sprintf('%sBinLister_%s.txt',BinListerDir, ERPsToExtract{1,ActualERP}),...
                      'ExportEL', sprintf('%sP%s_%s_Bins_%s.txt',...
                      EventListFilesDir, actualsbj_id, ExperimentName, ERPsToExtract{1,ActualERP}),...
                      'IndexEL',  1, 'SendEL2', 'All', 'Voutput', 'EEG' ); 
                  
%---------------------------------------------------------------------
% Create Bin-Based Epochs for each ERP component 
%---------------------------------------------------------------------     
EEG = pop_epochbin( EEG, str2num(EpochWindow{1,ActualERP}{1,1}),'none');

end                          
                                                      
%---------------------------------------------------------------------
% Save the dataset
%---------------------------------------------------------------------
ActualSetName = sprintf('P%s_%s_%s', actualsbj_id,...
                ExperimentName, ERPsToExtract{1,ActualERP});
mkdir(EpochedDataDir,sprintf('%s',ERPsToExtract{1,ActualERP}));
ActualERPDir = sprintf('%s%s\\',EpochedDataDir,ERPsToExtract{1,ActualERP}); 
pop_saveset( EEG(1),'filename',ActualSetName,'filepath',ActualERPDir);                           

end











