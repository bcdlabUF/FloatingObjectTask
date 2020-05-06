% % % % % % % % % % % % % % % % % % % % % % % % % % %   
%                                                   % 
% FOT preprocessing pipeline edited based off MADE  %
%  3/26/2020 ZP                                     %
%                                                   %
% Enter your directory path on lines 17, 20, and 23 %
%                                                   %
% Enter file name on line 109                       %
%                                                   %
% Includes eventlist, removes outer band, filtering %
%  ICA, rejects bad ICs                             %
%                                                   %
% % % % % % % % % % % % % % % % % % % % % % % % % % %

%%

clear % clear matlab workspace
clc % clear matlab command window
%addpath(genpath('C:\Users\Berger\Documents\eeglab13_4_4b'));% enter the path of the EEGLAB folder in this line
%addpath(genpath('C:\Users\Berger\Documents\eeglab13_4_4b'))
 
% CondArray = {'task','fixation'};
% [selectionIndex3, leftBlank] = listdlg('PromptString', 'Do you wan to process the task trials or fixation trials?:', 'SelectionMode', 'single', 'ListString', CondArray);
% Cond= CondArray{selectionIndex3};
 
% 1. Enter the path of the folder that has the raw data to be analyzed
rawdata_location = '/Volumes/Hard Drive/BEES fot/All set/';
 
% 2. Enter the path of the folder where you want to save the processed data
output_location = '/Volumes/Hard Drive/BEES fot/MADE/Split Condition/';
 
% 3. Enter the path of the channel location file
channel_locations = ['/Volumes/Hard Drive/BEES fot' filesep 'locsEEGLAB129HCL.mat'];
 
% 4. Do your data need correction for anti-aliasing filter and/or task related time offset?
adjust_time_offset = 0; % 0 = NO (no correction), 1 = YES (correct time offset)
% If your data need correction for time offset, initialize the offset time (in milliseconds)
filter_timeoffset = 0;     % anti-aliasing time offset (in milliseconds). 0 = No time offset
stimulus_timeoffset   = 0; % stimulus related time offset (in milliseconds). 0 = No time offset
response_timeoffset = 0;    % response related time offset (in milliseconds). 0 = No time offset
stimulus_markers = {'xxx', 'xxx'};      % enter the stimulus makers that need to be adjusted for time offset
respose_markers = {'xxx', 'xxx'};       % enter the response makers that need to be adjusted for time offset
 
% 5. Do you want to down sample the data?
down_sample = 0; % 0 = NO (no down sampling), 1 = YES (down sampling)
sampling_rate = 500; % set sampling rate (in Hz), if you want to down sample
 
% 6. Do you want to delete the outer layer of the channels? (Rationale has been described in MADE manuscript)
delete_outerlayer = 1; % 0 = NO (do not delete outer layer), 1 = YES (delete outerlayer);
% If you want to delete outer layer, make a list of channels to be deleted
outerlayer_channel = {'E17' 'E43' 'E48' 'E49' 'E56' 'E63' 'E68' 'E73' 'E81' 'E88' 'E94' 'E99' 'E107' 'E113' 'E119' 'E120' 'E125' 'E126' 'E127' 'E128'}; % list of channels
% 7. Initialize the filters
highpass = 1; % High-pass frequency
lowpass  = 30; % Low-pass frequency. We recommend low-pass filter at/below line noise frequency (see manuscript for detail)
 
% 8. Are you processing task-related or resting-state EEG data?
% task_eeg = 1; % 0 = resting, 1 = task
% if strcmp(Cond,'task')
%     task_event_markers = { 'cu04' 'iu02' 'un05' 'fa06' }; % enter all the event/condition markers
% else
%     task_event_markers = { 'fx07'};
% end
%  inmat3d = EEG.data; 
%  % Go through data matrix trial by trial and identify noisy channels,
%                 % replace those electrodes with the average of the closest 6, and 
%                 % apply the average reference
%                 for trial = 1:size(inmat3d,3)
% 
%                     %Get the data for one trial
%                     trialdata2d = inmat3d(:, :, trial); 
 
% 9. Do you want to epoch/segment your data?
% epoch_data = 1; % 0 = NO (do not epoch), 1 = YES (epoch data)
% if strcmp(Cond,'task')
%     task_epoch_length = [-.2 10]; % epoch length in second
% else
%     task_epoch_length = [-.2 5];
% end
% rest_epoch_length = 0; % for resting EEG continuous data will be segmented into consecutive epochs of a specified length (here 2 second) by adding dummy events
% overlap_epoch = 0;     % 0 = NO (do not create overlapping epoch), 1 = YES (50% overlapping epoch)
% dummy_events ={'d10'}; % enter dummy events name
 
% 10. Do you want to remove/correct baseline?
remove_baseline = 1; % 0 = NO (no baseline correction), 1 = YES (baseline correction)
baseline_window = [-200 0]; % baseline period in milliseconds (MS) [] = entire epoch
 
% 11. Do you want to remove artifact laden epoch based on voltage threshold?
voltthres_rejection = 1; % 0 = NO, 1 = YES
volt_threshold = [-400 400]; % lower and upper threshold (in ?V)
 
% 12. Do you want to perform epoch level channel interpolation for artifact laden epoch? (see manuscript for detail)
interp_epoch = 1; % 0 = NO, 1 = YES.
frontal_channels = {'E21', 'E25', 'E14', 'E9', 'E15'}; % If you set interp_epoch = 1, enter the list of frontal channels to check (see manuscript for detail)
 
% %13. Do you want to interpolate the bad channels that were removed from data?
% interp_channels = 1; % 0 = NO (Do not interpolate), 1 = YES (interpolate missing channels)
%  
% % 14. Do you want to rereference your data?
% rerefer_data = 1; % 0 = NO, 1 = YES
% reref=[]; % Enter electrode name/s or number/s to be used for rereferencing
% For channel name/s enter, reref = {'channel_name', 'channel_name'};
% For channel name/s enter, reref = [channel_number, channel_number];
% For average rereference enter, reref = []; default is average rereference
 
% 15. Do you want to save interim results?
save_interim_result = 1; % 0 = NO (Do not save) 1 = YES (save interim results)
 
% 16. How do you want to save your data? .set or .mat
output_format = 1; % 1 = .set (EEGLAB data structure), 2 = .mat (Matlab data structure)
 
% ********* no need to edit beyond this point for EGI .mff data **********
% ********* for non-.mff data format edit data import function ***********
% ********* below using relevant data import plugin from EEGLAB **********
 
%% Read files to analyses
cd(rawdata_location)
datafile_names=dir('BEES_PRE_101_9_*.set');
datafile_names=datafile_names(~ismember({datafile_names.name},{'.','..'}));
datafile_names={datafile_names.name};
[filepath,name,ext] = fileparts(char(datafile_names{1}));
 
%% Check whether EEGLAB and all necessary plugins are in Matlab path.
if exist('eeglab','file')==0
    error(['Please make sure EEGLAB is on your Matlab path. Please see EEGLAB' ...
        'wiki page for download and instalation instructions']);
end
 
if strcmp(ext, '.mff')==1
     if exist('mff_import', 'file')==0
        error(['Please make sure "mffmatlabio" plugin is in EEGLAB plugin folder and on Matlab path.' ...
             ' Please see EEGLAB wiki page for download and instalation instructions of plugins.' ...
            ' If you are not analysing EGI .mff data, edit the data import function below.']);
     end
 else
     warning('Your data are not EGI .mff files. Make sure you edit data import function before using this script');
 end
 
if exist('pop_firws', 'file')==0
    error(['Please make sure  "firfilt" plugin is in EEGLAB plugin folder and on Matlab path.' ...
        ' Please see EEGLAB wiki page for download and instalation instructions of plugins.']);
end
 
if exist('channel_properties', 'file')==0
    error(['Please make sure "FASTER" plugin is in EEGLAB plugin folder and on Matlab path.' ...
        ' Please see EEGLAB wiki page for download and instalation instructions of plugins.']);
end
 
if exist('ADJUST', 'file')==0
    error(['Please make sure you download modified "ADJUST" plugin from GitHub (link is in MADE manuscript)' ...
        ' and ADJUST is in EEGLAB plugin folder and on Matlab path.']);
end
 
%% Create output folders to save data
if save_interim_result ==1
    if exist([output_location filesep 'filtered_data'], 'dir') == 0
        mkdir([output_location filesep 'filtered_data'])
    end
    if exist([output_location filesep 'ica_data'], 'dir') == 0
        mkdir([output_location filesep 'ica_data'])
    end
end
if exist([output_location filesep 'processed_data'], 'dir') == 0
    mkdir([output_location filesep 'processed_data'])
end
 
%% Initialize output variables
ica_preparation_bad_channels=[]; % number of bad channel/s due to channel/s exceeding xx% of artifacted epochs
length_ica_data=[]; % length of data (in second) fed into ICA decomposition
total_ICs=[]; % total independent components (ICs)
ICs_removed=[]; % number of artifacted ICs

%% Filter and ICA
for subject=1:length(datafile_names)
    EEG=[];
     
    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', subject, datafile_names{subject});
     
    %% STEP 1: Import EGI data file and relevant information
%     EEG=mff_import([rawdata_location filesep datafile_names{subject}]);
    EEG = pop_loadset([rawdata_location filesep datafile_names{subject}]);
    EEG = eeg_checkset(EEG);
     
    % Edit this data import function and use appropriate plugin from EEGLAB
    % for non-.mff data. For example, to import biosemi data, use biosig plugin.
    % The example codes for 64 channels biosemi data:
%     EEG = pop_biosig([rawdata_location, filesep, datafile_names{subject}]);
%     EEG = eeg_checkset(EEG);
%     EEG = pop_select( EEG,'nochannel', 65:72); % delete redundant channels
     
    %% STEP 2: Import channel locations -- already loaded if using a setfile
%     EEG=pop_chanedit(EEG, 'load',{channel_locations 'filetype' 'autodetect'});
%     EEG = eeg_checkset( EEG );
%    
%     % Check whether the channel locations were properly imported. The EEG signals and channel numbers should be same.
%     if size(EEG.data, 1) ~= length(EEG.chanlocs)
%         error('The size of the data does not match with channel numbers.');
%     end
     
    %% STEP 3: Adjust anti-aliasing and task related time offset
    if adjust_time_offset==1
        % adjust anti-aliasing filter time offset
        if filter_timeoffset~=0
            for aafto=1:length(EEG.event)
                EEG.event(aafto).latency=EEG.event(aafto).latency+(filter_timeoffset/1000)*EEG.srate;
            end
        end
        % adjust stimulus time offset
        if stimulus_timeoffset~=0
            for sto=1:length(EEG.event)
                for sm=1:length(stimulus_markers)
                    if strcmp(EEG.event(sto).type, stimulus_markers{sm})
                        EEG.event(sto).latency=EEG.event(sto).latency+(stimulus_timeoffset/1000)*EEG.srate;
                    end
                end
            end
        end
        % adjust response time offset
        if response_timeoffset~=0
            for rto=1:length(EEG.event)
                for rm=1:length(response_markers)
                    if strcmp(EEG.event(rto).type, response_markers{rm})
                        EEG.event(rto).latency=EEG.event(rto).latency-(response_timeoffset/1000)*EEG.srate;
                    end
                end
            end
        end
    end
     
    %% STEP 4: Change sampling rate -- already done if using .set file
%     if down_sample==1
%         if floor(sampling_rate) > EEG.srate
%             error ('Sampling rate cannot be higher than recorded sampling rate');
%         elseif floor(sampling_rate) ~= EEG.srate
%             EEG = pop_resample( EEG, sampling_rate);
%             EEG = eeg_checkset( EEG );
%         end
%     end
%% eventlist
    EEG  = pop_creabasiceventlist( EEG , 'AlphanumericCleaning', 'on', ...
        'BoundaryNumeric', { -99 }, 'BoundaryString', { 'boundary' } );
    EEG = pop_editset(EEG, 'setname', strcat(datafile_names,'_fot_chan_elist'));  
    %% STEP 5: Delete outer layer of channels
    chans_labels=cell(1,EEG.nbchan);
    for i=1:EEG.nbchan
        chans_labels{i}= EEG.chanlocs(i).labels;
    end
    [chans,chansidx] = ismember(outerlayer_channel, chans_labels);
    outerlayer_channel_idx = chansidx(chansidx ~= 0);
    if delete_outerlayer==1
        if isempty(outerlayer_channel_idx)==1
            error(['None of the outer layer channels present in channel locations of data.'...
                ' Make sure outer layer channels are present in channel labels of data (EEG.chanlocs.labels).']);
        else
            EEG = pop_select( EEG,'nochannel', outerlayer_channel_idx);
            EEG = eeg_checkset( EEG );
        end
    end
     
    %% STEP 6: Filter data
    % Calculate filter order using the formula: m = dF / (df / fs), where m = filter order,
    % df = transition band width, dF = normalized transition width, fs = sampling rate
    % dF is specific for the window type. Hamming window dF = 3.3
     
    high_transband = highpass; % high pass transition band
    low_transband = 10; % low pass transition band
     
    hp_fl_order = 3.3 / (high_transband / EEG.srate);
    lp_fl_order = 3.3 / (low_transband / EEG.srate);
     
    % Round filter order to next higher even integer. Filter order is always even integer.
    if mod(floor(hp_fl_order),2) == 0
        hp_fl_order=floor(hp_fl_order);
    elseif mod(floor(hp_fl_order),2) == 1
        hp_fl_order=floor(hp_fl_order)+1;
    end
     
    if mod(floor(lp_fl_order),2) == 0
        lp_fl_order=floor(lp_fl_order)+2;
    elseif mod(floor(lp_fl_order),2) == 1
        lp_fl_order=floor(lp_fl_order)+1;
    end
     
    % Calculate cutoff frequency
    high_cutoff = highpass/2;
    low_cutoff = lowpass + (low_transband/2);
     
    % Performing high pass filtering
    EEG = eeg_checkset( EEG );
    EEG = pop_firws(EEG, 'fcutoff', high_cutoff, 'ftype', 'highpass', 'wtype', 'hamming', 'forder', hp_fl_order, 'minphase', 0);
    EEG = eeg_checkset( EEG );
     
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
     
    % pop_firws() - filter window type hamming ('wtype', 'hamming')
    % pop_firws() - applying zero-phase (non-causal) filter ('minphase', 0)
     
    % Performing low pass filtering
    EEG = eeg_checkset( EEG );
    EEG = pop_firws(EEG, 'fcutoff', low_cutoff, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder', lp_fl_order, 'minphase', 0);
    EEG = eeg_checkset( EEG );
     
    % pop_firws() - transition band width: 10 Hz
    % pop_firws() - filter window type hamming ('wtype', 'hamming')
    % pop_firws() - applying zero-phase (non-causal) filter ('minphase', 0)
     
    %% STEP 7: Run faster to find bad channels
    % First check whether reference channel (i.e. zeroed channels) is present in data
    % reference channel is needed to run faster
%     ref_chan=[]; FASTbadChans=[]; all_chan_bad_FAST=0;
%     ref_chan=find(any(EEG.data, 2)==0);
%     if numel(ref_chan)>1
%         error(['There are more than 1 zeroed channel (i.e. zero value throughout recording) in data.'...
%             ' Only reference channel should be zeroed channel. Delete the zeroed channel/s which is not reference channel.']);
%     elseif numel(ref_chan)==1
%         list_properties = channel_properties(EEG, 1:EEG.nbchan, ref_chan); % run faster
%         FASTbadIdx=min_z(list_properties);
%         FASTbadChans=find(FASTbadIdx==1);
%         FASTbadChans=FASTbadChans(FASTbadChans~=ref_chan);
%         reference_used_for_faster{subject}={EEG.chanlocs(ref_chan).labels};
%         EEG = pop_select( EEG,'nochannel', ref_chan);
%         EEG = eeg_checkset(EEG);
%         channels_analysed=EEG.chanlocs; % keep full channel locations to use later for interpolation of bad channels
%     elseif numel(ref_chan)==0
%         warning('Reference channel is not present in data. Cz channel will be used as reference channel.');
%         ref_chan=find(strcmp({EEG.chanlocs.labels}, 'Cz')); % find Cz channel index
        EEG_copy=[];
        EEG_copy=EEG; % make a copy of the dataset
%         EEG_copy = pop_reref( EEG_copy, ref_chan,'keepref','on'); % rerefer to Cz in copied dataset
        EEG_copy = eeg_checkset(EEG_copy);
%         list_properties = channel_properties(EEG_copy, 1:EEG_copy.nbchan, ref_chan); % run faster on copied dataset
%         FASTbadIdx=min_z(list_properties);
%         FASTbadChans=find(FASTbadIdx==1);
%         channels_analysed=EEG.chanlocs;
%         reference_used_for_faster{subject}={EEG.chanlocs(ref_chan).labels};
%     end
%      
%     % If FASTER identifies all channels as bad channels, save the dataset
%     % at this stage and ignore the remaining of the preprocessing.
%     if numel(FASTbadChans)==EEG.nbchan || numel(FASTbadChans)+1==EEG.nbchan
%         all_chan_bad_FAST=1;
%         warning(['No usable data for datafile', datafile_names{subject}]);
%         if output_format==1
%             EEG = eeg_checkset(EEG);
%             EEG = pop_editset(EEG, 'setname',  strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_channels'));
%             EEG = pop_saveset(EEG, 'filename', strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_channels.set'),'filepath', [output_location filesep 'processed_data' filesep ]); % save .set format
%         elseif output_format==2
%             save([[output_location filesep 'processed_data' filesep ] strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_channels.mat')], 'EEG'); % save .mat format
%         end
%     else
%         % Reject channels that are bad as identified by Faster
%         EEG = pop_select( EEG,'nochannel', FASTbadChans);
%         EEG = eeg_checkset(EEG);
%     end
%      
%     if numel(FASTbadChans)==0
%         faster_bad_channels{subject}='0';
%     else
%         faster_bad_channels{subject}=num2str(FASTbadChans');
%     end
%      
%     if all_chan_bad_FAST==1
%         faster_bad_channels{subject}='0';
%         ica_preparation_bad_channels{subject}='0';
%         length_ica_data(subject)=0;
%         total_ICs(subject)=0;
%         ICs_removed{subject}='0';
%         total_epochs_before_artifact_rejection=0;
%         total_epochs_after_artifact_rejection=0;
%         total_channels_interpolated=0;
%         continue % ignore rest of the processing and go to next subject
%     end
%      
    %% Save data after running filter and FASTER function, if saving interim results was preferred
    if save_interim_result ==1
        if output_format==1
            EEG = eeg_checkset( EEG );
            EEG = pop_editset(EEG, 'setname', strrep(datafile_names{subject}, ext, '_filtered_data'));
            EEG = pop_saveset( EEG,'filename',strrep(datafile_names{subject}, ext, '_filtered_data.set'),'filepath', [output_location filesep 'filtered_data' filesep]); % save .set format
        elseif output_format==2
            save([[output_location filesep 'filtered_data' filesep ] strrep(datafile_names{subject}, ext, '_filtered_data.mat')], 'EEG'); % save .mat format
        end
    end
     
    %% STEP 8: Prepare data for ICA
    EEG_copy=[];
    EEG_copy=EEG; % make a copy of the dataset
    EEG_copy = eeg_checkset(EEG_copy);
     
    % Perform 1Hz high pass filter on copied dataset
    transband = 1;
    fl_cutoff = transband/2;
    fl_order = 3.3 / (transband / EEG.srate);
     
    if mod(floor(fl_order),2) == 0
        fl_order=floor(fl_order);
    elseif mod(floor(fl_order),2) == 1
        fl_order=floor(fl_order)+1;
    end
     
    EEG_copy = pop_firws(EEG_copy, 'fcutoff', fl_cutoff, 'ftype', 'highpass', 'wtype', 'hamming', 'forder', fl_order, 'minphase', 0);
    EEG_copy = eeg_checkset(EEG_copy);
     
    % Create 1 second epoch
    EEG_copy=eeg_regepochs(EEG_copy,'recurrence', 1, 'limits',[0 1], 'rmbase', [NaN], 'eventtype', '999'); % insert temporary marker 1 second apart and create epochs
    EEG_copy = eeg_checkset(EEG_copy);
    
     
    % Find bad epochs and delete them from dataset
    vol_thrs = [-1000 1000]; % [lower upper] threshold limit(s) in mV.
    emg_thrs = [-100 30]; % [lower upper] threshold limit(s) in dB.
    emg_freqs_limit = [20 40]; % [lower upper] frequency limit(s) in Hz.
     
    % Find channel/s with xx% of artifacted 1-second epochs and delete them
    chanCounter = 1; ica_prep_badChans = [];
    numEpochs =EEG_copy.trials; % find the number of epochs
    all_bad_channels=0;
     
    for ch=1:EEG_copy.nbchan
        % Find artifaceted epochs by detecting outlier voltage
        EEG_copy = pop_eegthresh(EEG_copy,1, ch, vol_thrs(1), vol_thrs(2), EEG_copy.xmin, EEG_copy.xmax, 0, 0);
        EEG_copy = eeg_checkset( EEG_copy );
         
        % 1         : data type (1: electrode, 0: component)
        % 0         : display with previously marked rejections? (0: no, 1: yes)
        % 0         : reject marked trials? (0: no (but store the  marks), 1:yes)
         
        % Find artifaceted epochs by using thresholding of frequencies in the data.
        % this method mainly rejects muscle movement (EMG) artifacts
        EEG_copy = pop_rejspec( EEG_copy, 1,'elecrange',ch ,'method','fft','threshold', emg_thrs, 'freqlimits', emg_freqs_limit, 'eegplotplotallrej', 0, 'eegplotreject', 0);
         
        % method                : method to compute spectrum (fft)
        % threshold             : [lower upper] threshold limit(s) in dB.
        % freqlimits            : [lower upper] frequency limit(s) in Hz.
        % eegplotplotallrej     : 0 = Do not superpose rejection marks on previous marks stored in the dataset.
        % eegplotreject         : 0 = Do not reject marked trials (but store the  marks).
         
        % Find number of artifacted epochs
        EEG_copy = eeg_checkset( EEG_copy );
        EEG_copy = eeg_rejsuperpose( EEG_copy, 1, 1, 1, 1, 1, 1, 1, 1);
        artifacted_epochs=EEG_copy.reject.rejglobal;
         
        % Find bad channel / channel with more than 20% artifacted epochs
        if sum(artifacted_epochs) > (numEpochs*20/100)
            ica_prep_badChans(chanCounter) = ch;
            chanCounter=chanCounter+1;
        end
    end
     
    % If all channels are bad, save the dataset at this stage and ignore the remaining of the preprocessing.
    if numel(ica_prep_badChans)==EEG.nbchan || numel(ica_prep_badChans)+1==EEG.nbchan
        all_bad_channels=1;
        warning(['No usable data for datafile', datafile_names{subject}]);
        if output_format==1
            EEG = eeg_checkset(EEG);
            EEG = pop_editset(EEG, 'setname',  strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_channels'));
            EEG = pop_saveset(EEG, 'filename', strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_channels.set'),'filepath', [output_location filesep 'processed_data' filesep ]); % save .set format
        elseif output_format==2
            save([[output_location filesep 'processed_data' filesep ] strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_channels.mat')], 'EEG'); % save .mat format
        end
         
    else
        % Reject bad channel - channel with more than xx% artifacted epochs
        EEG_copy = pop_select( EEG_copy,'nochannel', ica_prep_badChans);
        EEG_copy = eeg_checkset(EEG_copy);
    end
     
    if numel(ica_prep_badChans)==0
        ica_preparation_bad_channels{subject}='0';
    else
        ica_preparation_bad_channels{subject}=num2str(ica_prep_badChans);
    end
     
    if all_bad_channels == 1
        length_ica_data(subject)=0;
        total_ICs(subject)=0;
        ICs_removed{subject}='0';
        total_epochs_before_artifact_rejection=0;
        total_epochs_after_artifact_rejection=0;
        total_channels_interpolated=0;
        continue % ignore rest of the processing and go to next datafile
    end
     
    % Find the artifacted epochs across all channels and reject them before doing ICA.
    EEG_copy = pop_eegthresh(EEG_copy,1, 1:EEG_copy.nbchan, vol_thrs(1), vol_thrs(2), EEG_copy.xmin, EEG_copy.xmax,1,0);
    EEG_copy = eeg_checkset(EEG_copy);
     
    % 1         : data type (1: electrode, 0: component)
    % 0         : display with previously marked rejections? (0: no, 1: yes)
    % 0         : reject marked trials? (0: no (but store the  marks), 1:yes)
     
    % Find artifaceted epochs by using power threshold in 20-40Hz frequency band.
    % This method mainly rejects muscle movement (EMG) artifacts.
    EEG_copy = pop_rejspec(EEG_copy, 1,'elecrange', 1:EEG_copy.nbchan, 'method', 'fft', 'threshold', emg_thrs ,'freqlimits', emg_freqs_limit, 'eegplotplotallrej', 0, 'eegplotreject', 1);
     
    % method                : method to compute spectrum (fft)
    % threshold             : [lower upper] threshold limit(s) in dB.
    % freqlimits            : [lower upper] frequency limit(s) in Hz.
    % eegplotplotallrej     : 0 = Do not superpose rejection marks on previous marks stored in the dataset.
    % eegplotreject         : 0 = Do not reject marked trials (but store the  marks).
     
    % Find the number of artifacted epochs and reject them
    EEG_copy = eeg_checkset(EEG_copy);
    EEG_copy = eeg_rejsuperpose(EEG_copy, 1, 1, 1, 1, 1, 1, 1, 1);
    reject_artifacted_epochs=EEG_copy.reject.rejglobal;
    EEG_copy = pop_rejepoch(EEG_copy, reject_artifacted_epochs, 0);
     
    %% STEP 9: Run ICA
    length_ica_data(subject)=EEG_copy.trials; % length of data (in second) fed into ICA
    EEG_copy = eeg_checkset(EEG_copy);
    EEG_copy = pop_runica(EEG_copy, 'icatype', 'runica', 'extended', 1, 'stop', 1E-7, 'interupt','off');
     
    % Find the ICA weights that would be transferred to the original dataset
    ICA_WINV=EEG_copy.icawinv;
    ICA_SPHERE=EEG_copy.icasphere;
    ICA_WEIGHTS=EEG_copy.icaweights;
    ICA_CHANSIND=EEG_copy.icachansind;
     
    % If channels were removed from copied dataset during preparation of ica, then remove
    % those channels from original dataset as well before transferring ica weights.
    EEG = eeg_checkset(EEG);
    EEG = pop_select(EEG,'nochannel', ica_prep_badChans);
     
    % Transfer the ICA weights of the copied dataset to the original dataset
    EEG.icawinv=ICA_WINV;
    EEG.icasphere=ICA_SPHERE;
    EEG.icaweights=ICA_WEIGHTS;
    EEG.icachansind=ICA_CHANSIND;
    EEG = eeg_checkset(EEG);
     
    %% STEP 10: Run adjust to find artifacted ICA components
    badICs=[]; EEG_copy =[];
    EEG_copy = EEG;
    EEG_copy =eeg_regepochs(EEG_copy,'recurrence', 1, 'limits',[0 1], 'rmbase', [NaN], 'eventtype', '999'); % insert temporary marker 1 second apart and create epochs
    EEG_copy = eeg_checkset(EEG_copy);
     
%     if save_interim_result==1
%             badICs = adjusted_ADJUST(EEG_copy, [[output_location filesep 'ica_data' filesep] strrep(datafile_names{subject}, ext, '_adjust_report')]);
%     else
%             badICs = adjusted_ADJUST(EEG_copy, [[output_location filesep 'processed_data' filesep] strrep(datafile_names{subject}, ext, '_adjust_report')]);
%     end
    close all;
        
    % Mark the bad ICs found by ADJUST
    for ic=1:length(badICs)
        EEG.reject.gcompreject(1, badICs(ic))=1;
        EEG = eeg_checkset(EEG);
    end
    total_ICs(subject)=size(EEG.icasphere, 1);
    if numel(badICs)==0
        ICs_removed{subject}='0';
    else
        ICs_removed{subject}=num2str(double(badICs));
    end
     
    %% Save dataset after ICA, if saving interim results was preferred
    if save_interim_result==1
        if output_format==1
            EEG = eeg_checkset(EEG);
            EEG = pop_editset(EEG, 'setname',  strrep(datafile_names{subject}, ext, '_ica_data'));
            EEG = pop_saveset(EEG, 'filename', strrep(datafile_names{subject}, ext, '_ica_data.set'),'filepath', [output_location filesep 'ica_data' filesep ]); % save .set format
        elseif output_format==2
            save([[output_location filesep 'ica_data' filesep ] strrep(datafile_names{subject}, ext, '_ica_data.mat')], 'EEG'); % save .mat format
        end
    end
     
    %% STEP 11: Remove artifacted ICA components from data
    all_bad_ICs=0;
    ICs2remove=find(EEG.reject.gcompreject); % find ICs to remove
     
    % If all ICs and bad, save data at this stage and ignore rest of the preprocessing for this subject.
    if numel(ICs2remove)==total_ICs(subject)
        all_bad_ICs=1;
        warning(['No usable data for datafile', datafile_names{subject}]);
        if output_format==1
            EEG = eeg_checkset(EEG);
            EEG = pop_editset(EEG, 'setname',  strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_ICs'));
            EEG = pop_saveset(EEG, 'filename', strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_ICs.set'),'filepath', [output_location filesep 'processed_data' filesep ]); % save .set format
        elseif output_format==2
            save([[output_location filesep 'processed_data' filesep ] strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_ICs.mat')], 'EEG'); % save .mat format
        end
    else
        EEG = eeg_checkset( EEG );
        EEG = pop_subcomp( EEG, ICs2remove, 0); % remove ICs from dataset
    end

    if all_bad_ICs==1
        total_epochs_before_artifact_rejection=0;
        total_epochs_after_artifact_rejection=0;
        total_channels_interpolated=0;
%         continue % ignore rest of the processing and go to next datafile
    end
end
