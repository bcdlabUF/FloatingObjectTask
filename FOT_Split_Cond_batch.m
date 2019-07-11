function [] = FOT_Split_Cond_batch(filemat, pathToFiles, dirFolder)
    % check for the filt folder and create it if it doesn't
    % exist
    if ~exist(strcat(pathToFiles, 'filt/'),'dir')
        mkdir(strcat(pathToFiles, 'filt/'))
    end

for j = 1:size(filemat,1)
    subject_string = deblank(filemat(j,:));
    Csubject = char(subject_string);
    C = strsplit(Csubject,'.');
    subject = char(C(1,1));
    file = strcat(dirFolder,subject)
    filename = strcat(dirFolder,subject, '.set');
    
    EEG = pop_loadset('filename', filename);
    EEG = pop_editset(EEG, 'setname', strcat(file, '_FOT_chan'));

    % Create Event List
    EEG  = pop_creabasiceventlist( EEG , 'AlphanumericCleaning', 'on', ...
        'BoundaryNumeric', { -99 }, 'BoundaryString', { 'boundary' } );
    EEG = pop_editset(EEG, 'setname', strcat(file,'_FOT_chan_elist'));
    
    % Bandpass filter from 1-50 Hz
    % desired freq/250 = value 
    dataAK=double(EEG.data); 
    [alow, blow] = butter(6, 0.20); 
    [ahigh, bhigh] = butter(3,0.004, 'high'); 

    dataAKafterlow = filtfilt(alow, blow, dataAK'); 
    dataAKafterhigh = filtfilt(ahigh, bhigh, dataAKafterlow)'; 

    EEG.data = single(dataAKafterhigh); 
    
    EEG = pop_editset(EEG, 'setname', strcat(file,'_FOT_chan_elist_filt'));
    %save the filtered file
    EEG = pop_saveset(EEG, 'filename', strcat(pathToFiles, 'filt/', subject, '_1_50_filt.set'))
    
    % check for the Split_Condition folder and create it if it doesn't
    % exist
    if ~exist(strcat(pathToFiles, 'Split_Condition/'),'dir')
        mkdir(strcat(pathToFiles, 'Split_Condition/'))
    end
    
    NEWpath= strcat(pathToFiles, 'Split_Condition/');
    
    %Assign bins via binlister for each condition separately and save each 
    %condition to file separately    
    Condition = 'iu';
        EEG_IU  = pop_binlister( EEG , 'BDF', strcat(pathToFiles, ...
         '/bin_files/', 'iu.txt'), 'IndexEL',  1, 'SendEL2',...
         'EEG', 'Voutput', 'EEG' );

        EEG_IU = pop_editset(EEG_IU, 'setname', strcat(pathToFiles, ...
         'DATA FILES/',subject,'_chan_elist_filt_bins'));

        % Create bin-based epochs -200ms to 10000ms using -200to0 for
        % baseline correction
        EEG_IU = pop_epochbin( EEG_IU , [-200.0  10000.0],  'pre'); 
        EEG_IU = pop_editset(EEG_IU, 'setname', strcat(pathToFiles, ...
            'DATA FILES/',subject,'_chan_elist_filt_bins_be'));

        EEG_IU = pop_saveset( EEG_IU, 'filename',strcat(NEWpath,subject,...
            '_',Condition,'.set'));

    Condition = 'cu';

        EEG_CU  = pop_binlister( EEG , 'BDF', strcat(pathToFiles, ...
            '/bin_files/', 'cu.txt'), 'IndexEL',  1, 'SendEL2',...
            'EEG', 'Voutput', 'EEG' );

        EEG_CU = pop_editset(EEG_CU, 'setname', strcat(pathToFiles, ...
            'DATA FILES/',subject,'_chan_elist_filt_bins'));

        % Create bin-based epochs
        EEG_CU = pop_epochbin( EEG_CU , [-200.0  10000.0],  'pre'); %change this to 10,000 for new labels that are 10 seconds each.
        EEG_CU = pop_editset(EEG_CU, 'setname', strcat(pathToFiles, ...
            'DATA FILES/',subject,'_chan_elist_filt_bins_be'));
       
        EEG_CU = pop_saveset( EEG_CU, 'filename',strcat(NEWpath,subject,...
            '_',Condition,'.set'));
        
    Condition = 'un';

        EEG_UN  = pop_binlister( EEG , 'BDF', strcat(pathToFiles, ...
            '/bin_files/', 'un.txt'), 'IndexEL',  1, 'SendEL2',...
            'EEG', 'Voutput', 'EEG' );

        EEG_UN = pop_editset(EEG_UN, 'setname', strcat(pathToFiles, ...
            'DATA FILES/',subject,'_chan_elist_filt_bins'));

        % Create bin-based epochs
        EEG_UN = pop_epochbin( EEG_UN , [-200.0  10000.0],  'pre'); %change this to 10,000 for new labels that are 10 seconds each.
        EEG_UN = pop_editset(EEG_UN, 'setname', strcat(pathToFiles, ...
            'DATA FILES/',subject,'_chan_elist_filt_bins_be'));

        EEG_UN = pop_saveset( EEG_UN, 'filename',strcat(NEWpath,subject,...
            '_',Condition,'.set'));

    Condition = 'fa';

        EEG_FA  = pop_binlister( EEG , 'BDF', strcat(pathToFiles, ...
            '/bin_files/', 'fa.txt'), 'IndexEL',  1, 'SendEL2',...
            'EEG', 'Voutput', 'EEG' );

        EEG_FA = pop_editset(EEG_FA, 'setname', strcat(pathToFiles, ...
            'DATA FILES/',subject,'_chan_elist_filt_bins'));
            
        % Create bin-based epochs
        EEG_FA = pop_epochbin( EEG_FA , [-200.0  10000.0],  'pre'); %change this to 10,000 for new labels that are 10 seconds each.
        EEG_FA = pop_editset(EEG_FA, 'setname', strcat(pathToFiles, ...
            'DATA FILES/',subject,'_chan_elist_filt_bins_be'));
        
        EEG_FA = pop_saveset( EEG_FA, 'filename',strcat(NEWpath,subject,...
            '_',Condition,'.set'));

    Condition = 'fx';

        EEG_FX  = pop_binlister( EEG , 'BDF', strcat(pathToFiles, ...
            '/bin_files/', 'fx.txt'), 'IndexEL',  1, 'SendEL2',...
            'EEG', 'Voutput', 'EEG' );

        EEG_FX = pop_editset(EEG_FX, 'setname', strcat(pathToFiles, ...
            'DATA FILES/',subject,'_chan_elist_filt_bins'));
            
        % Create bin-based epochs 
        % THESE ARE SHORTER THAN THE OTHERS!!
        EEG_FX = pop_epochbin( EEG_FX , [-200.0  5000.0],  'pre'); %change this to 10,000 for new labels that are 10 seconds each.
        EEG_FX = pop_editset(EEG_FX, 'setname', strcat(pathToFiles, ...
            'DATA FILES/',subject,'_chan_elist_filt_bins_be'));
        
        EEG_FX = pop_saveset( EEG_FX, 'filename',strcat(NEWpath,subject,...
            '_',Condition,'.set'));

end
         
end

