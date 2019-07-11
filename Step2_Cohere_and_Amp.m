%% Compute coherence and powmean
%-------------------------------------------------------------------%
%   Enter the pattern you want to use for searching on line17       %
%                                                                   %
%   This script divides each trial into 1sec segments and flags     %
%   and then removes segments greater than 1.5SD from median        %
%                                                                   %
%   Coherence is then computed across the 129 channels from 1-20Hz  %
%   And this file is saved as file_cohere.mat                       %
%                                                                   %
%   Amplitude is then computed and saved as file_pow.mat            %
%                                                                   %
%-------------------------------------------------------------------%

clear
clc
filematALL = dir('BEES_PRE_349_Adult_FOT_*_CLEAN.set'); % This loads a struct of files of a specific condition e.g. (Pre)    
filemat = {filematALL.name}'; % This takes the just the names from that struct and transposes the list so its in the correct format
pathToFiles = strcat(cd,'/');

 for j = 1:size(filemat,1)
    subject_string = deblank(filemat(j,:));
    Csubject = char(subject_string);
    C = strsplit(Csubject,'.');
    sub = char(C(1,1));
    D = strsplit(sub,'_');
    cond = char(D(1,6));
    filename = strcat(pathToFiles,Csubject);
    % load in a file
    EEG = pop_loadset('filename',filename);
    
    if cond== 'iu' | cond == 'cu' | cond== 'un' | cond == 'fa'
        % remove the baseline 
        data = EEG.data(:, 101:5100,:);
    elseif cond == 'fx'
        % remove the baseline-- fx trials are shorter 
        data = EEG.data(:, 101:2600,:);
    end
    
    % remove bad segments in each trial
    for trial=1:size(data,3);    
        % cut the trial into 5 (or 10) 1s segments
        if cond == 'fx'
            datamat = double(reshape(data(:,:,trial), [size(data,1) 500 5])); 
        else
            datamat = double(reshape(data(:,:,trial), [size(data,1) 500 10])); 
        end
        
        stdvec = squeeze(median(std(datamat)));
        % find segments greater than 1.5 stdev from the median
        badindex = find(stdvec > median(stdvec) .* 1.5) ;

        % remove those bad trials
        datamat(:, :, badindex) = [];

        disp('segments that are bad: ')
        disp(length(badindex))
        % assign data to a new variable
        v = genvarname(strcat('datamat',num2str(trial)));    
        eval([v '= datamat;']);
    end
    
    if size(data,3) == 2 % if there are 2 trials, stack the segments
        datamat=cat(3,datamat1,datamat2);
    elseif size(data,3) == 1 % if there's only 1 trial, use the segments from that one trial
        datamat=datamat1;
    elseif size(data,3) == 3
        datamat=cat(3,datamat1,datamat2,datamat3);
    elseif size(data,3) == 4
        datamat=cat(3,datamat1,datamat2,datamat3, datamat4);
    elseif size(data,3) > 4
        datamat=cat(3,datamat1,datamat2,datamat3, datamat4, datamat5);
    end
    
    %gives coherence from 0-20Hz as channels x channels x freq (129x129x21)
    tic
    [outmat, F] = cohere_baby(datamat);
    t = toc
    
    %calculates amplitude for each segment and creates and average
    for seg = 1:size(datamat,3); 
        if seg == 1; [powsum, phase, freqs] = FFT_spectrum(squeeze(datamat(:, :, seg)), 500);
        else powsum = powsum + FFT_spectrum(squeeze(datamat(:, :, seg)), 500);
        end
    end
    powmean = powsum./size(datamat,3);
    
    %save the files
    save(strcat(sub, '_pow.mat'), 'powmean');
    save(strcat(sub, '_cohere.mat'), 'outmat');
 end
 
