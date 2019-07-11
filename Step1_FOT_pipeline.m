%---------------------------------------------------%
% FOT STUDY EEG DATA PROCESSING PIPELINE            %
%                                                   %
%1. Create EEGLAB .set files for your data files    %
%                                                   %
%2. Downsample to 500Hz                             %
%                                                   %
%3. Set your current directory (cd- line 29) to the %
%   folder containing your .set files               %
%                                                   %
%4. Use the filematALL (line 31) to set the pattern %
%   you want to use to look for files to split      %
%                                                   %
%4. Use the filematALLSplit (line 47) to set the    % 
%   pattern you want to use to look for files to    %
%   clean (these files should be in the             %
%   Split_Condition folder from above)              %
%                                                   %
%---------------------------------------------------%
%% PART 1: filter, create event list, assign bins, %
%          segment by condition(-200ms to 10000ms),%
%          baseline correct,                       %
%          and save a separate file per condition. %
%          The resulting files will be in          %
%          pathToFiles/Split_Condition             %
%--------------------------------------------------%

clear
cd '/Users/BCDLAB1600/Desktop/BEES Study/ADULT/Floating Objects/DATA_FILES/'
filematALL = dir('BEES_PRE_301_Adult_FOT.set'); % This loads a struct of files of a specific condition e.g. (Pre)    
filemat = {filematALL.name}'; % This takes the just the names from that struct and transposes the list so its in the correct format
dirFolder = cd('../');dirFolder = strcat(dirFolder, '/');pathToFiles = cd(dirFolder);
pathToFiles = strcat(pathToFiles,'/');

FOT_Split_Cond_batch(filemat, pathToFiles, dirFolder);

%--------------------------------------------------%
%  PART 2: outer band of electrodes removed (if 1) %
%          bad channel replacement and avg ref     %
%          resulting files saved in                %
%          pathToFiles/Split_Condition/CLEAN CHAN  %
%          A list of interpolated channels for each%
%          trial saved as interpvec_filename       %
%--------------------------------------------------%

cd ../; cd 'Split_Condition/'
filematALLSplit = dir('BEES_PRE_301_Adult_FOT_*.set'); % This loads a struct of files of a specific condition e.g. (Pre)    
filematSplit = {filematALLSplit.name}'; % This takes the just the names from that struct and transposes the list so its in the correct format
pathToFilesSplit= strcat(pathToFiles,'Split_Condition/');

%(filemat, path, Remove outer band? 0 or 1)
BEES_clean_data_batch(filematSplit, pathToFilesSplit,0);

result_path = strcat(pathToFilesSplit,'CLEAN CHAN/')
cd(result_path)