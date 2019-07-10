%--------------------------------------------------%
% FOT STUDY EEG DATA PROCESSING PIPELINE           %
%                                                  %
%1. Create EEGLAB .set files for your data files   %
%                                                  %
%2. Downsample to 500Hz                            %
%                                                  %
%                                                  %
%3. Use the filematALL (line 33) to set the pattern%
%   you want to use to look for files and both cd  %
%   (line 32) &pathToFiles (line 35) to specify the%
%   directory                                      %
%                                                  %
%4. Use the filematALLSplit (line 48) to set the   % 
%   pattern you want to use to look for files and  %
%   both cd (line 47) & pathToFilesSplit (line 50) %
%   to specify the directory. (these files should  %
%   be in the Split_Condition folder from above)   %
%                                                  %
%--------------------------------------------------%
%% PART 1: filter, create event list, assign bins, %
%          segment by condition(-200ms to 10000ms),%
%          baseline correct,                       %
%          and save a separate file per condition. %
%          The resulting files will be in          %
%          pathToFiles/Split_Condition             %
%   Run Section to run PARTS 1&2                   %
%--------------------------------------------------%

clear
cd '/Users/BCDLAB1600/Desktop/BEES Study/ADULT/Floating Objects/DATA_FILES/'
filematALL = dir('BEES_POST_3*_Adult_FOT.set'); % This loads a struct of files of a specific condition e.g. (Pre)    
filemat = {filematALL.name}'; % This takes the just the names from that struct and transposes the list so its in the correct format
pathToFiles='/Users/BCDLAB1600/Desktop/BEES Study/ADULT/Floating Objects/'

FOT_Split_Cond_batch(filemat, pathToFiles);

%--------------------------------------------------%
%  PART 2: outer band of electrodes removed (if 1) %
%          bad channel replacement and avg ref     %
%          resulting files saved in                %
%          pathToFiles/Split_Condition/CLEAN CHAN  %
%          A list of interpolated channels for each%
%          trial saved as interpvec_filename       %
%--------------------------------------------------%

cd '/Users/BCDLAB1600/Desktop/BEES Study/ADULT/Floating Objects/Split_Condition/'
filematALLSplit = dir('BEES_POST_3*_Adult_FOT_*.set'); % This loads a struct of files of a specific condition e.g. (Pre)    
filematSplit = {filematALLSplit.name}'; % This takes the just the names from that struct and transposes the list so its in the correct format
pathToFilesSplit='/Users/BCDLAB1600/Desktop/BEES Study/ADULT/Floating Objects/Split_Condition/'

%(filemat, path, Remove outer band? 0 or 1)
BEES_clean_data_batch(filematSplit, pathToFilesSplit,0);

cd '/Users/BCDLAB1600/Desktop/BEES Study/ADULT/Floating Objects/Split_Condition/CLEAN CHAN/'