%% load the data 
[fn,fp] = uigetfile('*_DataTimeTables.mat');
load(fullfile(fp,fn));
disp(fn);
disp(DataTimeTables(:,[1,3,4])); 

%% get duration of baseline 
% baseline defined as initial baseline (1) only
disp(DataTimeTables{1,1})
disp('Duration of baseline recording is ')
disp(DataTimeTables{1,4} - DataTimeTables{1,3})

%% get stimulation info 
% stimulation defined as cortical only
disp(DataTimeTables{3,1})
disp('Duration of cortical stimulation is ')
disp(DataTimeTables{3,4} - DataTimeTables{3,3})

stimtbl = DataTimeTables{3,2};
stimevt = stimtbl.Properties.Events;
stimevt = stimevt(contains(lower(stimevt.EventLabels), 'stim'), :);
stimTime = stimevt.Time;

isi = seconds(diff(stimTime));
disp('inter-stim-interval min, median, mean, max (s): ')
disp([min(isi), median(isi), mean(isi), max(isi)])

[~,splidx] = max(isi);
splTime1 = stimTime(splidx); splTime2 = stimTime(splidx+1);
splTime = mean([splTime1, splTime2]);
disp('split into periods of duration: ')
disp([splTime-stimtbl.Time(1), stimtbl.Time(end)-splTime])