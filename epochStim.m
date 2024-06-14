function [EEGlist, EpochList] = epochStim(EEG_table, stimtype, epochT)

if nargin < 3
    epochT = [];
end
if isempty(epochT) | isnan(epochT)
    epochT = 1; % s
end

if strcmp(stimtype, 'PinPrick')
    stimnum = '11'; 
elseif strcmp(stimtype, 'TempStim')
    stimnum = '3';
end

evalstr = ...
    ['[EEG_table.',stimtype,'(''before experiment''), EEG_table.',stimtype,'(''after experiment'')];'];
EEGlist = eval(evalstr);
EEGlist = [EEGlist{:}]; 

    epochT = .5*epochT;

EpochList = cell(size(EEGlist));
for lstIdx = 1:length(EEGlist)
    eeg = EEGlist(lstIdx);
    if ~isempty(eeg)
        curEpochs3D = pop_epoch(eeg, {stimnum}, [-epochT, epochT]);
        t = 1:curEpochs3D.trials; 
        curEpochs = repmat(curEpochs3D, size(t));
        for idx = 1:length(t)
            curEpoch = curEpochs(idx); 
            curEpoch.trials = 1;
            curEpoch.data = curEpoch.data(:,:,idx);
            curEpochs(idx) = curEpoch;
        end
        EpochList{lstIdx} = curEpochs;
    end
end

end