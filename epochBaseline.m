function [EEGlist, EpochList] = ...
    epochBaseline(EEG_table, baselinetype, expphase, epoch_T, epoch_dt)

if nargin < 5
    epoch_dt = [];
    if nargin < 4
        epoch_T = [];
        if nargin < 3
            expphase = ''; 
            if nargin < 2
                baselinetype = '';
            end
        end
    end
end

if isempty(epoch_T)
    epoch_T = 30; % s
end
if isempty(epoch_dt)
    epoch_dt = 30; % s
end
if isempty(expphase)
    expphase = 'before experiment';
end
if isempty(baselinetype)
    baselinetype = 'BaselineOpen';
end

if strcmpi(expphase, 'both')
    evalstr = ...
        ['[EEG_table.',baselinetype,'(''before experiment''), EEG_table.',...
         baselinetype,'(''after experiment'')];'];
else
    expphase = ['''',expphase,''''];
    evalstr = ['EEG_table.',baselinetype,'(',expphase,');'];
end
EEGlist = eval(evalstr); 
EEGlist = [EEGlist{:}]; 

EpochList = cell(size(EEGlist));
for lstIdx = 1:length(EEGlist)
    eeg = EEGlist(lstIdx);
    if ~isempty(eeg)
        t = (eeg.xmin):epoch_dt:((eeg.xmax)-epoch_T);
        curEpochs = repmat(eeg, size(t));
        for idx = 1:length(t)
            if ~mod(idx/length(t), .05)
                disp(['Epoch ',num2str(idx),' of ',num2str(length(t)),...
                    ' (',num2str(100*idx/length(t),3),'%)'])
            end
            curEpoch = pop_select(eeg, 'time', t(idx)+[0,epoch_T]);
            curEpoch.xmin = curEpoch.xmin + t(idx);
            curEpoch.xmax = curEpoch.xmax + t(idx);
            curEpoch.times = curEpoch.times + t(idx)*1000;
            curEpochs(idx) = curEpoch;
        end
        EpochList{lstIdx} = curEpochs;
    end
end

end