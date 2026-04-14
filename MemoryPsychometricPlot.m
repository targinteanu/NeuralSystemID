% get behavior file 
% 'sternberg_***_***.mat' that contains var 'correct'
[fn, fp] = uigetfile('sternberg*.mat');
load(fullfile(fp,fn), 'correct');

%% get experiment info 
% subj ID 'PY**N***'
subjID = '';
FP = upper(fp);
K = strfind(FP,'PY');
while length(K) > 0
    k = K(1);
    if (k+7 <= length(FP)) && (FP(k+4) == 'N')
        subjID = FP(k:k+7); 
        break;
    else
        K = k(2:end);
    end
end

%% make accuracy plot
nTrlAvg = 7;
fig1 = figure; 
plot(cumsum(correct)./(1:length(correct)), 'LineWidth',2);
hold on; grid on;
plot(correct, 'o', 'LineWidth',1);
plot(movmean(correct, nTrlAvg), ':', 'LineWidth',1);
xlabel('trial'); ylabel('accuracy'); 
legend('cumulative', 'ind trl', [num2str(nTrlAvg),'-trl avg'], ...
    'Location','best');
title([subjID,' \Psi -metrics']);

%% save
tosave = questdlg('Save Figure?', 'Figure Saving');
tosave = strcmp(tosave, 'Yes');
if tosave
    saveas(fig1, fullfile(fp,'CorrectCurve'),'png');
end