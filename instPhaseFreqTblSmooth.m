function [phiTbl, fTbl] = instPhaseFreqTblSmooth(tbl, freqRng)
% 
% get instantaneous phase and frequency at each timepoint in a timetable 
% 

loco = freqRng(1); hico = freqRng(2);

fs = tbl.Properties.SampleRate; 
if isnan(fs)
    fs = 1/mean(seconds(diff(tbl.Time)));
end
[phi_block, f_block] = instPhaseFreq(tbl.Variables, fs);

% clop freq to range 
f_block(:) = min(hico, f_block(:));
f_block(:) = max(loco, f_block(:));

% smooth out deriv-based freq estimate 
Twin = (1/hico) * fs; % min-period-length window
Twin = floor(Twin);
f_block = smoothdata(f_block, "gaussian", Twin);

phiTbl = tbl; phiTbl.Variables = phi_block; 
for c = 1:width(phiTbl)
    phiTbl.Properties.VariableUnits{c} = 'Radians';
    phiTbl.Properties.VariableNames{c} = ...
        [tbl.Properties.VariableNames{c},' phase'];
end
phiTbl.Properties.Description = [tbl.Properties.Description,' phase'];

fTbl = tbl; fTbl.Variables = f_block; 
for c = 1:width(fTbl)
    fTbl.Properties.VariableUnits{c} = 'Hz';
    fTbl.Properties.VariableNames{c} = ...
        [tbl.Properties.VariableNames{c},' frequency'];
end
fTbl.Properties.Description = [tbl.Properties.Description,' frequency'];

end