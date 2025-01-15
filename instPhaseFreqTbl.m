function [phiTbl, fTbl] = instPhaseFreqTbl(tbl)

fs = tbl.Properties.SampleRate; 
[phi_block, f_block] = instPhaseFreq(tbl.Variables, fs);

phiTbl = tbl; phiTbl.Variables = phi_block; 
for c = 1:width(phiTbl)
    phiTbl.Properties.VariableUnits{c} = 'Radians';
end

fTbl = tbl; fTbl.Variables = f_block; 
for c = 1:width(fTbl)
    fTbl.Properties.VariableUnits{c} = 'Hz';
end

end