function [phi_block, f_block] = instPhaseFreq(blockData, fs)
% Provided a block of data sampled at constant rate fs, calculate the
% instantaneous phase phi and frequency f using the Hilbert transform. 

% use hilbert transform to determine inst. phase
H_block = hilbert(blockData); 
phi_block = angle(H_block); 

%%{
% calculate inst. freq. 
f_block = gradient(unwrap(phi_block)')' *fs/(2*pi);
%f_block = abs(f_block); % ?? why does it go negative??
%}

%{
% calculate inst. freq. recursively 
w = H_block(2:end,:) .* conj(H_block(1:(end-1),:));
w = angle(w);
w = [w(1,:); w]; % pad first element 
f_block = fs*w/(2*pi);
%}

end