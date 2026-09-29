function MIs = calcPACheatmap(x, phiBnd, ampBnd, fs, showheatmap)

if nargin < 5
    showheatmap = false;
end
if nargin < 4
    fs = [];
end
if isempty(fs) || isnan(fs)
    fs = 2; % freqs will be normalized
    frequnit = '';
else
    frequnit = ' (hz)';
end

if (numel(size(x))>2) || (~any(size(x)==1))
    error('x should be single channel.')
end
x = x(:); % make column 
L = length(x);
xPhi = repmat(x,1,length(phiBnd)-1);
xAmp = repmat(x,1,length(ampBnd)-1);

% filter to phase/amp ranges
disp('Filtering...')
for phi = 1:width(xPhi)
    BPF = buildFIRBPF(fs, phiBnd(phi),phiBnd(phi+1), 4, 1023);
    xPhi(:,phi) = filtfilt(BPF,1, xPhi(:,phi));
end
for amp = 1:width(xAmp)
    BPF = buildFIRBPF(fs, ampBnd(amp),ampBnd(amp+1), 4, 1023);
    xAmp(:,amp) = filtfilt(BPF,1, xAmp(:,amp));
end

% focus on the central 80% of data to avoid edge effects
t1 = floor(.1*L); t2 = ceil(.9*L);

% calc PAC 
disp('Calculating PAC...')
xAmp = abs(hilbert(xAmp));   xAmp = xAmp(t1:t2,:);
xPhi = angle(hilbert(xPhi)); xPhi = xPhi(t1:t2,:);
Amp = (ampBnd(2:end) + ampBnd(1:(end-1)))*.5;
Phi = (phiBnd(2:end) + phiBnd(1:(end-1)))*.5;
MIs = nan(width(xAmp), width(xPhi));
for amp = 1:width(xAmp)
    for phi = 1:width(xPhi)
        MIs(amp,phi) = calcPAChelper(xPhi(:,phi), xAmp(:,amp), [], false);
    end
end

if showheatmap
    figure; 
    img = imagesc(Phi,Amp,MIs); colorbar; 
    title('PAC MI'); 
    xlabel(['phase frequency',frequnit]);
    ylabel(['amp frequency',frequnit]);
    img.Parent.YDir = 'normal';
end

end