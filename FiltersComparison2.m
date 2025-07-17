fc = [13, 30]; 
srate = 1000;
N = 4096; % freq points to sample

% FIR
FIRorder = 2*fix(srate/fc(1));
bFIR = fir1(FIRorder, fc/(.5*srate)); aFIR = 1;
[hFIR,wFIR] = freqz(bFIR, aFIR, N);

% Butterworth 
Border = 2;
[bB,aB] = butter(Border, fc/(.5*srate), "bandpass");
[hB,wB] = freqz(bB, aB, N);

% Cheby II
C2order = 2;
[bC2,aC2] = cheby2(C2order,5,fc/(.5*srate),"bandpass");
[hC2,wC2] = freqz(bC2, aC2, N);

% mag/phase response and group delay 
H = [hFIR, hB, hC2];
W = [wFIR, wB, wC2]; 
W = W*srate/(2*pi);
Hm = abs(H); Hp = angle(H);
GD = -diff(unwrap(Hp))./diff(W);
GD = GD/(2*pi);
WGD = .5*(W(2:end,:) + W(1:(end-1),:));

% evaluate 
Wsel = (WGD >= fc(1)) & (WGD <= fc(2));
for c = 1:width(GD)
    gdsel = GD(Wsel(:,c),c);
    [mean(gdsel) std(gdsel)]
end

% plot
figure; 
ax(1) = subplot(311); 
plot(W,Hm); ylabel('magnitude gain'); grid on; 
legend('FIR','Butter','Cheby2');
ax(2) = subplot(312); 
plot(W,Hp); ylabel('phase (rad)'); grid on;
ax(3) = subplot(313);
plot(WGD,GD); ylabel('group delay (s)'); grid on;
xlabel('freq (Hz)');
linkaxes(ax,'x');