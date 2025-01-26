function yf = myFastForecastAR(arMdl, yPast, k)
% a faster version of forecast when arMdl is a no-input, discrete AR model
% with the same sampling rate as yPast, and yPast is a column vector/matrix. 

y = [yPast; nan(k,width(yPast))]; % [yPast; yf]

AA = arMdl.A; 
if ~iscell(AA)
    AA = {AA};
end

for chan = 1:height(AA)

    A = AA{chan,chan};
    A1 = A(1); A = A(2:end);
    N = length(A);
    
    for t = (height(yPast)+1):height(y)
        yp = y((t-1):-1:(t-N), chan);
        y(t,chan) = (-A*yp)/A1;
    end
    
end

yf = y((height(yPast)+1):end, :);

end