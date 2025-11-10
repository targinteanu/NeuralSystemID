function [fig, tl] = plotLTISS(sys)

A = sys.A; B = sys.B; C = sys.C; D = sys.D;
W1 = width(A); W2 = width(B); H1 = height(A); H2 = height(C);

% force min ratios to 1:N
N = 3; 
if W2/(W1+W2) < 1/N
    W2 = 1; W1 = N-1;
elseif W1/(W1+W2) < 1/N
    W1 = 1; W2 = N-1;
end
if H2/(H1+H2) < 1/N
    H2 = 1; H1 = N-1;
elseif H1/(H1+H2) < 1/N
    H1 = 1; H2 = N-1;
end

fig = figure; 
tl = tiledlayout(H1+H2, W1+W2);

% A
nexttile([H1, W1]); imagesc(A); title('A'); colorbar;
xlabel('state'); ylabel('state');
if numel(sys.StateName)
    if sum(cellfun(@numel,sys.StateName))
        xticks(1:length(sys.StateName));
        xticklabels(sys.StateName);
        yticks(1:length(sys.StateName));
        yticklabels(sys.StateName);
    end
end

% B
nexttile([H1, W2]); imagesc(B); title('B'); colorbar;
xlabel('input'); ylabel('state');
if numel(sys.InputName)
    if sum(cellfun(@numel,sys.InputName))
        xticks(1:length(sys.InputName));
        xticklabels(sys.InputName);
    end
end
if numel(sys.StateName)
    if sum(cellfun(@numel,sys.StateName))
        yticks(1:length(sys.StateName));
        yticklabels(sys.StateName);
    end
end

% C
nexttile([H2, W1]); imagesc(C); title('C'); colorbar;
xlabel('state'); ylabel('output');
if numel(sys.OutputName)
    if sum(cellfun(@numel,sys.OutputName))
        yticks(1:length(sys.OutputName));
        yticklabels(sys.OutputName);
    end
end
if numel(sys.StateName)
    if sum(cellfun(@numel,sys.StateName))
        xticks(1:length(sys.StateName));
        xticklabels(sys.StateName);
    end
end

% D
nexttile([H2, W2]); imagesc(D); title('D'); colorbar;
xlabel('input'); ylabel('output');
if numel(sys.OutputName)
    if sum(cellfun(@numel,sys.OutputName))
        yticks(1:length(sys.OutputName));
        yticklabels(sys.OutputName);
    end
end
if numel(sys.InputName)
    if sum(cellfun(@numel,sys.InputName))
        xticks(1:length(sys.InputName));
        xticklabels(sys.InputName);
    end
end

% title
if sys.Ts > 0
    % discrete 
    sgtitle(...
        {'\bf x[ \rm t+\Deltat \bf ] \rm = A \bf x[ \rm t \bf ] \rm + Bu[t]'; ...
        '\bf y[ \rm t \bf ] \rm = C \bf x[ \rm t \bf ] \rm + Du[t]'});
else
    % continuous 
    sgtitle(...
        {' ^d/_{dt} \bf x \rm = A \bf x \rm + Bu'; ...
        '\bf y \rm = C \bf x \rm + Du'});
end

end