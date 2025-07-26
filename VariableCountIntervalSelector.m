function [tRangesOut, commsOut] = VariableCountIntervalSelector(tRangesIn, commsIn)

if nargin < 1
    tRangesIn = [];
end
if nargin < 2
    commsIn = [];
end

if isempty(tRangesIn)
    tzone = '';
else
    tzone = tRangesIn.TimeZone;
end
if isempty(commsIn)
    commsIn = repmat("",size(tRangesIn));
end

tRangesOut = [];
commsOut = [];

fig = uifigure("WindowStyle","modal");
h0 = 300; 

donebtn = uibutton(fig, 'push', 'Position',[250,350,80,40], ...
    'Icon','success', 'Text','Done', ...
    'ButtonPushedFcn', @(src,evt) SaveTimesAndClose());

newpbtn = @(n) uibutton(fig, 'push', 'Position',[470,h0-(n-1)*60,20,20], ...
    'Text','+', 'ButtonPushedFcn', @(src,evt) addrow(n+1));
newmbtn = @(n) uibutton(fig, 'push', 'Position',[500,h0-(n-1)*60,20,20], ...
    'Text','-', 'ButtonPushedFcn', @(src,evt) remrow(n));

newsef = @(n, txt) uieditfield(fig, 'text', ...
    'Position',[100,h0-(n-1)*60,150,20], 'Value',txt);
neweef = @(n, txt) uieditfield(fig, 'text', ...
    'Position',[300,h0-(n-1)*60,150,20], 'Value',txt);

newcom = @(n, txt) uieditfield(fig, 'text', ...
    'Position',[100,h0-30-(n-1)*60,400,30], 'Value',txt);

btns = [newpbtn(1), newmbtn(1)]; 
etfs = [newsef(1,'NaT'), neweef(1,'NaT')];
coms = [];
if ~height(tRangesIn)
    tRangesIn = [NaT, NaT];
end
for n = 1:height(tRangesIn)
    btns(n,:) = [newpbtn(n), newmbtn(n)];
    etfs(n,:) = [newsef(n, string(tRangesIn(n,1))), ...
            neweef(n, string(tRangesIn(n,2)))];
    coms = [coms; newcom(n, commsIn(n))];
end

uiwait(fig);

function addrow(n)
if n > height(btns)
    btns(n,:) = [newpbtn(n), newmbtn(n)];
    etfs(n,:) = [newsef(n,'NaT'), neweef(n,'NaT')];
    coms = [coms; newcom(n, 'Enter comment here.')];
else
    for c = 1:width(btns)
        btns(n,c).Visible = true;
        etfs(n,c).Visible = true;
        coms(n,:).Visible = true;
    end
end
end

function remrow(n)
for c = 1:width(btns)
    btns(n,c).Visible = false;
    etfs(n,c).Visible = false;
    coms(n,:).Visible = false;
end
end

    function SaveTimesAndClose()
        for r = 1:height(etfs)
            if etfs(r,1).Visible && etfs(r,2).Visible
                tRangesOut = [tRangesOut; ...
                    datetime(etfs(r,1).Value, 'TimeZone',tzone), ...
                    datetime(etfs(r,2).Value, 'TimeZone',tzone)];
                commsOut = [commsOut; string(coms(r).Value)];
            end
        end
        pause(.01); drawnow; pause(.01);
        close(fig);
    end

end