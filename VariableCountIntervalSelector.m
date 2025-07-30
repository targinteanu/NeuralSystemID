function [tRangesOut, commsOut] = VariableCountIntervalSelector(tRangesIn, commsIn)
% Interactively select time ranges and assign comments in a modal window.
% Time ranges can be added, removed, and edited. 
% Starting values of time ranges (tRangesIn) and comments (commsIn) will
% determine how the window looks on startup. Time ranges are lists of time
% periods stacked vertically with [start, end] time horizontally. Comments
% are a vertical string array. 
% Values, whether input on startup or specified by user, will only be
% accepted if confirmed by the user. If the window is closed, output will
% be empty. 

%% handle input args 

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
    commsIn = repmat("",height(tRangesIn),1);
end

% handle overflowing inputs 
% TO DO: this is a bandaid fix. There should be a better way to handle
% this. 
maxNumRows = 5;
if height(tRangesIn) > maxNumRows
    % select the longest durations 
    durs = seconds(diff(tRangesIn,[],2));
    [~,longdurs] = sort(durs,'descend');
    longdurs = longdurs(1:maxNumRows);
    tRangesIn = tRangesIn(longdurs,:); commsIn = commsIn(longdurs,:);
end

%% setup 
% setup variables and (rules for new) buttons/UI elements 

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

% either initialize or apply input values 
btns = [newpbtn(1), newmbtn(1)]; 
etfs = [newsef(1,'NaT'), neweef(1,'NaT')];
coms = [];
if ~height(tRangesIn)
    tRangesIn = [NaT, NaT];
end
if ~height(commsIn)
    commsIn = "Enter comment here.";
end
for n = 1:height(tRangesIn)
    btns(n,:) = [newpbtn(n), newmbtn(n)];
    etfs(n,:) = [newsef(n, string(tRangesIn(n,1))), ...
            neweef(n, string(tRangesIn(n,2)))];
    coms = [coms; newcom(n, commsIn(n))];
end

uiwait(fig);

%% helper functions 
% define GUI behavior  

function addrow(n)
if n > height(btns)
    btns(n,:) = [newpbtn(n), newmbtn(n)];
    etfs(n,:) = [newsef(n,'NaT'), neweef(n,'NaT')];
    coms = [coms; newcom(n, "Enter comment here.")];
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