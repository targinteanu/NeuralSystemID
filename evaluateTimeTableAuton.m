function [x,y,xLin,yLin] = evaluateTimeTableAuton(nss,xtbl,shutoff)

if nargin < 3
    shutoff = false(1,height(xtbl));
end

X = table2array(xtbl)';
x2pred = zeros(width(xtbl),1); % column vector 
Xpred = nan(size(X));
Ypred = nan(length(nss.OutputName),height(xtbl));
XpredLin = Xpred; YpredLin = Ypred;

Th = xtbl.Properties.TimeStep; 
tPred = xtbl.Time + Th; % prediction is one step ahead
Th = seconds(Th);

for t = 1:height(xtbl)
    if ~shutoff(t)
        x = X(:,t); xLin = x;
        %{
        syslin = linearize(nss,x);
        syslin = c2d(syslin,Th); 
        %}
    else
        x = x2pred; xLin = x2predLin;
    end
    [x2pred,ypred] = evaluate(nss,x);
    %x2predLin = syslin.A*xLin; ypredLin = syslin.C*x2predLin;
    Xpred(:,t) = x2pred; Ypred(:,t) = ypred;
    %XpredLin(:,t) = x2predLin; YpredLin(:,t) = ypredLin;
end

x = array2timetable(Xpred',"RowTimes",tPred,"VariableNames",nss.StateName);
y = array2timetable(Ypred',"RowTimes",tPred,"VariableNames",nss.OutputName);
xLin = array2timetable(XpredLin',"RowTimes",tPred,"VariableNames",nss.StateName);
yLin = array2timetable(YpredLin',"RowTimes",tPred,"VariableNames",nss.OutputName);

end