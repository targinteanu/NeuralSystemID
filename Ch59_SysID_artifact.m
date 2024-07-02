load("Ch59.mat"); 
ind_bl_str = 4.5e5; ind_bl_end = 9.5e5;

Tr=NS2.Data(64,:);
Tr = Tr > 1e4;

dta = ns2timetable(NS2); 
dta = dta(:,1:63); % exclude stim (analog in)
dtaBL = dta(ind_bl_str:ind_bl_end,:);

tic
[predBL, predAll, trnEval, tstEval, A] = fitLTIauton(dtaBL,dta);
toc

tic
[adaptBL,adaptAll,adaptTrnEval,adaptTstEval] = ...
    AID_LTI_auton(dtaBL,dta,[],[],Tr);
toc