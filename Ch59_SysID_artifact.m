load("Ch59.mat"); 
ind_bl_str = 4.5e5; ind_bl_end = 9.5e5;

dta = ns2timetable(NS2);
dtaBL = dta(ind_bl_str:ind_bl_end,:);

[predBL, predAll, trnEval, tstEval, A] = fitLTIauton(dtaBL,dta);