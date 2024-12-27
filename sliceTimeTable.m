tBaseline = [datetime(2024,10,17,15,58,00), datetime(2024,10,17,15,58,59);...
             datetime(2024,10,17,16,15,00), datetime(2024,10,17,16,22,00);...
             datetime(2024,10,17,16,22,00), datetime(2024,10,17,16,32,00);...
             datetime(2024,10,17,16,45,00), datetime(2024,10,17,16,48,00)]

tStim = [datetime(2024,10,17,16,32,00), datetime(2024,10,17,16,40,00)]

tblBaseline = cell(size(tBaseline)); tblStim = cell(size(tStim));

for it = 1:height(tBaseline)
    t1 = tBaseline(it,1); t2 = tBaseline(it,2);
    tblBaseline{it,1} = NS2tbl((NS2tbl.Time > t1) & (NS2tbl.Time < t2), :);
    tblBaseline{it,2} = NS5tbl((NS5tbl.Time > t1) & (NS5tbl.Time < t2), :);
end

for it = 1:height(tStim)
    t1 = tStim(it,1); t2 = tStim(it,2);
    tblStim{it,1} = NS2tbl((NS2tbl.Time > t1) & (NS2tbl.Time < t2), :);
    tblStim{it,2} = NS5tbl((NS5tbl.Time > t1) & (NS5tbl.Time < t2), :);
end