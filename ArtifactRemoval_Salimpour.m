% ----------------------------------------------------------------------- %
% ----                                                                --- %
% ----                  Stimulus Artifact Removal                     --- %
% ----                    Bata Phase Detection                        --- %
% ----                                                                --- %
% ----------------------------------------------------------------------- %

% ----                        Triger times                           ---- %
 
Tr=NS2.Data(64,:);

b=diff(Tr);
bb=find(b > 20000);
bbb=diff(bb);
er=find(bbb<50)+1;
Triger=setdiff(bb, bb(er));
ss=zeros(1,length(Tr));
ss(Triger)=28000;

plot(Tr);
hold on
plot(ss);


% ----                       End of Triger m                         ---- %


% ----------------------------------------------------------------------- %
% ----                    Stimulus Artifact Removal                  ---- %

Ch59=double(NS2.Data(59,:));
figure;
plot(Ch59);
hold on
plot(ss);


Ch59_SAR=Ch59;

for i=1:length(Triger)

Sw=Ch59_SAR(Triger(i)-1000+1:Triger(i));
xar = iddata(Sw',[],0.001);
mb = ar(xar,50,'approach','yw');
p = forecast(mb,xar(end-50+1:end),15);
yp=p.OutputData';
Ch59_SAR(Triger(i)+1:Triger(i)+15)=yp;

end

figure;
plot(Ch59);
hold on
plot(Ch59_SAR);

% ----              End of Stimulus Artifact Removal                 ---- %
% ----------------------------------------------------------------------- %
% ----                    Phase @ Triger Times                       ---- %






