


%------------------------------------------------------------------------------iMM904---------------------------------------------------------------------------------------
%%{
% iMM904
% BIOMASS 1521
% OXYGEN 599
% GLUCOSE 508
clear;
diary MyDiary;
fprintf('test_ratio_1 ratioGene iMM904, date %s \n',datetime);
initCobraToolbox;
changeCobraSolver('gurobi');
tic;
load('data/iMM904.mat');
id_biomass=1521;
id_carbon=508;
id_oxygen=599;
model=iMM904;
model.lb(id_biomass)=0.05;
model.lb(id_carbon)=-15;
model.lb(id_oxygen)=-20;
model.ub(id_carbon)=0;
model.ub(id_oxygen)=0;
opt=optimizeCbModel(model);
TMGR=opt.f;
max_loop=1000;
numMultiStrat=10;
stat=zeros(1226,4);
kogene=cell(1226,1);


% glc__D_e 594 
% o2_e 934
parpool=(32);
parfor (i=1:1226,32)
    if i~=594 && i~=934
        tStart=tic;
        [new_model,id_target,TMPR] = introExchange(model,id_biomass,[id_carbon,id_oxygen],i);
        if TMPR>0
            [x_target,alpha,knockout] = ratioGene(new_model,id_biomass,id_target,TMGR,max_loop,20*TMPR/max_loop,numMultiStrat);
            stat(i,:)=[x_target,alpha,TMPR,toc(tStart)];
            kogene{i,1}=knockout;
        else
            x_target=0;
            alpha=0;
            stat(i,:)=[x_target,alpha,TMPR,toc(tStart)];
        end
    end
end
toc;
fprintf('success iMM904 ------------------ \nsuccess:%f \n',numel(find(stat(:,1)>0.001)));
diary off;
filename=sprintf('results/iMM904core_raioGene_1_date_%s',datetime('now','TimeZone','Asia/Tokyo','Format','yyyyMMdd'));
save(filename);
delete(gcp('nocreate'));
%}

%------------------------------------------------------------------------------iMM904---------------------------------------------------------------------------------------







%-----------------------------------------------------------------------------iND750----------------------------------------------------------------------------------------

%{
% iND750
% BIOMASS 1265
% OXYGEN 6
% GLUCOSE 437
clear;
diary MyDiary;
fprintf('test_ratio_1 ratioGene iND750, date %s \n',datetime);
initCobraToolbox;
changeCobraSolver('gurobi');
tic;
load('data/iND750.mat');
id_biomass=1265;
id_carbon=437;
id_oxygen=6;
model=iND750;
model.lb(id_biomass)=0.05;
model.lb(id_carbon)=-10;
model.lb(id_oxygen)=-25;
model.ub(id_carbon)=0;
model.ub(id_oxygen)=0;
opt=optimizeCbModel(model);
TMGR=opt.f;
max_loop=1000;
numMultiStrat=10;
stat=zeros(1059,4);

% glc__D_e 55 o2_e 95
parfor (i=1:1059,32)
    if i~=55 && i~=95
        tStart=tic;
        [new_model,id_target,TMPR] = introExchange(model,id_biomass,[id_carbon,id_oxygen],i);
        if TMPR>0
            [x_target,alpha] = ratioGene(new_model,id_biomass,id_target,TMGR,max_loop,20*TMPR/max_loop,numMultiStrat);
            stat(i,:)=[x_target,alpha,TMPR,toc(tStart)];
        else
            x_target=0;
            alpha=0;
            stat(i,:)=[x_target,alpha,TMPR,toc(tStart)];
        end
    end
end


toc;
diary off;
filename=sprintf('results/iND750_raioGene_1_date_%s',datetime('now','TimeZone','Asia/Tokyo','Format','yyyyMMdd'));
save(filename);
%}
%-----------------------------------------------------------------------------iND750----------------------------------------------------------------------------------------






%-----------------------------------------------------------------------------iJR904----------------------------------------------------------------------------------------

%%{
% iJR904
% BIOMASS 269
% OXYGEN 500
% GLUCOSE 398
clear;
diary MyDiary;
fprintf('test_ratio_1 ratioGene iJR904, date %s \n',datetime);
initCobraToolbox;
changeCobraSolver('gurobi');
tic;
load('data/iJR904');
%load('data/iJR904.mat');
id_biomass=269;
id_carbon=398;
id_oxygen=500;
model=iJR904;
model.lb(id_biomass)=0.05;
model.lb(id_carbon)=-15;
model.lb(id_oxygen)=-15;
model.ub(id_carbon)=0;
model.ub(id_oxygen)=0;
opt=optimizeCbModel(model);
TMGR=opt.f;
max_loop=1000;
numMultiStrat=10;
stat=zeros(761,4);
kogene=cell(761,1);

% glc__D_e 436 o2_e 579
parpool(32);
parfor (i=1:761,32)
    if i~=436 && i~=579
        tStart=tic;
        [new_model,id_target,TMPR] = introExchange(model,id_biomass,[id_carbon,id_oxygen],i);
        if TMPR>0
            [x_target,alpha,knockout] = ratioGene(new_model,id_biomass,id_target,TMGR,max_loop,20*TMPR/max_loop,numMultiStrat);
            stat(i,:)=[x_target,alpha,TMPR,toc(tStart)];
            kogene{i,1}=knockout;
        else
            x_target=0;
            alpha=0;
            stat(i,:)=[x_target,alpha,TMPR,toc(tStart)];
        end
    end
end


toc;
fprintf('success iJR904 ------------------ \nsuccess:%f \n',numel(find(stat(:,1)>0.001)));
diary off;
filename=sprintf('results/iJR904_raioGene_1_date_%s',datetime('now','TimeZone','Asia/Tokyo','Format','yyyyMMdd'));
save(filename);
delete(gcp('nocreate'));
%}
%-----------------------------------------------------------------------------iJR904----------------------------------------------------------------------------------------

%---------------------------------------------------------------------------iML1515--------------------------------------------------------------------------------------
%{
% iML1515
% BIOMASS 2669
% OXYGEN 1982
% GLUCOSE 181
clear;
diary MyDiary;
fprintf('test_ratio_4 ratioGene iML1515, date %s \n',datetime);
initCobraToolbox;
changeCobraSolver('gurobi');
tic;
load('data/iML1515.mat');
id_biomass=2669;
id_oxygen=1982;
id_carbon=181;
model=iML1515;
model.lb(id_biomass)=0.05;
model.lb(id_carbon)=-15;
model.lb(id_oxygen)=-20;
model.ub(id_carbon)=0;
model.ub(id_oxygen)=0;
opt=optimizeCbModel(model);
TMGR=opt.f;
max_loop=1000;
numMultiStrat=10;
stat=zeros(1877,4);
kogene=cell(1877,1);



% glc__D_e 1203 o2_e 1250
parfor (i=1:1877,32)
    if i~=1203 && i~=1250
        tStart=tic;
        [new_model,id_target,TMPR] = introExchange(model,id_biomass,[id_carbon,id_oxygen],i);
        if TMPR>0
            [x_target,alpha,knockout] = ratioGene(new_model,id_biomass,id_target,TMGR,max_loop,20*TMPR/max_loop,numMultiStrat);
            stat(i,:)=[x_target,alpha,TMPR,toc(tStart)];
            kogene{i,1}=knockout;
        else
            x_target=0;
            alpha=0;
            stat(i,:)=[x_target,alpha,TMPR,toc(tStart)];
        end
    end
end

toc;
fprintf('success iML1515 ------------------ \nsuccess rate:%f \n',size(find(stat(:,1)>0.001),1)*100/size(find(stat(:,3)>0.001),1));
diary off;
filename=sprintf('results/iML1515_raioGene_5_date_%s',datetime('now','TimeZone','Asia/Tokyo','Format','yyyyMMdd'));
save(filename);
%}


%---------------------------------------------------------------------------iML1515--------------------------------------------------------------------------------------



%---------------------------------------------------------------------------iAF1260--------------------------------------------------------------------------------------

%{
% iAF1260
% BIOMASS 926
% OXYGEN 864
% GLUCOSE 975
clear;
diary MyDiary;
fprintf('test_ratio_1 ratioGene iAF1260, date %s \n',datetime);
initCobraToolbox;
changeCobraSolver('gurobi');
tic;
load('data/iAF1260.mat');
id_biomass=926;
id_oxygen=864;
id_carbon=975;
model=iAF1260;
model.lb(id_biomass)=0.05;
model.lb(id_carbon)=-15;
model.lb(id_oxygen)=-15;
model.ub(id_carbon)=0;
model.ub(id_oxygen)=0;
opt=optimizeCbModel(model);
TMGR=opt.f;
max_loop=1000;
numMultiStrat=10;
stat=zeros(1668,4);



% glc__D_e 878 o2_e 1369
parfor (i=1:1668,32)
    if i~=878 && i~=1369
        tStart=tic;
        [new_model,id_target,TMPR] = introExchange(model,id_biomass,[id_carbon,id_oxygen],i);
        if TMPR>0
            [x_target,alpha] = ratioGene(new_model,id_biomass,id_target,TMGR,max_loop,20*TMPR/max_loop,numMultiStrat);
            stat(i,:)=[x_target,alpha,TMPR,toc(tStart)];
        else
            x_target=0;
            alpha=0;
            stat(i,:)=[x_target,alpha,TMPR,toc(tStart)];
        end
    end
end

toc;
fprintf('------------------ success iAF1260 ------------------ \nsuccess rate:%f \n',size(find(stat(:,1)>0.001),1)*100/size(find(stat(:,3)>0.001),1));
diary off;
filename=sprintf('results/iAF1260_raioGene_1_date_%s',datetime('now','TimeZone','Asia/Tokyo','Format','yyyyMMdd'));
save(filename);
%}

%---------------------------------------------------------------------------iAF1260--------------------------------------------------------------------------------------


%---------------------------------------------------------------------------iJO1366--------------------------------------------------------------------------------------

%{
% iJO1366
% BIOMASS 19
% OXYGEN 185
% GLUCOSE 12
clear;
diary MyDiary;
fprintf('test_ratio_1 ratioGene iJO1366, date %s \n',datetime);
initCobraToolbox;
changeCobraSolver('gurobi');
tic;
load('data/iJO1366.mat');
id_biomass=19;
id_oxygen=185;
id_carbon=12;
model=iJO1366;
model.lb(id_biomass)=0.05;
model.lb(id_carbon)=-15;
model.lb(id_oxygen)=-15;
model.ub(id_carbon)=0;
model.ub(id_oxygen)=0;
opt=optimizeCbModel(model);
TMGR=opt.f;
max_loop=1000;
numMultiStrat=10;
stat=zeros(1805,4);



% glc__D_e 1300 o2_e 1353
parfor (i=1:1805,32)
    if i~=1300 && i~=1353
        tStart=tic;
        [new_model,id_target,TMPR] = introExchange(model,id_biomass,[id_carbon,id_oxygen],i);
        if TMPR>0
            [x_target,alpha] = ratioGene(new_model,id_biomass,id_target,TMGR,max_loop,20*TMPR/max_loop,numMultiStrat);
            stat(i,:)=[x_target,alpha,TMPR,toc(tStart)];
        else
            x_target=0;
            alpha=0;
            stat(i,:)=[x_target,alpha,TMPR,toc(tStart)];
        end
    end
end

toc;
fprintf('---------------- success iJO1366 ------------------ \nsuccess rate:%f \n',size(find(stat(:,1)>0.001),1)*100/size(find(stat(:,3)>0.001),1));
diary off;
filename=sprintf('results/iJO1366_raioGene_1_date_%s',datetime('now','TimeZone','Asia/Tokyo','Format','yyyyMMdd'));
save(filename);
%}

%---------------------------------------------------------------------------iAF1260--------------------------------------------------------------------------------------


