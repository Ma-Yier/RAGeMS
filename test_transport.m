%------------------------------------------------------------------------ecoli_core-------------------------------------------------------------------------------------

% e_coli_core model
%{
tic;
clear;
load('e_coli_core.mat');
id_biomass=25;
id_carbon=52;
id_oxygen=60;
model=e_coli_core;
model.lb(id_carbon)=-10;
model.lb(id_oxygen)=-10;
model.ub(id_carbon)=-10;
model.ub(id_oxygen)=-10;
stat=zeros(72,3);
max_loop=1000;

for i=2:72
    if i~=25
        [new_model,id_target,TMPR] = introExchange(model,id_biomass,[id_carbon,id_oxygen],i);
        [x_target,alpha] = ratioMethod(new_model,id_biomass,id_target,max_loop,TMPR/max_loop);
        stat(i,1)=x_target;
        stat(i,2)=alpha;
        stat(i,3)=TMPR;
    else
        continue;
    end
end
toc;
%}

% e_coli_core model rationGapMethod
%{
tic;
clear;
load('e_coli_core.mat');
id_biomass=25;
id_carbon=52;
id_oxygen=60;
model=e_coli_core;
model.lb(id_carbon)=-10;
model.lb(id_oxygen)=-10;
model.ub(id_carbon)=-10;
model.ub(id_oxygen)=-10;
stat=zeros(72,3);
max_loop=1000;

for i=2:72
    if i~=25
        [new_model,id_target,TMPR] = introExchange(model,id_biomass,[id_carbon,id_oxygen],i);
        [x_target,alpha] = ratioGapMethod(new_model,id_biomass,id_target,max_loop,TMPR/max_loop);
        stat(i,1)=x_target;
        stat(i,2)=alpha;
        stat(i,3)=TMPR;
    else
        continue;
    end
end
toc;
%}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              

%------------------------------------------------------------------------ecoli_core-------------------------------------------------------------------------------------

%------------------------------------------------------------------------iJO1366----------------------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
%{
% iJO1366 model
clear;
load('iJO1366.mat');
id_biomass=19;
id_carbon=12;
id_oxygen=185;
model=iJO1366;
model.lb(id_carbon)=-10;
model.lb(id_oxygen)=-100;
model.ub(id_carbon)=0;
model.ub(id_oxygen)=0;
stat1=zeros(1805,3);
max_loop=1000;

% 1300 glc__D_e; 1353 o2_e
tic;
parfor i=1:1805
    if i~=1300 && i~=1353
        [new_model,id_target,TMPR] = introExchange(model,id_biomass,[id_carbon,id_oxygen],i);
        [x_target,alpha] = ratioMethod(new_model,id_biomass,id_target,max_loop,TMPR/max_loop);
        stat1(i,:)=[x_target,alpha,TMPR];
    else
        continue;
    end
end
toc;
%}



%{
% iJO1366 model ratioGapMethod
tic;
clear;
load('iJO1366.mat');
id_biomass=19;
id_carbon=12;
id_oxygen=185;
model=iJO1366;
model.lb(id_carbon)=-10;
model.lb(id_oxygen)=-100;
model.ub(id_carbon)=-10;
model.ub(id_oxygen)=-100;
stat=zeros(1805,3);
max_loop=1000;

% 1300 glc__D_e; 1353 o2_e
for i=1:20
    if i~=1300 && i~=1353
        [new_model,id_target,TMPR] = introExchange(model,id_biomass,[id_carbon,id_oxygen],i);
        [x_target,alpha] = ratioGapMethod(new_model,id_biomass,id_target,max_loop,TMPR/max_loop);
        stat(i,1)=x_target;
        stat(i,2)=alpha;
        stat(i,3)=TMPR;
    else
        continue;
    end
end
toc;
%}




%%{
% iJO1366 model ratioMethod
clear;
diary MyDiary;
fprintf('test_transport_1 ratioMethod iJO1366, date %s \n',datetime);
initCobraToolbox;
changeCobraSolver('gurobi');
tic;
load('data/iJO1366.mat');
id_biomass=19;
id_carbon=12;
id_oxygen=185;
model=iJO1366;
model.lb(id_carbon)=-10;
model.lb(id_oxygen)=-100;
model.ub(id_carbon)=0;
model.ub(id_oxygen)=0;
model.lb(id_biomass)=0.05;
stat=zeros(1805,4);
max_loop=1000;

% 1300 glc__D_e; 1353 o2_e
parfor (i=1:1805,32)
    if i~=878 && i~=1369
        tStart=tic;
        [new_model,id_target,TMPR] = introExchange(model,id_biomass,[id_carbon,id_oxygen],i);
        if TMPR>0
            [x_target,alpha] = ratioMethod(new_model,id_biomass,id_target,max_loop,20*TMPR/max_loop);
            stat(i,:)=[x_target,alpha,TMPR,toc(tStart)];
        else
            x_target=0;
            alpha=0;
            stat(i,:)=[x_target,alpha,TMPR,toc(tStart)];
        end
    else
        continue;
    end
end
toc;
diary off;
filename=sprintf('results/iJO1366_raioMethod_1_date_%s',datetime('now','TimeZone','Asia/Tokyo','Format','yyyyMMdd'));
save(filename);
%}

%------------------------------------------------------------------------iJO1366----------------------------------------------------------------------------------------

%------------------------------------------------------------------------iJR904-----------------------------------------------------------------------------------------


%{
% iJR904 model

clear;
load('iJR904.mat');
id_biomass=269;
id_carbon=398;
id_oxygen=500;
model=iJR904;
model.lb(id_carbon)=-10;
model.lb(id_oxygen)=-20;
model.ub(id_carbon)=0;
model.ub(id_oxygen)=0;
model.lb(id_biomass)=0.05;
stat=zeros(1075,3);
max_loop=1000;

% 436 glc__D_e; 579 o2_e
tic;
parfor (i=1:761,6)
    if i~=436 && i~=579
        [new_model,id_target,TMPR] = introExchange(model,id_biomass,[id_carbon,id_oxygen],i);
        [x_target,alpha] = ratioMethod(new_model,id_biomass,id_target,max_loop,TMPR/max_loop);
        stat(i,:)=[x_target,alpha,TMPR];
    else
        continue;
    end
end
toc;
%}

%------------------------------------------------------------------------iJR904-----------------------------------------------------------------------------------------

%------------------------------------------------------------------------iAF1260----------------------------------------------------------------------------------------

% iAF1260 model
%rxns: biomass 926; glc 975; o2 864%
%%{
clear;
diary MyDiary;
fprintf('test_transport_1 ratioMethod iAF1260, date %s \n',datetime);
initCobraToolbox;
changeCobraSolver('gurobi');
tic;
load('data/iAF1260.mat');
id_biomass=926;
id_carbon=975;
id_oxygen=864;
model=iAF1260;
model.lb(id_carbon)=-10;
model.lb(id_oxygen)=-5;
model.ub(id_carbon)=0;
model.ub(id_oxygen)=0;
opt1=optimizeCbModel(model);
TMGR=opt1.f;
model.lb(id_biomass)=0.05;
stat=zeros(1668,4);
max_loop=1000;

%glc__D_e 878 o2_e 1369
parfor (i=1:1668,32)
    if i~=878 && i~=1369
        tStart=tic;
        [new_model,id_target,TMPR] = introExchange(model,id_biomass,[id_carbon,id_oxygen],i);
        if TMPR>0
            [x_target,alpha] = ratioMethod(new_model,id_biomass,id_target,max_loop,20*TMPR/max_loop);
            stat(i,:)=[x_target,alpha,TMPR,toc(tStart)];
        else
            x_target=0;
            alpha=0;
            stat(i,:)=[x_target,alpha,TMPR,toc(tStart)];
        end
    else
        continue;
    end
end
toc;
diary off;
filename=sprintf('results/iAF1260_raioMethod_1_date_%s',datetime('now','TimeZone','Asia/Tokyo','Format','yyyyMMdd'));
save(filename);
%}


%------------------------------------------------------------------------iAF1260----------------------------------------------------------------------------------------


%------------------------------------------------------------------------iMM904----------------------------------------------------------------------------------------

% iMM904
% BIOMASS 1521
% OXYGEN 599
% GLUCOSE 508
%%{
clear;
diary MyDiary;
fprintf('test_transport_1 ratioMethod iMM904, date %s \n',datetime);
initCobraToolbox;
changeCobraSolver('gurobi');
tic;
load('data/iMM904.mat');
id_biomass=1521;
id_carbon=508;
id_oxygen=599;
model=iMM904;
model.lb(id_carbon)=-15;
model.lb(id_oxygen)=-20;
model.ub(id_carbon)=0;
model.ub(id_oxygen)=0;
opt1=optimizeCbModel(model);
TMGR=opt1.f;
model.lb(id_biomass)=0.05;
stat=zeros(1226,4);
max_loop=1000;

%glc__D_e 594 o2_e 934
parfor (i=1:1226,32)
    if i~=594 && i~=934
        tStart=tic;
        [new_model,id_target,TMPR] = introExchange(model,id_biomass,[id_carbon,id_oxygen],i);
        if TMPR>0
            [x_target,alpha] = ratioMethod(new_model,id_biomass,id_target,max_loop,20*TMPR/max_loop);
            stat(i,:)=[x_target,alpha,TMPR,toc(tStart)];
        else
            x_target=0;
            alpha=0;
            stat(i,:)=[x_target,alpha,TMPR,toc(tStart)];
        end
    else
        continue;
    end
end
toc;
diary off;
filename=sprintf('results/iMM904_raioMethod_1_date_%s',datetime('now','TimeZone','Asia/Tokyo','Format','yyyyMMdd'));
save(filename);
%}


%------------------------------------------------------------------------iMM904----------------------------------------------------------------------------------------




%------------------------------------------------------------------------universe---------------------------------------------------------------------------------------


%{
% universeModel model ratioMethod
clear;
diary MyDiary;
sprintf('test6 ratio universeModel, date %s',datetime);
initCobraToolbox;
changeCobraSolver('ibm_cplex');
tic;
load('universeModel.mat');
id_biomass=269;
id_carbon=398;
id_oxygen=500;
model=universeModel;
model.lb(id_carbon)=-10;
model.lb(id_oxygen)=-20;
model.ub(id_carbon)=0;
model.ub(id_oxygen)=0;
model.lb(id_biomass)=0.05;
stat=zeros(1245,4);
max_loop=1000;

% 436 glc__D_e; 579 o2_e
parfor (i=1:1245,6)
    if i~=436 && i~=579
        tStart=tic;
        [new_model,id_target,TMPR] = introExchange(model,id_biomass,[id_carbon,id_oxygen],i);
        if TMPR>0
            [x_target,alpha] = ratioMethod(new_model,id_biomass,id_target,max_loop,20*TMPR/max_loop);
            stat(i,:)=[x_target,alpha,TMPR,toc(tStart)];
        else
            x_target=0;
            alpha=0;
            stat(i,:)=[x_target,alpha,TMPR,toc(tStart)];
        end
    else
        continue;
    end
end
toc;
diary off;
%}


%------------------------------------------------------------------------universe---------------------------------------------------------------------------------------




















