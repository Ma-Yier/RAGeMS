%%{
% universeModel model ratioMethod
clear;
diary MyDiary;
sprintf('test1 ratio universal_model, date %s',datetime);
initCobraToolbox;
changeCobraSolver('ibm_cplex');
tic;
load('universal_model.mat');
id_biomass=269;
id_carbon=398;
id_oxygen=500;
model=universal_model;
model.lb(id_carbon)=-10;
model.lb(id_oxygen)=-20;
model.ub(id_carbon)=0;
model.ub(id_oxygen)=0;
model.lb(id_biomass)=0.05;
stat=zeros(1543,4);
max_loop=1000;

% 436 glc__D_e; 579 o2_e
parfor (i=1:1543,6)
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
