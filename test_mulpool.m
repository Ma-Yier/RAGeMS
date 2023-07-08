
% iMM904
% BIOMASS 1521
% OXYGEN 599
% GLUCOSE 508
clear;
%diary MyDiary;
%fprintf('test_ratio_1 ratioGene iMM904, date %s \n',datetime);
%initCobraToolbox;
changeCobraSolver('gurobi');
%tic;
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
max_loop=200;
numMultiStrat=5;
stat=zeros(1226,4);

%ticBytes(pcb)
% glc__D_e 594 
% o2_e 934
i=4;
%parfor (i=1:1226,32)
    if i~=594 && i~=934
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
%end
%toc;
%fprintf('success iMM904 ------------------ \nsuccess rate:%f \n',size(find(stat(:,1)),1)/size(find(stat(:,3)),1));
%diary off;
%filename=sprintf('results/iMM904_raioGene_1_date_%s',datetime('now','TimeZone','Asia/Tokyo','Format','yyyyMMdd'));
%save(filename);

