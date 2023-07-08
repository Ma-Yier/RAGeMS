%---------------------------------------------------------------------------iML1515--------------------------------------------------------------------------------------
%%{
% iML1515
% BIOMASS 2669
% OXYGEN 1982
% GLUCOSE 181
%clear;
%diary MyDiary;
%fprintf('test_ratio_4 ratioGene iML1515, date %s \n',datetime);
%initCobraToolbox;
changeCobraSolver('ibm_cplex');
%tic;
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

% glc__D_e 1203 o2_e 1250
%for i=1:1877
for i=101:500
    if stat(i,1)>0
        [new_model,id_target,TMPR] = introExchange(model,id_biomass,[id_carbon,id_oxygen],i);
        if TMPR>0
            [preNumKnockouts,afterNumKnockouts,newKnockouts]=cutKnockouts(new_model,kogene{i,1},id_target,id_biomass);
            disp('========');
            disp(preNumKnockouts);
            disp(i);
            disp(afterNumKnockouts);
            disp('========');
        end
    end
end

%toc;
%fprintf('success iML1515 ------------------ \nsuccess rate:%f \n',size(find(stat(:,1)>0.001),1)*100/size(find(stat(:,3)>0.001),1));
%diary off;
%filename=sprintf('results/iML1515_raioGene_5_date_%s',datetime('now','TimeZone','Asia/Tokyo','Format','yyyyMMdd'));
%save(filename);
%}


%---------------------------------------------------------------------------iML1515--------------------------------------------------------------------------------------

