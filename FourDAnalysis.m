%% Comment Block
%{
Analysis Code: Written by Margaret Elise Bullock, last updated 4/29/2022

Will read .mat files generated from runParameterScan, and analyze/visualize

Dependencies: 
Hartigan Dip Test
HeatMaps
Violin
%}

%% Heatmaps of 4D parameter space

%load the parameter space used for runParameterScan
load('important.mat')

Time = ismember(savedItems,'Time');
Average = zeros(length(parameterz),length(parametery),length(parameterm),length(parameterx),length(savedItems)-4);
Variance = Average;
Median = Average;
Dip = Average;
ptest = Average;
bimodal=Average;
xlow=Average;
xup = Average;
m=1;

for u=1:length(parameterm) %PPRR
    for p=1:length(parameterz) %BTR
        for i=1:length(parametery) %BIR
            for j=1:length(parameterx) %PBR
                load(['run_', paramNames{1},'=',strrep(num2str(parameterx(j)),'.','_'),...
                    '_',paramNames{2},'=',strrep(num2str(parametery(i)),'.','_'),...
                    '_',paramNames{3},'=',strrep(num2str(parameterz(p)),'.','_'),...
                    '_',paramNames{4},'=',strrep(num2str(parameterm(u)),'.','_')]);
                
                index = find(allValues(1,Time,:)==endtime);
                
                for k=1:length(savedItems)-1
                    Average(p,i,u,j,k) = mean(allValues(:,k+1,index));
                    Variance(p,i,u,j,k) = var(allValues(:,k+1,index));
                    Median(p,i,u,j,k) = median(allValues(:,k+1,index));
                    [Dip(p,i,u,j,k), ptest(p,i,u,j,k),xlow(p,i,u,j,k),xup(p,i,u,j,k)] = hartigansdipsigniftest(allValues(:,k+1,index),10000);
                    bimodal(p,i,u,j,k) = double(Dip(p,i,u,j,k) > 0.05 & ptest(p,i,u,j,k) < 0.15);
                    
                    if bimodal(p,i,u,j,k)==0
                        if xlow(p,i,u,j,k)~=0 && xup(p,i,u,j,k)~=0
                            bimodal(p,i,u,j,k)=2;
                        end
                    end
                end
            end
        end
    end
end

fano = Variance./Average;
cv = sqrt(Variance)./Average;
skew = (Average-Median)./sqrt(Variance);
     
save('stats.mat','Average','Variance','Median','fano','cv','skew',...
        'Dip','ptest','bimodal');
 
    
 % plot the heatmaps as tiled maps

createTiledHeatmap(Average,'Average',paramNames,savedItems,parameterx, ...
    parametery, parameterz, parameterm,k,0,endtime)

createTiledHeatmap(Variance,'Variance',paramNames,savedItems,parameterx, ...
    parametery, parameterz, parameterm,k,0,endtime)

createTiledHeatmap(Median,'Median',paramNames,savedItems,parameterx, ...
    parametery, parameterz, parameterm,k,0,endtime)

createTiledHeatmap(Dip,'Dip',paramNames,savedItems,parameterx, ...
    parametery, parameterz, parameterm,k,0,endtime)

createTiledHeatmap(ptest,'Ptest',paramNames,savedItems,parameterx, ...
    parametery, parameterz, parameterm,k,0,endtime)

createTiledHeatmap(fano,'Fano',paramNames,savedItems,parameterx, ...
    parametery, parameterz, parameterm,k,0,endtime)

createTiledHeatmap(cv,'CV',paramNames,savedItems,parameterx, ...
    parametery, parameterz, parameterm,k,0,endtime)

createTiledHeatmap(skew,'Skew',paramNames,savedItems,parameterx, ...
    parametery, parameterz, parameterm,k,0,endtime)

createTiledHeatmap(bimodal,'Bimodal Logical',paramNames,savedItems,parameterx, ...
    parametery, parameterz, parameterm,k,0,endtime)

%%