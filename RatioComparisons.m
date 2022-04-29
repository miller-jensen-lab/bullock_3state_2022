
clear all
close all

load('run_PBR=10_BIR=0_1_BTR=0_05_PPRR=10.mat')
load('important.mat')
figure('Name','All Trajectories');
for p=1:size(allValues,1)
    x(:)=allValues(p,1,:);
    y(:)=allValues(p,6,:);
    hold on;
    plot(x,y)
    
    xlim([0 endtime])
    ylim([0 50])
end
set(gca,'xtick',[],'ytick',[])
box on;

figure('Name','Pie Chart');
index = allValues(1,1,:)==endtime;
starting_values = sum(allValues(:,2:4,1));
ending_values = sum(allValues(:,2:4,index));

x1 = starting_values';
labels = {'UP','AP','BP'};

for h = 1:3
    if x1(h) == 0
        labels{h} = ' ';
    end
end

x1(x1 == 0)=realmin('double');
x1 = x1./sum(x1);
p=pie(x1,labels);



figure('Name','Trajectories of one');
for p=1:size(allValues,1)
    if allValues(p,2,1)==1
        subplot(1,3,1)
        x(:)=allValues(p,1,:);
        y(:)=allValues(p,6,:);
        hold on;
        plot(x,y)
        
        xlim([0 endtime])
        ylim([0 50])
    elseif allValues(p,3,1)==1
        subplot(1,3,2)
        x(:)=allValues(p,1,:);
        y(:)=allValues(p,6,:);
        hold on;
        plot(x,y)
        
        xlim([0 endtime])
        ylim([0 50])
    else
        subplot(1,3,3)
        x(:)=allValues(p,1,:);
        y(:)=allValues(p,6,:);
        hold on;
        plot(x,y)
        xlim([0 endtime])
        ylim([0 50])
    end
end
set(gca,'xtick',[],'ytick',[])
box on;

%%
Time = 1;
UP = 2;
AP = 3;
BP = 4;
RNA = 5;
Protein =6;
nboot=10000; alpha = 0.05;
numOfRuns = 100;
endtime = [0 2 4];
index(1,:) = allValues(1,Time,:)==endtime(1);
index(2,:) = allValues(1,Time,:)==endtime(2);
index(3,:) = allValues(1,Time,:)==endtime(3);
%go through time points
for k=2:4
    for j=1:length(endtime)
        averageRNA(k,j) = mean(allValues((allValues(:,k,1)==1),RNA,index(j,:)));
        ci = bootci(nboot,{@mean,allValues((allValues(:,k,1)==1),RNA,index(j,:)) },'alpha',alpha);
        erravRNA(k,j,1) = ci(1)-averageRNA(k,j);
        erravRNA(k,j,2) = ci(2)-averageRNA(k,j);
        
        cvRNA(k,j) = cv(allValues((allValues(:,k,1)==1),RNA,index(j,:)));
        ci = bootci(nboot,{@cv,allValues((allValues(:,k,1)==1),RNA,index(j,:)) },'alpha',alpha);
        errcvRNA(k,j,1) = ci(1)-cvRNA(k,j);
        errcvRNA(k,j,2) = ci(2)-cvRNA(k,j);
        
        fanoRNA(k,j) = fano(allValues((allValues(:,k,1)==1),RNA,index(j,:)));
        ci = bootci(nboot,{@fano,allValues((allValues(:,k,1)==1),RNA,index(j,:)) },'alpha',alpha);
        errfanoRNA(k,j,1) = ci(1)-fanoRNA(k,j);
        errfanoRNA(k,j,2) = ci(2)-fanoRNA(k,j);
    end
end

%%
%Mean
outputs1=averageRNA(2:4,:);


figure('Name',['Average mRNA']);
b=bar([1:3],outputs1);
hu(:,:,:)=[erravRNA(2:4,:,:)];
hold on;

upbound=[hu(:,:,1)]';
downbound=[hu(:,:,2)]';

% Calculate the number of bars in each group
nbars = size(outputs1, 2);
% Get the x coordinate of the bars
x = [];
for i = 1:nbars
    x = [x ; b(i).XEndPoints];
end
% Plot the errorbars
errorbar(x',outputs1,upbound,downbound,'k','linestyle','none')
set(gca, 'XTickLabel', ['UP'; 'AP'; 'BP']);
xtickangle(45)
hold off

legend('Basal','2 Hr','4 Hr','Location','Northwest')
title(['Average mRNA'])

%%
%CV
outputs1=cvRNA(2:4,:);


figure('Name',['mRNA CV']);
b=bar([1:3],outputs1);
hu(:,:,:)=[errcvRNA(2:4,:,:)];
hold on;

upbound=[hu(:,:,1)]';
downbound=[hu(:,:,2)]';

% Calculate the number of bars in each group
nbars = size(outputs1, 2);
% Get the x coordinate of the bars
x = [];
for i = 1:nbars
    x = [x ; b(i).XEndPoints];
end
% Plot the errorbars
errorbar(x',outputs1,upbound,downbound,'k','linestyle','none')
set(gca, 'XTickLabel', ['UP'; 'AP'; 'BP']);
xtickangle(45)
hold off

legend('Basal','2 Hr','4 Hr','Location','Northeast')
title(['mRNA CV'])

%%
%FANO
outputs1=fanoRNA(2:4,:);


figure('Name',['mRNA Fano']);
b=bar([1:3],outputs1);
hu(:,:,:)=[errfanoRNA(2:4,:,:)];
hold on;

upbound=[hu(:,:,1)]';
downbound=[hu(:,:,2)]';

% Calculate the number of bars in each group
nbars = size(outputs1, 2);
% Get the x coordinate of the bars
x = [];
for i = 1:nbars
    x = [x ; b(i).XEndPoints];
end
% Plot the errorbars
errorbar(x',outputs1,upbound,downbound,'k','linestyle','none')
set(gca, 'XTickLabel', ['UP'; 'AP'; 'BP']);
xtickangle(45)
hold off

legend('Basal','2 Hr','4 Hr','Location','Northeast')
title(['mRNA Fano'])

%%

parameterx = [0.01 0.05 0.1 0.5 1];
parametery = [0.01 0.05 0.1 0.5 1];
paramNames = {'BIR','BTR'}; %x, y

for i=1:length(parametery) %y
    for j=1:length(parameterx) %x
        
        load(['run_PBR=10',...
                        '_BIR=',strrep(num2str(parameterx(j)),'.','_'),...
                        '_BTR=',strrep(num2str(parametery(i)),'.','_'),...
                        '_PPRR=10']);
        index = allValues(1,1,:)==endtime;
        starting_values = sum(allValues(:,2:4,1));
        ending_values = sum(allValues(:,2:4,index));
        figure(7);
        ax1 = axes('Position',[(j)/(length(parameterx)+2) ((length(parametery)+1)-i)/(length(parametery)+2) 1/(length(parameterx)+2) 1/(length(parametery)+2)]); %[left bottom width height]

        x1 = starting_values';
        labels = {' ',' ',' '};
        savemean(i,j)=mean(allValues(:,6,25));
        for h = 1:3
            if x1(h) == 0
                labels{h} = ' '; 
            end
        end

        x1(x1 == 0)=realmin('double');
        x1 = x1./sum(x1);
        p=pie(ax1,x1,labels);
        
        fanoblock(i,j)=fano(allValues(:,5,1));
        averageblock(i,j)=mean(allValues(:,5,1));
        
        averageprotein24block(i,j)=mean(allValues(:,6,25));
        
                figure(8);
        ax2 = axes('Position',[(j)/(length(parameterx)+2) ((length(parametery)+1)-i)/(length(parametery)+2) 1/(length(parameterx)+3) 1/(length(parametery)+3)]); %[left bottom width height]
        data=allValues(:,6,end)+1;
        [y1,x1] = ksdensity(data,'Support','positive');
        %semilogx(x,y)
        plot(ax2,x1,y1)
        xlim([0 400])
        ylim([0 0.05])
    end
end

figure(7)
ax0 = axes('Position',[0 0 1 1], 'Visible', 'off');

for i=1:length(parameterx)
    text((.5+i)/(length(parameterx)+2), 0.085, num2str(parameterx(i)), 'HorizontalAlignment','right')
end
 H=findobj(gca,'Type','text');
 set(H,'Rotation',45); % tilt
 text(5.5/11, .02,paramNames{1}, 'HorizontalAlignment','center');
text(0.5,10.5/11,'Initialization Ratio', 'HorizontalAlignment','center','FontSize',12,'FontWeight', 'bold');

for i = 1:length(parametery)
   text(0.09, (.5+i)/(length(parametery)+2), num2str(parametery(length(parametery)+1-i)), 'HorizontalAlignment','right')
end

h=text(0.02, 5.5/11,paramNames{2}, 'HorizontalAlignment','center');
set(h,'Rotation',90);

figure('Name',['Fano RNA Basal']);
 heatmaps(fanoblock, parameterx,parametery,'%0.0f','TickAngle',45,'Colorbar',true)
    xlabel(paramNames{1});ylabel(paramNames{2})
    title('Fano RNA Basal')
    
figure('Name',['Mean RNA Basal']);
 heatmaps(averageblock, parameterx,parametery,'%0.0f','TickAngle',45,'Colorbar',true)
    xlabel(paramNames{1});ylabel(paramNames{2})
    title('Mean RNA Basal')
    
figure('Name',['Average Protein 24hr']);
 heatmaps(averageprotein24block, parameterx,parametery,'%0.0f','TickAngle',45,'Colorbar',true)
    xlabel(paramNames{1});ylabel(paramNames{2})
    title('Average Protein 24hr')
    
%%
clear all;
load('important.mat')

load('run_PBR=10_BIR=0_1_BTR=1_PPRR=10.mat')
%lowfano=allValues(1:1000,:,:);
lowfano=allValues((allValues(:,5,1)<=3),:,:);
%lowfano=lowfano(1:2000,:,:);
load('run_PBR=10_BIR=0_05_BTR=0_5_PPRR=10.mat')%313
%medfano=allValues(1:1000,:,:);
medfano=allValues((allValues(:,5,1)<=3),:,:);
%medfano=medfano(1:2000,:,:);
load('run_PBR=10_BIR=0_01_BTR=0_1_PPRR=10.mat')%313
%highfano=allValues(1:1000,:,:);
highfano=allValues((allValues(:,5,1)<=3),:,:);
%highfano=highfano(1:2000,:,:);
RNA = ismember(savedItems,'RNA');

figure('Name','PPRR 10-fold Activation, States at Basal');
t=tiledlayout(1,3);
nexttile
title(t, 'PPRR 10-fold Activation, States at Basal')
starting_values = sum(medfano(:,2:4,1));
x1 = starting_values';
labels = {'UP','AP','BP'};

for h = 1:3
    if x1(h) == 0
        labels{h} = ' ';
    end
end

x1(x1 == 0)=realmin('double');
x1 = x1./sum(x1);

pie(x1,labels);
title('MedFano')

nexttile
index = allValues(1,1,:)==endtime;
starting_values = sum(highfano(:,2:4,1));
x1 = starting_values';
labels = {'UP','AP','BP'};

for h = 1:3
    if x1(h) == 0
        labels{h} = ' ';
    end
end

x1(x1 == 0)=realmin('double');
x1 = x1./sum(x1);


pie(x1,labels);
title('HighFano')
nexttile;
starting_values = sum(lowfano(:,2:4,1));
x1 = starting_values';
labels = {'UP','AP','BP'};

for h = 1:3
    if x1(h) == 0
        labels{h} = ' ';
    end
end

x1(x1 == 0)=realmin('double');
x1 = x1./sum(x1);

pie(x1,labels);
title('LowFano')


figure('Name','PPRR 10-fold Activation, States at 2 hrs');
t=tiledlayout(1,3);
nexttile
title(t, 'PPRR 10-fold Activation, States at 2 hrs')

starting_values = sum(medfano(:,2:4,allValues(1,1,:)==2));
x1 = starting_values';
labels = {'UP','AP','BP'};

for h = 1:3
    if x1(h) == 0
        labels{h} = ' ';
    end
end

x1(x1 == 0)=realmin('double');
x1 = x1./sum(x1);

pie(x1,labels);
title('MedFano')

nexttile
index = allValues(1,1,:)==endtime;
starting_values = sum(highfano(:,2:4,allValues(1,1,:)==2));
x1 = starting_values';
labels = {'UP','AP','BP'};

for h = 1:3
    if x1(h) == 0
        labels{h} = ' ';
    end
end

x1(x1 == 0)=realmin('double');
x1 = x1./sum(x1);


pie(x1,labels);
title('HighFano')
nexttile;
starting_values = sum(lowfano(:,2:4,allValues(1,1,:)==2));
x1 = starting_values';
labels = {'UP','AP','BP'};

for h = 1:3
    if x1(h) == 0
        labels{h} = ' ';
    end
end

x1(x1 == 0)=realmin('double');
x1 = x1./sum(x1);

pie(x1,labels);
title('LowFano')
%%
figure('Name','Trajectories');
endtime=24;
tiledlayout(1,3)
nexttile
title('MedFano')
for p=1:length(medfano)
    x(:)=medfano(p,1,:);
    y(:)=medfano(p,RNA,:);
    hold on;
    plot(x,y)
    xlim([0 endtime])
    ylim([0 50])
    if p==100
        break;
    end
end

nexttile
title('HighFano')
for p=1:length(highfano)
    x(:)=highfano(p,1,:);
    y(:)=highfano(p,RNA,:);
    hold on;
    plot(x,y)
    xlim([0 endtime])
    ylim([0 50])
    if p==100
        break;
    end
end

nexttile
title('LowFano')
for p=1:length(lowfano)
    x(:)=lowfano(p,1,:);
    y(:)=lowfano(p,RNA,:);
    hold on;
    plot(x,y)
    xlim([0 endtime])
    ylim([0 50])
    if p==100
        break;
    end
end
%%
hour(1).name='hr0';
hour(2).name='hr1';
hour(3).name='hr2';
hour(4).name='hr4';
hour(5).name='hr24';
hour(1).fano(2) = get_mvstd_ci(medfano(:,RNA,allValues(1,1,:)==0));
hour(2).fano(2) = get_mvstd_ci(medfano(:,RNA,allValues(1,1,:)==1));
hour(3).fano(2) = get_mvstd_ci(medfano(:,RNA,allValues(1,1,:)==2));
hour(4).fano(2) = get_mvstd_ci(medfano(:,RNA,allValues(1,1,:)==4));
hour(5).fano(2) = get_mvstd_ci(medfano(:,RNA,allValues(1,1,:)==24));
hour(1).fano(3) = get_mvstd_ci(highfano(:,RNA,allValues(1,1,:)==0));
hour(2).fano(3) = get_mvstd_ci(highfano(:,RNA,allValues(1,1,:)==1));
hour(3).fano(3) = get_mvstd_ci(highfano(:,RNA,allValues(1,1,:)==2));
hour(4).fano(3) = get_mvstd_ci(highfano(:,RNA,allValues(1,1,:)==4));
hour(5).fano(3) = get_mvstd_ci(highfano(:,RNA,allValues(1,1,:)==24));
hour(1).fano(1) = get_mvstd_ci(lowfano(:,RNA,allValues(1,1,:)==0));
hour(2).fano(1) = get_mvstd_ci(lowfano(:,RNA,allValues(1,1,:)==1));
hour(3).fano(1) = get_mvstd_ci(lowfano(:,RNA,allValues(1,1,:)==2));
hour(4).fano(1) = get_mvstd_ci(lowfano(:,RNA,allValues(1,1,:)==4));
hour(5).fano(1) = get_mvstd_ci(lowfano(:,RNA,allValues(1,1,:)==24));

%%
%Mean
for i=1:5
    outputs1(:,i)=[hour(i).fano.mean];
    hu(:,i,:)=[hour(i).fano.mci]';
    upbound(:,i)=[-hu(:,i,1)+outputs1(:,i)];
    downbound(:,i)=[hu(:,i,2)-outputs1(:,i)];
end

figure('Name',['Average mRNA']);
b=bar([1:5],outputs1');
hold on;

% Calculate the number of bars in each group
nbars = size(outputs1', 2);
% Get the x coordinate of the bars
x = [];
for i = 1:nbars
    x = [x ; b(i).XEndPoints];
end
% Plot the errorbars
errorbar(x',outputs1',upbound',downbound','k','linestyle','none')
set(gca, 'XTickLabel', {hour.name});
xtickangle(45)
hold off

legend({'LowFano','MedFano','HighFano'})
title(['Average mRNA'])

%fano
for i=1:5
    outputs1(:,i)=[hour(i).fano.fano];
    hu(:,i,:)=[hour(i).fano.fci]';
    upbound(:,i)=[-hu(:,i,1)+outputs1(:,i)];
    downbound(:,i)=[hu(:,i,2)-outputs1(:,i)];
end

figure('Name',['Average Fano']);
b=bar([1:5],outputs1');
hold on;

% Calculate the number of bars in each group
nbars = size(outputs1', 2);
% Get the x coordinate of the bars
x = [];
for i = 1:nbars
    x = [x ; b(i).XEndPoints];
end
% Plot the errorbars
errorbar(x',outputs1',upbound',downbound','k','linestyle','none')
set(gca, 'XTickLabel', {hour.name});
xtickangle(45)
hold off

legend({'LowFano','MedFano','HighFano'})
title(['Average Fano'])

%%
clear all;
load('important.mat')

load('run_PBR=10_BIR=1_BTR=0_01_PPRR=10.mat')
lowUP=allValues((allValues(:,5,1)<=6),:,:);
%lowUP=allValues(1:1000,:,:);
load('run_PBR=10_BIR=0_1_BTR=0_05_PPRR=10.mat')%313
%medUP=allValues(1:1000,:,:);
medUP=allValues((allValues(:,5,1)<=3),:,:);
%medUP=medUP(1:1918,:,:);
load('run_PBR=10_BIR=0_01_BTR=0_1_PPRR=10.mat')%313
%highUP=allValues(1:1000,:,:);
highUP=allValues((allValues(:,5,1)<=3),:,:);
%highUP=highUP(1:2000,:,:);
RNA = ismember(savedItems,'RNA');

figure('Name','PPRR 10-fold Activation, States at Basal')
t=tiledlayout(1,3);
nexttile
title(t, 'PPRR 10-fold Activation, States at Basal')
starting_values = sum(medUP(:,2:4,1));
x1 = starting_values';
labels = {'UP','AP','BP'};

for h = 1:3
    if x1(h) == 0
        labels{h} = ' ';
    end
end

x1(x1 == 0)=realmin('double');
x1 = x1./sum(x1);

pie(x1,labels);
title('MedUP')

nexttile
index = allValues(1,1,:)==endtime;
starting_values = sum(highUP(:,2:4,1));
x1 = starting_values';
labels = {'UP','AP','BP'};

for h = 1:3
    if x1(h) == 0
        labels{h} = ' ';
    end
end

x1(x1 == 0)=realmin('double');
x1 = x1./sum(x1);


pie(x1,labels);
title('HighUP')
nexttile;
starting_values = sum(lowUP(:,2:4,1));
x1 = starting_values';
labels = {'UP','AP','BP'};

for h = 1:3
    if x1(h) == 0
        labels{h} = ' ';
    end
end

x1(x1 == 0)=realmin('double');
x1 = x1./sum(x1);

pie(x1,labels);
title('LowUP')


figure('Name','PPRR 10-fold Activation, States at 2 hrs');
t=tiledlayout(1,3);
nexttile
title(t, 'PPRR 10-fold Activation, States at 2 hrs')

starting_values = sum(medUP(:,2:4,allValues(1,1,:)==2));
x1 = starting_values';
labels = {'UP','AP','BP'};

for h = 1:3
    if x1(h) == 0
        labels{h} = ' ';
    end
end

x1(x1 == 0)=realmin('double');
x1 = x1./sum(x1);

pie(x1,labels);
title('MedUP')

nexttile
index = allValues(1,1,:)==endtime;
starting_values = sum(highUP(:,2:4,allValues(1,1,:)==2));
x1 = starting_values';
labels = {'UP','AP','BP'};

for h = 1:3
    if x1(h) == 0
        labels{h} = ' ';
    end
end

x1(x1 == 0)=realmin('double');
x1 = x1./sum(x1);


pie(x1,labels);
title('HighUP')
nexttile;
starting_values = sum(lowUP(:,2:4,allValues(1,1,:)==2));
x1 = starting_values';
labels = {'UP','AP','BP'};

for h = 1:3
    if x1(h) == 0
        labels{h} = ' ';
    end
end

x1(x1 == 0)=realmin('double');
x1 = x1./sum(x1);

pie(x1,labels);
title('LowUP')
%%
figure;
endtime=24;
tiledlayout(1,3)
nexttile
title('MedUP')
for p=1:size(medUP,1)
    x(:)=medUP(p,1,:);
    y(:)=medUP(p,RNA,:);
    hold on;
    plot(x,y)
    xlim([0 endtime])
    ylim([0 50])
    if p==100
        break;
    end
end

nexttile
title('HighUP')
for p=1:length(highUP)
    x(:)=highUP(p,1,:);
    y(:)=highUP(p,RNA,:);
    hold on;
    plot(x,y)
    xlim([0 endtime])
    ylim([0 50])
    if p==100
        break;
    end
end

nexttile
title('LowUP')
for p=1:size(lowUP,1)
    x(:)=lowUP(p,1,:);
    y(:)=lowUP(p,RNA,:);
    hold on;
    plot(x,y)
    xlim([0 endtime])
    ylim([0 50])
    if p==100
        break;
    end
end
%%
hour(1).name='hr0';
hour(2).name='hr1';
hour(3).name='hr2';
hour(4).name='hr4';
hour(5).name='hr24';
hour(1).fano(2) = get_mvstd_ci(medUP(:,RNA,allValues(1,1,:)==0));
hour(2).fano(2) = get_mvstd_ci(medUP(:,RNA,allValues(1,1,:)==1));
hour(3).fano(2) = get_mvstd_ci(medUP(:,RNA,allValues(1,1,:)==2));
hour(4).fano(2) = get_mvstd_ci(medUP(:,RNA,allValues(1,1,:)==4));
hour(5).fano(2) = get_mvstd_ci(medUP(:,RNA,allValues(1,1,:)==24));
hour(1).fano(3) = get_mvstd_ci(highUP(:,RNA,allValues(1,1,:)==0));
hour(2).fano(3) = get_mvstd_ci(highUP(:,RNA,allValues(1,1,:)==1));
hour(3).fano(3) = get_mvstd_ci(highUP(:,RNA,allValues(1,1,:)==2));
hour(4).fano(3) = get_mvstd_ci(highUP(:,RNA,allValues(1,1,:)==4));
hour(5).fano(3) = get_mvstd_ci(highUP(:,RNA,allValues(1,1,:)==24));
hour(1).fano(1) = get_mvstd_ci(lowUP(:,RNA,allValues(1,1,:)==0));
hour(2).fano(1) = get_mvstd_ci(lowUP(:,RNA,allValues(1,1,:)==1));
hour(3).fano(1) = get_mvstd_ci(lowUP(:,RNA,allValues(1,1,:)==2));
hour(4).fano(1) = get_mvstd_ci(lowUP(:,RNA,allValues(1,1,:)==4));
hour(5).fano(1) = get_mvstd_ci(lowUP(:,RNA,allValues(1,1,:)==24));

%%
%Mean
for i=1:5
    outputs1(:,i)=[hour(i).fano.mean];
    hu(:,i,:)=[hour(i).fano.mci]';
    upbound(:,i)=[-hu(:,i,1)+outputs1(:,i)];
    downbound(:,i)=[hu(:,i,2)-outputs1(:,i)];
end

figure('Name',['Average mRNA']);
b=bar([1:5],outputs1');
hold on;

% Calculate the number of bars in each group
nbars = size(outputs1', 2);
% Get the x coordinate of the bars
x = [];
for i = 1:nbars
    x = [x ; b(i).XEndPoints];
end
% Plot the errorbars
errorbar(x',outputs1',upbound',downbound','k','linestyle','none')
set(gca, 'XTickLabel', {hour.name});
xtickangle(45)
hold off

legend({'LowUP','MedUP','HighUP'})
title(['Average mRNA'])

%fano
for i=1:5
    outputs1(:,i)=[hour(i).fano.fano];
    hu(:,i,:)=[hour(i).fano.fci]';
    upbound(:,i)=[-hu(:,i,1)+outputs1(:,i)];
    downbound(:,i)=[hu(:,i,2)-outputs1(:,i)];
end

figure('Name',['Average Fano']);
b=bar([1:5],outputs1');
hold on;

% Calculate the number of bars in each group
nbars = size(outputs1', 2);
% Get the x coordinate of the bars
x = [];
for i = 1:nbars
    x = [x ; b(i).XEndPoints];
end
% Plot the errorbars
errorbar(x',outputs1',upbound',downbound','k','linestyle','none')
set(gca, 'XTickLabel', {hour.name});
xtickangle(45)
hold off

legend({'LowUP','MedUP','HighUP'})
title(['Average Fano'])
