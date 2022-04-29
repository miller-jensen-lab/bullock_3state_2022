

load('important.mat')
Time = ismember(savedItems,'Time');
RNA = ismember(savedItems,'RNA');
nboot=10000; alpha = 0.05;
exp_time = [0 2 4];

load(['run_PBR=50_BIR=0_1_BTR=0_5_PPRR=1.mat'])
index(1,:) = allValues(1,Time,:)==exp_time(1);
index(2,:) = allValues(1,Time,:)==exp_time(2);
index(3,:) = allValues(1,Time,:)==exp_time(3);

load('expData.mat')
averageDataRNA_chx = zeros(4,3);
    cvDataRNA_chx = averageDataRNA_chx;
    fanoDataRNA_chx= averageDataRNA_chx;
    
    errDataavRNA_chx = zeros(4,3,2);
    errDatacvRNA_chx = errDataavRNA_chx;
    errDatafanoRNA_chx = errDataavRNA_chx;
   
averageDataRNA_fb = zeros(4,3);
    cvDataRNA_fb = averageDataRNA_fb;
    fanoDataRNA_fb = averageDataRNA_fb;
    
    errDataavRNA_fb = zeros(4,3,2);
    errDatacvRNA_fb = errDataavRNA_fb;
    errDatafanoRNA_fb = errDataavRNA_fb;
    
m=1;
for i = 1:4
    for j = 1:3
    averageDataRNA_chx(i,j) = CHX(2,m);
    errDataavRNA_chx(i,j,1) = CHX(3,m);
    errDataavRNA_chx(i,j,2) = CHX(4,m);
    
    cvDataRNA_chx(i,j) = CHX(8,m);
    errDatacvRNA_chx(i,j,1) = CHX(9,m);
    errDatacvRNA_chx(i,j,2) = CHX(10,m);
    
    fanoDataRNA_chx(i,j) = CHX(14,m);
    errDatafanoRNA_chx(i,j,1) = CHX(15,m);
    errDatafanoRNA_chx(i,j,2) = CHX(16,m);
    
    averageDataRNA_fb(i,j) = Feedback(2,m);
    errDataavRNA_fb(i,j,1) = Feedback(3,m);
    errDataavRNA_fb(i,j,2) = Feedback(4,m);
    
    cvDataRNA_fb(i,j) = Feedback(8,m);
    errDatacvRNA_fb(i,j,1) = Feedback(9,m);
    errDatacvRNA_fb(i,j,2) = Feedback(10,m);
    
    fanoDataRNA_fb(i,j) = Feedback(14,m);
    errDatafanoRNA_fb(i,j,1) = Feedback(15,m);
    errDatafanoRNA_fb(i,j,2) = Feedback(16,m);
    
    m=m+1;
    end
end

actual84(1,:)=averageDataRNA_fb(1,:); %average
actual84(2,:)=fanoDataRNA_fb(1,:); %fano
actual84(3,:)=cvDataRNA_fb(1,:); %cv

actual44(1,:)=averageDataRNA_fb(2,:);
actual44(2,:)=fanoDataRNA_fb(2,:);
actual44(3,:)=cvDataRNA_fb(2,:);

actual66(1,:)=averageDataRNA_fb(3,:);
actual66(2,:)=fanoDataRNA_fb(3,:);
actual66(3,:)=cvDataRNA_fb(3,:);

actual104(1,:)=averageDataRNA_fb(4,:);
actual104(2,:)=fanoDataRNA_fb(4,:);
actual104(3,:)=cvDataRNA_fb(4,:);

for u=1:length(parameterm) %PPRR
    for p=1:length(parameterz) %BTR
        for i=1:length(parametery) %BIR
            for j=1:length(parameterx) %PBR

                load(['run_', paramNames{1},'=',strrep(num2str(parameterx(j)),'.','_'),...
                        '_',paramNames{2},'=',strrep(num2str(parametery(i)),'.','_'),...
                        '_',paramNames{3},'=',strrep(num2str(parameterz(p)),'.','_'),...
                        '_',paramNames{4},'=',strrep(num2str(parameterm(u)),'.','_')]);

                allValues1=allValues((allValues(:,5,1)<=3),:,:);
                
                stored(1,1) = mean(allValues1(:,5,index(1,:)));
                stored(2,1) = var(allValues1(:,5,index(1,:)))./stored(1,1);
                stored(3,1) =sqrt(var(allValues1(:,5,index(1,:))))./stored(1,1);

                stored(1,2) = mean(allValues1(:,5,index(2,:)));
                stored(2,2) = var(allValues1(:,5,index(2,:)))./stored(1,2);
                stored(3,2) =sqrt(var(allValues1(:,5,index(2,:))))./stored(1,2);

                stored(1,3) = mean(allValues1(:,5,index(3,:)));
                stored(2,3) = var(allValues1(:,5,index(3,:)))./stored(1,3);
                stored(3,3) =sqrt(var(allValues1(:,5,index(3,:))))./stored(1,3);

                residual84(p,i,u,j)=sum(sqrt(((actual84-stored).^2))./actual84,'all');
                residual44(p,i,u,j)=sum(sqrt(((actual44-stored).^2))./actual44,'all');
                residual66(p,i,u,j)=sum(sqrt(((actual66-stored).^2))./actual66,'all');
                residual106(p,i,u,j)=sum(sqrt(((actual104-stored).^2))./actual104,'all');
                
                ratiopprr(p,i,u,j)=parameterm(u)/parameterx(j);
                ratiobtr(p,i,u,j)=parameterz(p)/parametery(i);
            end
        end
    end
end

amount=50;
[M84,I84]=mink(residual84(:),amount);
[M44,I44]=mink(residual44(:),amount);
[M66,I66]=mink(residual66(:),amount);
[M106,I106]=mink(residual106(:),amount);

%%
m=1;
V = zeros(amount,amount,amount,amount);
for i=1:amount
    for j=1:amount
        for k=1:amount
            for l=1:amount
               V(i,j,k,l)=double(ratiopprr(I84(i))<ratiopprr(I66(k)));
               V(i,j,k,l)=double(ratiopprr(I84(i))<ratiopprr(I106(l)))+V(i,j,k,l);
               V(i,j,k,l)=double(ratiobtr(I84(i))>ratiobtr(I66(k)))+V(i,j,k,l);
               V(i,j,k,l)=double(ratiobtr(I84(i))>ratiobtr(I106(l)))+V(i,j,k,l);
               V(i,j,k,l)=double(ratiopprr(I44(j))<ratiopprr(I66(k)))+V(i,j,k,l);
               V(i,j,k,l)=double(ratiopprr(I44(j))<ratiopprr(I106(l)))+V(i,j,k,l);
               V(i,j,k,l)=double(ratiobtr(I44(j))>ratiobtr(I66(k)))+V(i,j,k,l);
               V(i,j,k,l)=double(ratiobtr(I44(j))>ratiobtr(I106(l)))+V(i,j,k,l);               
               if V(i,j,k,l)==8
                   holdthese(m,:)=[i j k l];
                   m=m+1;
               end
            end
        end
    end
end

[M,I]=sort(sum(holdthese(:,:),2));

indexes = [I84(holdthese(I(:),1)) I44(holdthese(I(:),2)) I66(holdthese(I(:),3)) I106(holdthese(I(:),4))];
for i=1:length(indexes)
    if indexes(i,3)==indexes(i,4)
        indexes(i,:)=[0 0 0 0];
    end
end       
    
chosen = find(indexes(:,1)~=0, 1, 'first');

%%

% for u=1:length(parameterm) %PPRR
%     for p=1:length(parameterz) %BTR
%         for i=1:length(parametery) %BIR
%             for j=1:length(parameterx) %PBR
%                 
%             residual84(p,i,u,j)
            
            %this order because the residual matrixes are made BTR BIR PPRR PBR but
            %the labels are in PBR BIR BTR PPRR
[I_84(3),I_84(2),I_84(4),I_84(1)]=ind2sub(size(residual84),I84(holdthese(I(chosen),1)));
parameters_84 = [parameterx(I_84(1)) parametery(I_84(2)) parameterz(I_84(3)) parameterm(I_84(4))];

[I_44(3),I_44(2),I_44(4),I_44(1)]=ind2sub(size(residual44),I44(holdthese(I(chosen),2)));
parameters_44 = [parameterx(I_44(1)) parametery(I_44(2)) parameterz(I_44(3)) parameterm(I_44(4))];

[I_66(3),I_66(2),I_66(4),I_66(1)]=ind2sub(size(residual66),I66(holdthese(I(chosen),3)));
parameters_66 = [parameterx(I_66(1)) parametery(I_66(2)) parameterz(I_66(3)) parameterm(I_66(4))];

[I_106(3),I_106(2),I_106(4),I_106(1)]=ind2sub(size(residual106),I106(holdthese(I(chosen),4)));
parameters_106 = [parameterx(I_106(1)) parametery(I_106(2)) parameterz(I_106(3)) parameterm(I_106(4))];

save('parameters','parameters_84','parameters_44','parameters_66','parameters_106')

%% Picked clones        

clones{1,1} =  ['run_', paramNames{1},'=',strrep(num2str(parameters_84(1)),'.','_'),...
                        '_',paramNames{2},'=',strrep(num2str(parameters_84(2)),'.','_'),...
                        '_',paramNames{3},'=',strrep(num2str(parameters_84(3)),'.','_'),...
                        '_',paramNames{4},'=',strrep(num2str(parameters_84(4)),'.','_')]; %8.4 - b
clones{2,1} = ['run_', paramNames{1},'=',strrep(num2str(parameters_44(1)),'.','_'),...
                        '_',paramNames{2},'=',strrep(num2str(parameters_44(2)),'.','_'),...
                        '_',paramNames{3},'=',strrep(num2str(parameters_44(3)),'.','_'),...
                        '_',paramNames{4},'=',strrep(num2str(parameters_44(4)),'.','_')]; %4.4 - b
clones{3,1} = ['run_', paramNames{1},'=',strrep(num2str(parameters_66(1)),'.','_'),...
                        '_',paramNames{2},'=',strrep(num2str(parameters_66(2)),'.','_'),...
                        '_',paramNames{3},'=',strrep(num2str(parameters_66(3)),'.','_'),...
                        '_',paramNames{4},'=',strrep(num2str(parameters_66(4)),'.','_')]; %6.6 - b
clones{4,1} = ['run_', paramNames{1},'=',strrep(num2str(parameters_106(1)),'.','_'),...
                        '_',paramNames{2},'=',strrep(num2str(parameters_106(2)),'.','_'),...
                        '_',paramNames{3},'=',strrep(num2str(parameters_106(3)),'.','_'),...
                        '_',paramNames{4},'=',strrep(num2str(parameters_106(4)),'.','_')]; %10.6 - b

averageRNA = zeros(4,3);
cvRNA = averageRNA;
fanoRNA = averageRNA;

erravRNA = zeros(4,3,2);
errcvRNA = erravRNA;
errfanoRNA = erravRNA;

figure;
%go through clones
for k=1:length(clones)
    load([clones{k,1}])
    allValues1=allValues((allValues(:,5,1)<=3),:,:);
    points{k} = allValues1(:,6,allValues(1,1,:)==24);
    basalpoints{k} = allValues1(:,6,allValues(1,1,:)==0);

    %go through time points
    for j=1:length(exp_time)    
        averageRNA(k,j) = mean(allValues1(:,RNA,index(j,:)));
        ci = bootci(nboot,{@mean,allValues1(:,RNA,index(j,:)) },'alpha',alpha);
        erravRNA(k,j,1) = ci(1)-averageRNA(k,j);
        erravRNA(k,j,2) = ci(2)-averageRNA(k,j);

        cvRNA(k,j) = cv(allValues1(:,RNA,index(j,:)));
        ci = bootci(nboot,{@cv,allValues1(:,RNA,index(j,:)) },'alpha',alpha);
        errcvRNA(k,j,1) = ci(1)-cvRNA(k,j);
        errcvRNA(k,j,2) = ci(2)-cvRNA(k,j);

        fanoRNA(k,j) = fano(allValues1(:,RNA,index(j,:)));
        ci = bootci(nboot,{@fano,allValues1(:,RNA,index(j,:)) },'alpha',alpha);
        errfanoRNA(k,j,1) = ci(1)-fanoRNA(k,j);
        errfanoRNA(k,j,2) = ci(2)-fanoRNA(k,j);
    end
    
    data=allValues(:,6,stepsize*hour+1)+1;
    
    [y,x] = ksdensity(data,'Support','positive');
    semilogx(x,y)
    %plot(x,y)
    hold on

    
end


for k=1:length(clones)
    load([clones{k,1}])
    allValues1=allValues((allValues(:,5,1)<=3),:,:);
    allValues2=allValues((allValues(:,6,1)<=50),:,:);

    %go through time points
    for j=1:25    
        averageprotein1(k,j) = mean(allValues2(:,6,j));
        ci = bootci(nboot,{@mean,allValues2(:,6,j) },'alpha',alpha);
        erravprotein1(k,j,1) = ci(1)-averageprotein1(k,j);
        erravprotein1(k,j,2) = ci(2)-averageprotein1(k,j);

        cvprotein1(k,j) = cv(allValues2(:,6,j));
        ci = bootci(nboot,{@cv,allValues2(:,6,j) },'alpha',alpha);
        errcvprotein1(k,j,1) = ci(1)-cvprotein1(k,j);
        errcvprotein1(k,j,2) = ci(2)-cvprotein1(k,j);

        fanoprotein1(k,j) = fano(allValues1(:,6,j));
        ci = bootci(nboot,{@fano,allValues1(:,6,j) },'alpha',alpha);
        errfanoprotein1(k,j,1) = ci(1)-fanoprotein1(k,j);
        errfanoprotein1(k,j,2) = ci(2)-fanoprotein1(k,j);
        
        averagerna1(k,j) = mean(allValues2(:,5,j));
        ci = bootci(nboot,{@mean,allValues2(:,5,j) },'alpha',alpha);
        erravrna1(k,j,1) = ci(1)-averagerna1(k,j);
        erravrna1(k,j,2) = ci(2)-averagerna1(k,j);

        cvrna1(k,j) = cv(allValues1(:,5,j));
        ci = bootci(nboot,{@cv,allValues1(:,5,j) },'alpha',alpha);
        errcvrna1(k,j,1) = ci(1)-cvrna1(k,j);
        errcvrna1(k,j,2) = ci(2)-cvrna1(k,j);

        fanorna1(k,j) = fano(allValues1(:,5,j));
        ci = bootci(nboot,{@fano,allValues1(:,5,j) },'alpha',alpha);
        errfanorna1(k,j,1) = ci(1)-fanorna1(k,j);
        errfanorna1(k,j,2) = ci(2)-fanorna1(k,j);
    end
    
end

    xlabel('Protein/cell')
    ylabel('PDF')
    %xlim([0,inf])
title('Protein at 24 hrs')
legend('8.4','4.4','6.6','10.6')

figure('Name',['TNF Distributions Violin 24 hours']);
violin(points,'xlabel',{'8.4','4.4','6.6','10.6'},'bw',75)
title('Protein Distributions (24 Hours)')
ylabel('# of proteins')
set(gcf,'renderer','Painters')
yl = ylim;

figure('Name',['TNF Distributions Violin Basal ']);
violin(basalpoints,'xlabel',{'8.4','4.4','6.6','10.6'},'bw',75)
title('Protein Distributions (Basal)')
ylabel('# of proteins')
ylim(yl);
set(gcf,'renderer','Painters')

for i=1:4
percentages(i,2)=sum(points{i}>300)*100/length(points{i});
percentages(i,1)=sum(basalpoints{i}>300)*100/length(basalpoints{i});
end

%%
figure; hold on;
colortypes=['g';'c';'r';'b'];
for k=1:length(clones)
        hold on;
        shadedErrorBar([0:12],averageprotein1(k,1:13),erravprotein1(k,1:13,1),'lineprops',colortypes(k))
        
        ylabel('Average Protein')
        xlabel('Time (hours)')
        title('Protein Trajectories over time')

        
end
        legend('8.4','4.4','6.6','10.6')

figure; hold on;

colortypes=['g';'c';'r';'b'];
for k=1:length(clones)
        hold on;
        shadedErrorBar([0:12],averagerna1(k,1:13),erravrna1(k,1:13,1),'lineprops',colortypes(k))
        
        ylabel('Average mRNA')
        xlabel('Time (hours)')
        title('RNA Trajectories over time')

    
end
        legend('8.4','4.4','6.6','10.6')
%%
figure('Name',['mRNA Average']); 
color = {'xr' '*r' 'or';...
    'xc' '*c' 'oc';...
    'xm' '*m' 'om';...
    'xb' '*b' 'ob'};
hold on;

for m=1:4

    for k = 1:3
        errorbar(averageDataRNA_fb(m,k),averageRNA(m,k),erravRNA(m,k,1),erravRNA(m,k,2),errDataavRNA_fb(m,k,1),errDataavRNA_fb(m,k,2),color{m,k})
    end

end

ylabel('Simulation'); xlabel('Experimental')
title(['mRNA Average']);
getlimit = axis;
axis([0 max(getlimit) 0 max(getlimit)])

[P, S] = polyfit(averageDataRNA_fb,averageRNA,1);
hold on;
plot(averageDataRNA_fb,polyval(P,averageDataRNA_fb))
regression = 1 - (S.normr/norm(averageRNA(:) - mean(averageRNA(:))))^2;
%text(max(getlimit)*0.25,max(getlimit)*0.75,['simul = ',num2str(P(1)),'*exp + ',num2str(P(2))])
text(max(getlimit)*0.75,max(getlimit)*0.25,['R^2 = ',num2str(regression)])
legend('Clone 1 Basal','Clone 1 2hr','Clone 1 4hr','Clone 2 Basal','Clone 2 2hr','Clone 2 4hr','Clone 3 Basal','Clone 3 2hr','Clone 3 4hr','Clone 4 Basal','Clone 4 2hr','Clone 4 4hr','Location','best')


figure('Name',['mRNA CV']); 
color = {'xr' '*r' 'or';...
    'xc' '*c' 'oc';...
    'xm' '*m' 'om';...
    'xb' '*b' 'ob'};
hold on;

for m=1:4

    for k = 1:3
        errorbar(cvDataRNA_fb(m,k),cvRNA(m,k),errcvRNA(m,k,1),errcvRNA(m,k,2),errDatacvRNA_fb(m,k,1),errDatacvRNA_fb(m,k,2),color{m,k})
    end

end

ylabel('Simulation'); xlabel('Experimental')
title(['mRNA CV']);
getlimit = axis;
axis([0 max(getlimit) 0 max(getlimit)])

[P, S] = polyfit(cvDataRNA_fb,cvRNA,1);
hold on;
plot(cvDataRNA_fb,polyval(P,cvDataRNA_fb))
regression = 1 - (S.normr/norm(cvRNA(:) - mean(cvRNA(:))))^2;
%text(max(getlimit)*0.25,max(getlimit)*0.75,['simul = ',num2str(P(1)),'*exp + ',num2str(P(2))])
text(max(getlimit)*0.75,max(getlimit)*0.25,['R^2 = ',num2str(regression)])
legend('Clone 1 Basal','Clone 1 2hr','Clone 1 4hr','Clone 2 Basal','Clone 2 2hr','Clone 2 4hr','Clone 3 Basal','Clone 3 2hr','Clone 3 4hr','Clone 4 Basal','Clone 4 2hr','Clone 4 4hr','Location','best')

figure('Name',['mRNA Fano']); 
color = {'xr' '*r' 'or';...
    'xc' '*c' 'oc';...
    'xm' '*m' 'om';...
    'xb' '*b' 'ob'};
hold on;

for m=1:4

    for k = 1:3
        errorbar(fanoDataRNA_fb(m,k),fanoRNA(m,k),errfanoRNA(m,k,1),errfanoRNA(m,k,2),errDatafanoRNA_fb(m,k,1),errDatafanoRNA_fb(m,k,2),color{m,k})
    end

end

ylabel('Simulation'); xlabel('Experimental')
title(['mRNA Fano']);
getlimit = axis;
axis([0 max(getlimit) 0 max(getlimit)])

[P, S] = polyfit(fanoDataRNA_fb,fanoRNA,1);
hold on;
plot(fanoDataRNA_fb,polyval(P,fanoDataRNA_fb))
regression = 1 - (S.normr/norm(fanoRNA(:) - mean(fanoRNA(:))))^2;
%text(max(getlimit)*0.25,max(getlimit)*0.75,['simul = ',num2str(P(1)),'*exp + ',num2str(P(2))])
text(max(getlimit)*0.75,max(getlimit)*0.25,['R^2 = ',num2str(regression)])
legend('Clone 1 Basal','Clone 1 2hr','Clone 1 4hr','Clone 2 Basal','Clone 2 2hr','Clone 2 4hr','Clone 3 Basal','Clone 3 2hr','Clone 3 4hr','Clone 4 Basal','Clone 4 2hr','Clone 4 4hr','Location','best')


%%

load([clones{1,1}])
c84=allValues;
load([clones{2,1}])%313
c44=allValues;
load([clones{3,1}])%313
c66=allValues;
load([clones{4,1}])%313
c106=allValues;
        

figure;
t=tiledlayout(1,4);
nexttile
title(t, 'PPRR 10-fold Activation, States at Basal') 
        

starting_values = sum(c84(:,2:4,1));
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
title('8.4')    





nexttile;
starting_values = sum(c44(:,2:4,1));
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
title('4.4')



nexttile;
starting_values = sum(c66(:,2:4,1));
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
title('6.6')
        
        
nexttile;
starting_values = sum(c106(:,2:4,1));
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
title('10.6')
        



figure;
t=tiledlayout(1,4);
nexttile
title(t, 'PPRR 10-fold Activation, States at 2 hrs') 

starting_values = sum(c84(:,2:4,allValues(1,1,:)==2));
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
title('8.4')
        


nexttile
index = allValues(1,1,:)==exp_time;
starting_values = sum(c44(:,2:4,allValues(1,1,:)==2));
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
title('4.4')



nexttile;
starting_values = sum(c66(:,2:4,allValues(1,1,:)==2));
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
title('6.6')

nexttile;
starting_values = sum(c106(:,2:4,allValues(1,1,:)==2));
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
title('10.6')

divisionmean = @(a,b,c) (mean(a+b)/mean(c)); 
ach3(1) = mean(c106(:,3,1)+c106(:,4,1))/mean(c106(:,2,1));
ci_(1,:) = bootci(nboot,{divisionmean, c106(:,3,1),c106(:,4,1),c106(:,2,1) },'alpha',alpha);

ach3(2) = mean(c66(:,3,1)+c66(:,4,1))/mean(c66(:,2,1));
ci_(2,:) = bootci(nboot,{divisionmean, c66(:,3,1),c66(:,4,1),c66(:,2,1) },'alpha',alpha);

ach3(3) = mean(c44(:,3,1)+c44(:,4,1))/mean(c44(:,2,1));
ci_(3,:) = bootci(nboot,{divisionmean, c44(:,3,1),c44(:,4,1),c44(:,2,1) },'alpha',alpha);     

ach3(4) = mean(c84(:,3,1)+c84(:,4,1))/mean(c84(:,2,1));
ci_(4,:) = bootci(nboot,{divisionmean, c84(:,3,1),c84(:,4,1),c84(:,2,1) },'alpha',alpha); 

figure('Name',['Simulated AcH3H3 Ratio']);
b=bar([1:4],ach3');
hold on;

for i=1:4
    upbound(:,i)=[-ci_(i,1)+ach3(i)];
    downbound(:,i)=[ci_(i,2)-ach3(i)];
end

% Calculate the number of bars in each group
nbars = size(ach3', 2);
% Get the x coordinate of the bars
x = [];
for i = 1:nbars
    x = [x ; b(i).XEndPoints];
end
% Plot the errorbars
errorbar(x',ach3',upbound',downbound','k','linestyle','none')
set(gca, 'XTickLabel', {'10.6','6.6','4.4','8.4'});
xtickangle(45)
hold off

title(['Simulated AcH3/H3 Ratio'])