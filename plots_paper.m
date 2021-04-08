%Figure 1D
%plots of DNA content
cd('simulations/16C_R1_S05')
figure
for i=1:100
    load(strcat('iter',num2str(i),'_C16_R1'))
    hold on
    plot(evolution(:,1),evolution(:,2)./12040487+1,'Color','k')
end

xlabel('Time (min)')
ylabel('Total DNA content (C)')
plot([0 300],[1 1],':k')
plot([0 300],[2 2],':k')
plot([0 300],[4 4],':k')
plot([0 300],[8 8],':k')
plot([0 300],[16 16],':k')
box on
set(gca,'fontsize',12,'ytick',[1 2 4 8 16])
%%
%Figure 1E
figure
hold on
cd(folder_16C_R1_S05)
for i=1:100
    load(strcat('iter',num2str(i),'_C16_R1'))
    plot(evolution(:,1),evolution(:,5),'Color',[204/255 204/255 1])
    plot(evolution(:,1),evolution(:,9),'Color',[1 1 150/255])
    plot(evolution(:,1),evolution(:,6)+evolution(:,7)+evolution(:,8),'Color',[0.76 0.87 0.78])
end
xlabel('Time (min)')
ylabel('Number of origins')
%legend('ORIs in PreR','ORIs in PassR','ORIs in RB,RR,RL')
box on
set(gca,'xlim',[0 300],'ylim',[0 18000],'fontsize',12)

%%
%Figure 2
cd('simulations/peaks')
load exp_peaks.mat
%denoise and find peaks from experimental data
x=80;
YS1 = mslowess(loc1_exp,int1,'order',2,'span',x,'SHOWPLOT',false);
YS2 = mslowess(loc2_exp,int2,'order',2,'span',x,'SHOWPLOT',false);
YS3 = mslowess(loc3_exp,int3,'order',2,'span',x,'SHOWPLOT',false);
[P1, pfwhh1] = mspeaks(loc1_exp,YS1,'denoising',true,'HeightFilter',1,'OverSegmentationFilter',40000,'SHOWPLOT',0);
[P2, pfwhh2] = mspeaks(loc2_exp,YS2,'denoising',true,'HeightFilter',1,'OverSegmentationFilter',40000,'SHOWPLOT',0);
[P3, pfwhh3] = mspeaks(loc3_exp,YS3,'denoising',true,'HeightFilter',1,'OverSegmentationFilter',40000,'SHOWPLOT',0);

subplot(2,3,1)
hold on;box on;
plot(loc1_exp,int1,'b')
plot(loc1_exp,YS1,'k')
plot([P1(:,1) P1(:,1)],[0 3],':k')
set(gca,'xlim',[0 6*10^6],'ylim',[0.5 3])

subplot(2,3,2)
hold on;box on;
plot(loc2_exp,int2,'b')
plot(loc2_exp,YS2,'k')
plot([P2(:,1) P2(:,1)],[0 3],':k')
set(gca,'xlim',[0 6*10^6],'ylim',[0.5 3])

subplot(2,3,3)
hold on;box on;
plot(loc3_exp,int3,'b')
plot(loc3_exp,YS3,'k')
plot([P3(:,1) P3(:,1)],[0 3],':k')
set(gca,'xlim',[0 6*10^6],'ylim',[0.5 3])

load sim_peaks.mat
P1s = mspeaks(loc1,norm_mean1,'denoising',false,'multiplier',1,'HeightFilter',1,'OverSegmentationFilter',40000,'SHOWPLOT',0);
P2s = mspeaks(loc2,norm_mean2,'denoising',false,'multiplier',1,'HeightFilter',1,'OverSegmentationFilter',40000,'SHOWPLOT',0);
P3s = mspeaks(loc3,norm_mean3,'denoising',false,'multiplier',2,'HeightFilter',1,'OverSegmentationFilter',40000,'SHOWPLOT',0);

subplot(2,3,4)
hold on
plot(loc1,norm_mean1,'k')
plot([P1s(:,1) P1s(:,1)],[0 10],':k')
xlabel('Chromosome 1')
set(gca,'xlim',[0 6*10^6],'ylim',[0 5])

subplot(2,3,5)
hold on
plot(loc2,norm_mean2,'k')
plot([P2s(:,1) P2s(:,1)],[0 10],':k')
xlabel('Chromosome 2')
set(gca,'xlim',[0 6*10^6],'ylim',[0 5])

subplot(2,3,6)
hold on
plot(loc3,norm_mean3,'k')
plot([P3s(:,1) P3s(:,1)],[0 10],':k')
xlabel('Chromosome 3')
set(gca,'xlim',[0 6*10^6],'ylim',[0 10])
%plot common peaks
subplot(2,3,1)
hold on;box on;
plot([P1(cntr1e==1,1) P1(cntr1e==1,1)],[0 5],'-.r')

subplot(2,3,4)
hold on;box on;
plot([P1s(cntr1==1,1) P1s(cntr1==1,1)],[0 5],'-.r')

subplot(2,3,2)
hold on;box on;
plot([P2(cntr2e==1,1) P2(cntr2e==1,1)],[0 5],'-.r')

subplot(2,3,5)
hold on;box on;
plot([P2s(cntr2==1,1) P2s(cntr2==1,1)],[0 5],'-.r')

subplot(2,3,3)
hold on; box on;
plot([P3s(cntr3==1,1) P3s(cntr3==1,1)],[0 5],'-.r')

subplot(2,3,6)
hold on; box on;
plot([P3s(cntr3==1,1) P3s(cntr3==1,1)],[0 10],'-.r')

%%
%data for Figure 3
cd('simulations/2C_16C_R1_S05')
files=dir('*.mat');
for i=1:size(files,1)
    load(files(i).name)
    for j=1:903
      noc2{1,1}(i,j)=numel(find(OS2(:,j)==0));
      noc16{1,1}(i,j)=numel(find(OS(:,j)==0));
    end
end

%%
%Figure 3A-B
cd('simulations')
load('ori_locations.mat')
rawI_16=noc16{1,1}(:,[2:409]);
rawII_16=noc16{1,1}(:,[412:516 519:725]);
rawIII_16=noc16{1,1}(:,[728:810 813:902]);

rawI_2=noc2{1,1}(:,[2:409]);
rawII_2=noc2{1,1}(:,[412:516 519:725]);
rawIII_2=noc2{1,1}(:,[728:810 813:902]);

p=[39 72 73 95];
figure;
for i=1:4
    subplot(1,4,i)
    hold on;box on;
    plot(loc1,rawI_2(p(i),:)./2,'.-')
    %title(sprintf('Simulation %s',num2str(p(i))))
    %xlabel('Chromosome I');ylabel('Signal ratio')
    set(gca,'xlim',[3.5*10^6 4.5*10^6])
    set(gca,'ylim',[0 20])
end
figure;
for i=1:4
    subplot(1,4,i)
    hold on;box on;
    plot(loc1,rawI_16(p(i),:)./16,'-')
    %title(sprintf('Simulation %s',num2str(p(i))))
    %xlabel('Chromosome I');ylabel('Signal ratio')
    set(gca,'xlim',[3.5*10^6 4.5*10^6])
    set(gca,'ylim',[0 20])
end

%%
%Figure 3C

rawI_16=noc16{1,1}(:,[2:409]);
histogram(rawI_16(:,272),20,'facecolor','blue')
[f,x]=ksdensity(rawI_16(:,272),'support','positive');
hold on;area(x,f*800,'facecolor','blue')
alpha(0.5)
%%
%Figure 3D-E
rawII_16=noc16{1,1}(:,[412:516 519:725]);

subplot(2,1,1)
histogram(rawII_16(:,132))
[f,x]=ksdensity(rawII_16(:,132));
hold on;plot(x,f*160,'b','linewidth',2)
subplot(2,1,2)
histogram(rawII_16(:,153))
[f,x]=ksdensity(rawII_16(:,153));
hold on;plot(x,f*160,'b','linewidth',2)

figure;hold on;
plot(loc2,mean(rawII_16))
plot(loc2,(rawII_16(70,:)))

%%
%Figure 4
subplot(2,2,1)
[f1,x1]=ksdensity(nuclei(:,1));
[f2,x2]=ksdensity(nuclei(1:1234,2));
hold on;box on;
area(x1,f1)
area(x2,f2,'facecolor','yellow')
alpha(0.5)
%set(gca,'xlim',[0 0.06],'ytick',[],'fontsize',12)
h=legend('Normal replication','Re-replication');
set(h,'fontsize',12, 'box','off')
ylabel('Probability density')

subplot(2,2,2)
gscatter(foci_all, nuclei_all,ind,'by')

subplot(2,2,4)
[f1,x1]=ksdensity(foci(:,1));
[f2,x2]=ksdensity(foci(1:1234,2));
hold on;box on;
area(x1,f1)
area(x2,f2,'facecolor','yellow')
alpha(0.5)
set(gca,'xlim',[0 2.5*10^7],'ytick',[],'fontsize',12)
h=legend('Normal replication','Re-replication');
set(h,'fontsize',12, 'box','off')
ylabel('Probability density')
%%
%Figure 5
%data for Figure 5
cd('simulations/2C_16C_R1_S05')
files=dir('*.mat');
for i=1:size(files,1)
    load(files(i).name)
    for j=1:903
       nop_C16_R1_S05(i,j)=numel(find(OS(:,j)==4));
       nof_C16_R1_S05(i,j)=numel(find(OS(:,j)==2 | OS(:,j)==1 | OS(:,j)==-1 | OS(:,j)==3));
       noc_C16_R1_S05(i,j)=numel(find(OS(:,j)==0));
    end
end
clear i j files index copies redistr evolution evolution2 lambdacurrent fname file

noc={noc_C16_R1_S05};
for j=1:size(noc,1)
    raw1{j}=noc{j}(:,2:409);
    raw2{j}=noc{j}(:,[412:516 519:725]);
    raw3{j}=noc{j}(:,[728:810 813:902]);
end
nof_all=nof_C16_R1_S05;
nof_all(:,[1 410 411 517 518 726 727 811 812 903])=[];
noc_all=noc_C16_R1_S05;
noc_all(:,[1 410 411 517 518 726 727 811 812 903])=[];
nop_all=nop_C16_R1_S05;
nop_all(:,[1 410 411 517 518 726 727 811 812 903])=[];

cd('simulations')
load efficiencies.mat
eff_all=[eff1; eff2; eff3];
%%
%Figure 5A
[f1,x1]=ksdensity(raw2{1,1}(:,45),'support','positive');
[f2,x2]=ksdensity(raw2{1,1}(:,54),'support','positive');
figure;hold on;box on;
area(x1,f1)
area(x2,f2,'facecolor','yellow')
alpha(0.5)
set(gca,'xlim',[-5 200],'ytick',[],'fontsize',12)
h=legend('Ori II-45 FE=62%, median=36','Ori II-54 FE=9%, median=10');
set(h,'fontsize',12, 'box','off')
xlabel('Number of copies')
ylabel('Probability density')

%%
%Figure 5B
scatter(eff_all,mean(nof_all),15,(std(nof_all)./mean(nof_all)));
h2=text(60,8,'\rho = 0.4');
set(h2,'fontsize',12)
h=gca;
h.YScale='log';
set(h,'xlim',[-1 100],'fontsize',12)
box on;
c = colorbar;
c.Label.String = 'Coefficient of variation';
c.Label.FontSize = 12;
xlabel('Efficiency');ylabel('Number of fires (mean)')
%%
%Figure 5C
scatter(eff_all,mean(noc_all),15,(std(noc_all)./mean(noc_all)));
h2=text(60,8,'\rho = 0.4');
set(h2,'fontsize',12)
h=gca;
h.YScale='log';
set(h,'xlim',[-1 100],'fontsize',12)
box on;
c = colorbar;
c.Label.String = 'Coefficient of variation';
c.Label.FontSize = 12;
xlabel('Efficiency');ylabel('Number of copies (mean)')

%%
%Figure 5D
scatter(eff_all,mean(nop_all),15,(std(nop_all)./mean(nop_all)));
h2=text(80,100,'\rho = 0.17');
set(h2,'fontsize',12)
h=gca;
h.YScale='log';
set(h,'xlim',[-1 100],'fontsize',12)
box on;
c = colorbar;
c.Label.String = 'Coefficient of variation';
c.Label.FontSize = 12;
xlabel('Efficiency');ylabel('Number of passive replications (mean)')
%%
%Figure 5E
[f1,x1]=ksdensity(raw1{1,1}(:,287),'support','positive');
figure;hold on;box on;
area(x1,f1)
alpha(0.5)
set(gca,'xlim',[-5 200],'ytick',[],'fontsize',12)
h=legend('Ori I-287, FE=9%, median=48');
set(h,'fontsize',12, 'box','off')
xlabel('Number of copies')
ylabel('Probability density')

%%
%Figure 5G
cd('simulations')
load efficiencies.mat
load ori_locations.mat

cd('simulations/16C_R1_S05')
files=dir('*.mat');
for i=1:size(files,1)
    load(files(i).name)
    for j=1:903
       noc{1,1}(i,j)=numel(find(OS(:,j)==0));
    end
end
clear i j files index copies redistr evolution evolution2 lambdacurrent fname file

cd('simulations/16C_R1_S1')
files=dir('*.mat');
for i=1:size(files,1)
    load(files(i).name)
    for j=1:903
       noc{2,1}(i,j)=numel(find(OS(:,j)==0));
    end
end
clear i j files index copies redistr evolution evolution2 lambdacurrent fname file

cd('simulations/16C_R1_S3')
files=dir('*.mat');
for i=1:size(files,1)
    load(files(i).name)
    for j=1:903
       noc{3,1}(i,j)=numel(find(OS(:,j)==0));
    end
end
clear i j files index copies redistr evolution evolution2 lambdacurrent fname file

for j=1:size(noc,1)
    raw1{j}=noc{j,1}(:,2:409);
    raw2{j}=noc{j,1}(:,[412:516 519:725]);
    raw3{j}=noc{j,1}(:,[728:810 813:902]);
end

loc_all=[loc1;loc2;loc3];
eff_all=[eff1;eff2;eff3];
c=colormap(parula(3));
for k=1:3
    raw_all=[raw1{k} raw2{k} raw3{k}];
    cx=corr(raw_all,'type','spearman');
    [id1,id2]=sort(eff_all,'descend');
    for i=1:10
        hold on
        locx=loc_all-loc_all(id2(i));
        plot(locx(id2(i)-35:id2(i)+35),cx(id2(i)-35:id2(i)+35,id2(i)),'.-','Color',c(k,:))
    end
end
set(gca,'xlim',[-300000 300000],'ylim',[-0.4 1.1],'fontsize',12)
box on
xlabel('Distance from flanking origins (in bases)')
ylabel('Spearman correlation coefficient')
%%
%Figure 5H
[f1,x1]=ksdensity(nof_all(:,287),'bandwidth',0.5);
[f2,x2]=ksdensity(nof_all(:,40),'bandwidth',0.5);
figure;hold on;box on;
area(x1,f1)
area(x2,f2,'facecolor','yellow')
alpha(0.5)
set(gca,'xlim',[0 8],'ytick',[],'fontsize',12)
h=legend('Ori I-287, FE=9%, median=2','Ori I-40, FE=9%, median=0');
set(h,'fontsize',12, 'box','off')
xlabel('Number of fires')
ylabel('Probability density')
%%
%data for Figure 5J
cd('simulations')
load efficiencies.mat
load ori_locations.mat

cd('simulations/16C_R1_S05')
files=dir('*.mat');
for i=1:size(files,1)
    load(files(i).name)
    for j=1:903
       nof{1,1}(i,j)=numel(find(OS(:,j)==2 | OS(:,j)==1 | OS(:,j)==-1 | OS(:,j)==3));
    end
end
clear i j files index copies redistr evolution evolution2 lambdacurrent fname file

cd('simulations/16C_R1_S1')
files=dir('*.mat');
for i=1:size(files,1)
    load(files(i).name)
    for j=1:903
       nof{2,1}(i,j)=numel(find(OS(:,j)==2 | OS(:,j)==1 | OS(:,j)==-1 | OS(:,j)==3));
    end
end
clear i j files index copies redistr evolution evolution2 lambdacurrent fname file

cd('simulations/16C_R1_S3')
files=dir('*.mat');
for i=1:size(files,1)
    load(files(i).name)
    for j=1:903
       nof{3,1}(i,j)=numel(find(OS(:,j)==2 | OS(:,j)==1 | OS(:,j)==-1 | OS(:,j)==3));
    end
end
clear i j files index copies redistr evolution evolution2 lambdacurrent fname file

for j=1:size(nof,1)
    raw1{j}=nof{j,1}(:,2:409);
    raw2{j}=nof{j,1}(:,[412:516 519:725]);
    raw3{j}=nof{j,1}(:,[728:810 813:902]);
end

loc_all=[loc1;loc2;loc3];
eff_all=[eff1;eff2;eff3];
c=colormap(parula(3));
for k=1:3
    raw_all=[raw1{k} raw2{k} raw3{k}];
    cx=corr(raw_all,'type','spearman');
    [id1,id2]=sort(eff_all,'descend');
    for i=1:10
        hold on
        locx=loc_all-loc_all(id2(i));
        plot(locx(id2(i)-35:id2(i)+35),cx(id2(i)-35:id2(i)+35,id2(i)),'.-','Color',c(k,:))
    end
end
set(gca,'xlim',[-300000 300000],'ylim',[-0.4 1.1],'fontsize',12)
box on
xlabel('Distance from flanking origins (in bases)')
ylabel('Spearman correlation coefficient')
%%
%Figure 6A
%load data for Figure 3 (see above)
raw16=noc16{1,1}(:,[2:409 412:516 519:725 728:810 813:902]);
[coeff,score] = pca(raw16);
figure;
biplot(coeff(:,1:2),'scores',score(:,1:2));

%Figure 6B-D
cd('simulations')
load('data_Fig6.mat')

[idx_final]=kmeans(raw16,3,'replicates',100,'display','final');
[ii,jj]=sort(idx_final);
figure;
imagesc(lognorm_all16(jj,:))
figure;
imagesc(lognorm_all2(jj,:))
%%
% statistical analysis for Figure 6B-D
% gap statistic
eva = evalclusters(raw2,'kmeans','gap','KList',[1:20],'B',100,'searchmethod','firstmaxSE');

% Adjusted Rand Index
for i=1:100;[idx_2(:,i)]=kmeans(raw2,3,'replicates',50,'display','off');end
for l=1:100
    for m=1:100
        ri_2(l,m) = rand_index(idx_2(:,l), idx(:,m), 'adjusted');
    end
end
for i=1:100;[idx_16(:,i)]=kmeans(raw16,3,'replicates',50,'display','off');end
for l=1:100
    for m=1:100
        ri_16(l,m) = rand_index(idx_16(:,l), idx(:,m), 'adjusted');
    end
end

boxplot([ri_2(:) ri_16(:)])

for i=1:100;[idx_2(:,i)]=kmeans(raw2,3,'replicates',50,'display','off');end
for i=1:100;[idx_16(:,i)]=kmeans(raw16,3,'replicates',50,'display','off');end
for l=1:100
    for m=1:100
        ri_across(l,m) = rand_index(idx_2(:,l), idx_16(:,m), 'adjusted');
    end
end

%%




