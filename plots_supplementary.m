%%
%Fig.S2 - A
cd('simulations/16C_R1_S05')

figure
for i=1:100
    load(strcat('iter',num2str(i),'_C16_R1'))
    hold on
    plot(evolution(:,1),evolution(:,2)./12040487+1,'Color','k')
end

cd('simulations/16C_R1_S1')
for i=1:100
    load(strcat('iter',num2str(i),'_C16_R1'))
    hold on
    plot(evolution(:,1),evolution(:,2)./12040487+1,'Color','r')
end

cd('simulations/16C_R1_S3')
for i=1:100
    load(strcat('iter',num2str(i),'_C16_R1'))
    hold on
    plot(evolution(:,1),evolution(:,2)./12040487+1,'Color',[1 140/255 0])
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
%Fig.S2-B
%computation of completion times
for i=1:100
    load(strcat('iter',num2str(i),'_C16_R1'))
    c=evolution(:,2)./12040487+1;
    [i2,j2]=(min(abs(c-2)));
    [i4,j4]=(min(abs(c-4)));
    [i8,j8]=(min(abs(c-8)));
    [i16,j16]=(min(abs(c-16)));
    mean_all(i,:)=[evolution(j2,1) evolution(j4,1) evolution(j8,1) evolution(j16,1)];
end
CompTime=mean(mean_all);
%%
%Fig.S2-C
figure
subplot(1,4,1)
hold on
cd('simulations/16C_R1_S05')
for i=1:100
    load(strcat('iter',num2str(i),'_C16_R1'))
    plot(evolution(:,1),evolution(:,5),'Color',[204/255 204/255 1])
    plot(evolution(:,1),evolution(:,9),'Color',[1 1 150/255])
    plot(evolution(:,1),evolution(:,6)+evolution(:,7)+evolution(:,8),'Color',[0.76 0.87 0.78])
end
box on
set(gca,'ylim',[0 22000],'fontsize',12)

subplot(1,4,2)
hold on
cd('simulations/16C_R1_S1')
for i=1:100
    load(strcat('iter',num2str(i),'_C16_R1'))
    plot(evolution(:,1),evolution(:,5),'Color',[204/255 204/255 1])
    plot(evolution(:,1),evolution(:,9),'Color',[1 1 150/255])
    plot(evolution(:,1),evolution(:,6)+evolution(:,7)+evolution(:,8),'Color',[0.76 0.87 0.78])
end
box on
set(gca,'ylim',[0 22000],'fontsize',12)

subplot(1,4,3)
hold on
cd('simulations/16C_R1_S3')
for i=1:100
    load(strcat('iter',num2str(i),'_C16_R1'))
    plot(evolution(:,1),evolution(:,5),'Color',[204/255 204/255 1])
    plot(evolution(:,1),evolution(:,9),'Color',[1 1 150/255])
    plot(evolution(:,1),evolution(:,6)+evolution(:,7)+evolution(:,8),'Color',[0.76 0.87 0.78])
end
box on
set(gca,'ylim',[0 22000],'fontsize',12)

subplot(1,4,4)
hold on
cd('simulations/16C_R0_S3')
for i=1:100
    load(strcat('iter',num2str(i),'_C16_R0'))
    plot(evolution(:,1),evolution(:,5),'Color',[204/255 204/255 1])
    plot(evolution(:,1),evolution(:,9),'Color',[1 1 150/255])
    plot(evolution(:,1),evolution(:,6)+evolution(:,7)+evolution(:,8),'Color',[0.76 0.87 0.78])
end
box on
set(gca,'ylim',[0 22000],'fontsize',12)

%%
%%
%Fig.S2-D
figure
subplot(1,2,1)
cd('simulations/16C_R1_S05')

for i=1:100
    load(strcat('iter',num2str(i),'_C16_R1'))
    hold on;box on
    plot(evolution(:,1),evolution(:,3),'k')
end

cd('simulations/16C_R0_S3')
for i=1:100
    load(strcat('iter',num2str(i),'_C16_R0'))
    hold on;box on
    plot(evolution(:,1),evolution(:,3),'b')
end
%ylabel('Total Firing Propensity')


subplot(1,2,2)
cd('simulations/16C_R1_S05')
for i=1:100
    load(strcat('iter',num2str(i),'_C16_R1'))
    hold on;box on
    plot(evolution(:,1),evolution(:,4),'k')
end
cd('simulations/16C_R0_S3')
for i=1:100
    load(strcat('iter',num2str(i),'_C16_R0'))
    hold on;box on
    plot(evolution(:,1),evolution(:,4),'b')
end

set(gca,'fontsize',12)

%%
%Fig.S2 - E
cd('simulations/16C_R1_S05')
figure
for i=1:100
    load(strcat('iter',num2str(i),'_C16_R1'))
    hold on
    plot(evolution(:,1),evolution(:,2)./12040487+1,'Color','k')
end

cd('simulations/16C_R0_S3')
for i=1:100
    load(strcat('iter',num2str(i),'_C16_R0'))
    hold on
    plot(evolution(:,1),evolution(:,2)./12040487+1,'Color','c')
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
%Figure S3

cd('simulations/16C_R1_S05')
files=dir('*.mat');
for i=1:size(files,1)
    load(strcat('iter',num2str(i),'_C16_R1'))
    for j=1:903
        noc{1,1}(i,j)=numel(find(OS(:,j)==0));
    end
end

cd('simulations/16C_R1_S1')
files=dir('*.mat');
for i=1:size(files,1)
    load(strcat('iter',num2str(i),'_C16_R1'))
    for j=1:903
        noc{2,1}(i,j)=numel(find(OS(:,j)==0));
    end
end

cd('simulations/16C_R1_S3')
files=dir('*.mat');
for i=1:size(files,1)
    load(strcat('iter',num2str(i),'_C16_R1'))
    for j=1:903
        noc{3,1}(i,j)=numel(find(OS(:,j)==0));
    end
end

cd('simulations/16C_R0_S3')
files=dir('*.mat');
for i=1:size(files,1)
    load(strcat('iter',num2str(i),'_C16_R0'))
    for j=1:903
        noc{4,1}(i,j)=numel(find(OS(:,j)==0));
    end
end

%plots of normalized means
for j=1:size(noc,1)
    gen_mean(j)=mean(noc{j}(:));    
    raw_mean1(j,:)=mean(noc{j}(:,2:409));
    raw_mean2(j,:)=mean(noc{j}(:,[412:516 519:725]));
    raw_mean3(j,:)=mean(noc{j}(:,[728:810 813:902]));
    norm_mean1(j,:)=raw_mean1(j,:)./gen_mean(j);
    norm_mean2(j,:)=raw_mean2(j,:)./gen_mean(j);
    norm_mean3(j,:)=raw_mean3(j,:)./gen_mean(j);
end

cd('simulations')
load('ori_locations.mat')

figure
subplot(1,3,1)
hold on;box on
plot(loc1,norm_mean1(1:3,:))
xlabel('Chromosome 1')
set(gca,'xlim',[0 6*10^6],'ylim',[0 5])

subplot(1,3,2)
hold on; box on
plot(loc2,norm_mean2(1:3,:))
xlabel('Chromosome 2')
set(gca,'xlim',[0 6*10^6],'ylim',[0 5])

subplot(1,3,3)
hold on;box on;
plot(loc3,norm_mean3(1:3,:))
xlabel('Chromosome 3')
set(gca,'xlim',[0 6*10^6],'ylim',[0 12])

figure
subplot(1,3,1)
hold on;box on
plot(loc1,norm_mean1([1 4],:))
xlabel('Chromosome 1')
set(gca,'xlim',[0 6*10^6],'ylim',[0 5])

subplot(1,3,2)
hold on; box on
plot(loc2,norm_mean2([1 4],:))
xlabel('Chromosome 2')
set(gca,'xlim',[0 6*10^6],'ylim',[0 5])

subplot(1,3,3)
hold on;box on;
plot(loc3,norm_mean3([1 4],:))
xlabel('Chromosome 3')
set(gca,'xlim',[0 6*10^6],'ylim',[0 12])

%%
%Figure S4

cd('simulations/2C_16C_R1_S05')
files=dir('*.mat');
for i=1:size(files,1)
    load(files(i).name)
    for j=1:903
      noc2{1,1}(i,j)=numel(find(OS2(:,j)==0));
      noc16{1,1}(i,j)=numel(find(OS(:,j)==0));
    end
end

rawI_16=noc16{1,1}(:,[2:409]);
rawII_16=noc16{1,1}(:,[412:516 519:725]);
rawIII_16=noc16{1,1}(:,[728:810 813:902]);

rawI_2=noc2{1,1}(:,[2:409]);
rawII_2=noc2{1,1}(:,[412:516 519:725]);
rawIII_2=noc2{1,1}(:,[728:810 813:902]);

%p=randi(100,4,1);
p=[39 72 73 95];
ind=[1;4;7;10];
for i=1:4
    subplot(4,3,ind(i))
    hold off;box on;
    plot(loc1,rawI_16(p(i),:),'.-')
    %plot(loc1,rawI_2(p(i),:),'.-')
    %title(sprintf('Simulation %s',num2str(p(i))))
    %xlabel('Chromosome I');ylabel('Signal ratio')
    set(gca,'xlim',[0 6*10^6])
    set(gca,'ylim',[0 200])
end

ind=[2;5;8;11];
for i=1:4
    subplot(4,3,ind(i))
    hold off;box on;
    plot(loc2,rawII_16(p(i),:),'.-')
    %plot(loc2,rawII_2(p(i),:),'.-')
    %title(sprintf('Simulation %s',num2str(p(i))))
    %xlabel('Chromosome I');ylabel('Signal ratio')
    set(gca,'xlim',[0 6*10^6])
    set(gca,'ylim',[0 200])
end

ind=[3;6;9;12];
for i=1:4
    subplot(4,3,ind(i))
    hold off;box on;
    plot(loc3,rawIII_16(p(i),:),'.-')
    %plot(loc3,rawIII_2(p(i),:),'.-')
    %title(sprintf('Simulation %s',num2str(p(i))))
    %xlabel('Chromosome I');ylabel('Signal ratio')
    set(gca,'xlim',[0 6*10^6])
    set(gca,'ylim',[0 300])
end










