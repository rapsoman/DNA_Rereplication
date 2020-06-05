function [Tfire,TPR,TSR,TSL,OS,evolution,lambdacurrent]=rereplicationalg(copies,redistr)
% -------INTPUTS--------------------
% copies: genome amplification parameter (C)
% redistr: binary variable indicating which model variant to use (1: LF, 0: UF)
% -------OUTPUTS--------------------
% Tfire -> origin firing time
% TPR -> time that an origin is passively replicated
% TSR/TSL -> time that an origin stops to replicate in its Right/Left
% OS -> Origin State (PreR->0, RB->2, RR->1, RL->-1, PostR->3, PassR->4, not appear->5)
% evolution -> [timepoints, total replicated DNA, Total Firing Propensity, average TFP, number of origins in PreR, RB, RR, RL, PassR, PostR]
warning off
load InputData
tic
Location=alldata(:,1)'; 
Speed=alldata(:,2)'; %1000 b/min
Efficiency=alldata(:,3)'; %origin efficiencies
oris=length(Efficiency);
GenomeLength=12040487;
for i=1:oris
    if Efficiency(i)~=0 %efficiency=0 for dummy origins 
        lambda(i)=(-1/20)*log(1-Efficiency(i)); %firing propensity
    else
        lambda(i)=0;
    end
end

Tfire=FiringTimes(lambda); %compute firing time

%initialize variables
TPR=inf*ones(size(Tfire));
TSR=inf*ones(size(Tfire));
TSL=inf*ones(size(Tfire));
OS=zeros(size(Tfire));
OS(Speed==0)=2; %for dummy origins
Tfire(Speed==0)=0;%same
OSnew=OS;
forkid(1)=1;
lambdacurrent=lambda;
DNAcontent=0;
endtimeall=0; %time an event will happen
evolution=[];

disp('Running re-replication model...')
while DNAcontent/GenomeLength<(copies-1)
    DNAcontentnew=0;
    [forks, oris]=size(Tfire);
    eventfork=[]; %
    for fork=1:forks %for all forks
        %first compute origin states 
        preroris=find(OS(fork,:)==0);
        nonappear=find(OS(fork,:)==5);
        RBoris=find(OS(fork,:)==2);
        RLoris=find(OS(fork,:)==-1);
        RRoris=find(OS(fork,:)==1);
        nonactive=sort([preroris nonappear]);
        active=sort([RBoris RLoris RRoris]);
               
        %transitions: PreR|NonAppear -> PassR
        for j=1:length(nonactive) %for all non-active origins (PreR or NonAppear)
            ori=nonactive(j);
            prev=active(find(active<ori,1,'last')); %find the previous active 
            next=active(find(active>ori,1,'first')); %find the next active
            TPRleft=((Location(ori)-Location(prev))/Speed(prev))+...
                Tfire(fork,prev); %TPR from the previous
            TPRright=((Location(next)-Location(ori))/Speed(next))+...
                Tfire(fork,next); %TPR from the next
            TPR(fork,ori)=min(TPRleft,TPRright); %Time of Passive Replication is the minimum of that
        end
        
        %transitions: RB -> RR & RL
        for j=2:length(RBoris)-1 %for all RB origins (except the dummy ones)
            ori=RBoris(j);
            prev=active(max(find(active<ori))); %find previous active
            next=active(min(find(active>ori))); %find next active
            %Time it will stop replicating to the left
            TSL(fork,ori)=(Location(ori)-Location(prev)+Tfire(fork,ori)*Speed(ori)+...
                Tfire(fork,prev)*Speed(prev))/(Speed(ori)+Speed(prev));
            %Time it will stop replicating to the right
            TSR(fork,ori)=(Location(next)-Location(ori)+Tfire(fork,next)*Speed(next)+...
                Tfire(fork,ori)*Speed(ori))/(Speed(ori)+Speed(next));
        end
        
        %transition: RR -> PostR 
        for j=1:length(RRoris)
            ori=RRoris(j);
            next=active(min(find(active>ori)));
            TSR(fork,ori)=(Location(next)-Location(ori)+Tfire(fork,next)*Speed(next)+Tfire(fork,ori)*Speed(ori))/(Speed(ori)+Speed(next));
        end
        %transition: RL -> PostR
        for j=1:length(RLoris)
            ori=RLoris(j);
            prev=active(max(find(active<ori)));
            TSL(fork,ori)=(Location(ori)-Location(prev)+Tfire(fork,ori)*Speed(ori)+Tfire(fork,prev)*Speed(prev))/(Speed(ori)+Speed(prev));
        end
        
        %choose possible events to happen next
        oristofire=preroris(find(Tfire(fork,preroris)>endtimeall));%candidates for firing
        nexttoPassR=nonactive(find(TPR(fork,nonactive)>endtimeall));%candidates for passive replication
        RR=sort([RBoris RRoris]);
        nexttoTSR=RR(find(TSR(fork,RR)>endtimeall));%candidates to stop on the right
        RL=sort([RBoris RLoris]);
        nexttoTSL=RL(find(TSL(fork,RL)>endtimeall));%candidates to stop on the left
        %candidate time an event might happen
        eventtimes=[TPR(fork,nexttoPassR) TSR(fork,nexttoTSR) TSL(fork,nexttoTSL) Tfire(fork,oristofire)];
        eventtype=[ones(size(TPR(fork,nexttoPassR))) 2*ones(size(TSR(fork,nexttoTSR))) 3*ones(size(TSL(fork,nexttoTSL))) 4*ones(size(Tfire(fork,oristofire)))];
        origins=[nexttoPassR nexttoTSR nexttoTSL oristofire];
        endtime(fork)=min(eventtimes); %time of first possible event 
        temp=find(eventtimes==endtime(fork));%index of events
        %eventfork=[which_fork which_origin which_event]
        for i=1:length(temp)
            eventfork=[eventfork; [fork origins(temp(i))  eventtype(temp(i))]];
        end
    end
    endtimeallprev=endtimeall;
    endtimeall=min(endtime);%update minimum time
    forkevent=find(endtime==endtimeall);%fork where the event will take place
    eventservice=eventfork(find(eventfork(:,1)==forkevent),:);%what will happen

    %update DNA content
    for fork=1:forks %for all forks
        for oriid=1:oris %and all origins
            %if origin is in RB, RR, RL or PostR
            if ((OS(fork,oriid)==2 || OS(fork,oriid)==1 || OS(fork,oriid)==-1 || OS(fork,oriid)==3) && Speed(oriid)~=0)
                %update DNA content
                DNAcontentnew=DNAcontentnew + min(Speed(oriid)*(endtimeall-Tfire(fork,oriid)),Speed(oriid)*(TSL(fork,oriid)-Tfire(fork,oriid)))*(endtimeall>Tfire(fork,oriid))+...
                                            + min(Speed(oriid)*(endtimeall-Tfire(fork,oriid)),Speed(oriid)*(TSR(fork,oriid)-Tfire(fork,oriid)))*(endtimeall>Tfire(fork,oriid));
            end
        end
    end
    
    if endtimeall<endtimeallprev
        error('Back Time')
    end
    if DNAcontentnew<DNAcontent
        error('Back Distance')
    end
    DNAcontent=DNAcontentnew;

    [events, data]=size(eventservice); %how many events will happen
    
    %perform events and update states 
    for i=1:events %for all events
        fork=eventservice(i,1); %which fork
        ori=eventservice(i,2); %which origin
        forkID=forkid(fork);%assign id to fork where the event happens
        child1=find(forkid==2*forkID);%the children of that
        child2=find(forkid==2*forkID+1);
        
        if eventservice(i,3)==1 %if passive replication
            OSnew(fork,ori)=4; %new OS state = PassR
            TPR(fork,ori)=endtimeall; %update time
            %update states of children to PreR
            OSnew(child1,ori)=0*(OS(child1,ori)==5)+4*(OS(child1,ori)==4);
            OSnew(child2,ori)=0*(OS(child2,ori)==5)+4*(OS(child2,ori)==4);
            %update lambda
            lambdanew=redistribution(lambdacurrent,lambda,OS,OSnew,redistr);
            Tfirenew=FiringTimes(lambdanew); %new firing times
            [forks, ori2alter]=find(OSnew==0); %find new PreR origins
            for j=1:length(forks) %for all PreR origins
                %update Tfire 
                Tfire(forks(j),ori2alter(j))=endtimeall+Tfirenew(forks(j),ori2alter(j));
            end
            lambdacurrent=lambdanew;
        end
        
        if eventservice(i,3)==4 %if firing (PreR->RB)
            OSnew(fork,ori)=2; %new OS state = RB
            TPR(fork,ori)=-1; %TPR=-1 -> origin has fired
            if ~isempty([child1 child2]) %if there are children
                %update children state to PreR
                OSnew(child1,ori)=0*(OS(child1,ori)==5)+4*(OS(child1,ori)==4);
                OSnew(child2,ori)=0*(OS(child2,ori)==5)+4*(OS(child2,ori)==4);
            else
                %add new lines to state and time matrices
                forkid(end+1)=2*forkID;%for first child
                OS(end+1,:)=5; %all "non appear"
                OSnew(end+1,:)=5;
                %new times - all Inf
                TPR(end+1,:)=inf;
                TSR(end+1,:)=inf;
                TSL(end+1,:)=inf;
                Tfire(end+1,:)=inf;
                lambdacurrent(end+1,:)=0;
                OSnew(end,ori)=0;%OS of fired ori=PreR
                forkid(end+1)=2*forkID+1; %for second child
                OS(end+1,:)=5; 
                OSnew(end+1,:)=5;
                TPR(end+1,:)=inf;
                TSR(end+1,:)=inf;
                TSL(end+1,:)=inf;
                Tfire(end+1,:)=inf;
                lambdacurrent(end+1,:)=0;
                OSnew(end,ori)=0;
            end
            
                lambdanew=redistribution(lambdacurrent,lambda,OS,OSnew,redistr);

           %update firing times according to new lambda
           Tfirenew=FiringTimes(lambdanew);
           [forks, ori2alter]=find(OSnew==0);
            %update firing time of children
            for j=1:length(forks)
                Tfire(forks(j),ori2alter(j))=endtimeall+Tfirenew(forks(j),ori2alter(j));
            end
            lambdacurrent=lambdanew;
        end
        
        if eventservice(i,3)==2 %if event = stop to the right
            OSnew(fork,ori)=-1*(OS(fork,ori)==2)+3*(OS(fork,ori)==1);
        end
        
        if eventservice(i,3)==3 %if event = stop to the left
            OSnew(fork,ori)=1*(OS(fork,ori)==2)+3*(OS(fork,ori)==-1);

        end
    end
    %check the events
    changes=find(OSnew~=OS); %find location of events
    A=[find(OS(changes)==5 & OSnew(changes)~=0);
    find(OS(changes)==0 & (OSnew(changes)~=2 & OSnew(changes)~=4));
    find(OS(changes)==4  | OS(changes)==3);
    find(OS(changes)==2 & (OSnew(changes)~=1 & OSnew(changes)~=-1));
    find(OS(changes)==1 & OSnew(changes)~=3);
    find(OS(changes)==-1 & OSnew(changes)~=3)];  
    if ~isempty(A)
        error('Wrong Transition')
    end
    
    OSnew(1:end,find(Speed==0))=2;
    Tfire(1:end,find(Speed==0))=0;
    OS=OSnew; %update state
    
   
    %--------------- Compute Outputs --------------------------------------
    %count OS states
    NoPreR=length(find(OS==0));
    NoRB=length(find(OS==2));
    NoRR=length(find(OS==1));
    NoRL=length(find(OS==-1));
    NoPassR=length(find(OS==4));
    NoPostR=length(find(OS==3));
    %update TF
    [forks oris]=size(Tfire);
    
    preroris=find(OS==0);
    TFP=sum(lambdacurrent(preroris));
    meanTFP=TFP/NoPreR;
    endtime;
    evolution=[evolution;endtimeall DNAcontent TFP meanTFP NoPreR NoRB NoRR NoRL NoPassR NoPostR];
    if DNAcontent/GenomeLength<(copies-1)
    end
end
toc
