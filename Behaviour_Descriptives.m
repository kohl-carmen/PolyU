%% Descriptives: Behaviour
% based on: BasicBeh_for_Publication.m 

% -  looks at average ACC and RT for non-tms trials.
% -  plots Figure 2 a/b


clearvars 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


beh_data_folder='D:\PolyU\TMS\Data\Beh';

Partic=1:31;
rt_min=0;
rt_max=1500;

%% load and sort data
ACCURACY=[];
HV_V=[];
D_V=[];
LV_V=[];

for partic=1:length(Partic)
    partic_char=num2str(Partic(partic));
    if Partic(partic)<10
        partic_char=strcat('0',num2str(Partic(partic)));
    end

    MIP=load(strcat(beh_data_folder,'\7',partic_char,'1\transformed\data_behavior_tms'));
    V5=load(strcat(beh_data_folder,'\7',partic_char,'0\transformed\data_behavior_tms'));     
   
     for file={'MIP', 'V5'};
         eval([file{1} '=' file{1} '.data.behavior;'])
     end
    
     
    %combine and make a data matrix that I'm more comfi with
    %     session         1=MIP   0=V5
    %     trial           trial number
    %     rt              RT in ms
    %     tms             0=no 1=stim
    %     tms_cond        0=no 1=MIP 2=V5
    %     accuracy_bolton 1=corr 0=err nan=empty/distractor
    %     accuracy_abs    1=corr 0=error (including emppty/distractor)
    %     accuracy_def    1=corr 0=err 2=distractor 3=empty
    %     d_v             distractor value(prob*rew)
    %     lv_v            low value option value(prob*rew)
    %     hv_v            high value option value(prob*rew)
    %     all_v           value(prob*rew)of each option - 1st column: HV, 2nd column: LV, 3rd column: D
    %     all_prob        probability of each option - 1st column: HV, 2nd column: LV, 3rd column: D
    %     all_mag         magnitude of each option - 1st column: HV, 2nd column: LV, 3rd column: D
    %     position        positions of stims on screen - 1st column: HV, 2nd column: LV, 3rd column: D - 1=top left, 2=top right, 3=bottom left, 4=bottom right
    
    %% caps ones for all partics
    session=[ones(length(MIP.trial{:}),1); zeros(length(V5.trial{:}),1)];%1=MIP, 0=V5;
    trial=[MIP.trial{:};V5.trial{:}];  
    tms=[MIP.tms{:};V5.tms{:}];  
    rt=[MIP.RT{:};V5.RT{:}];
    accuracy_bolton=[MIP.accuracy{:};V5.accuracy{:}]; 
    accuracy_abs=accuracy_bolton; 
    accuracy_abs(isnan(accuracy_abs))=0;  
    d_v=[MIP.D{:};V5.D{:}]; 
    lv_v=[MIP.LV{:};V5.LV{:}]; 
    hv_v=[MIP.HV{:};V5.HV{:}];
    all_v=[MIP.vals{:};V5.vals{:}];
    all_prob=[MIP.probs{:};V5.probs{:}];
    all_mag=[MIP.rews{:};V5.rews{:}];
    position=[MIP.pos_rews{:};V5.pos_rews{:}];   
    
    %for distractor location
     locations={[1 3], [2 4]}; %ipsi=[1 3]=left
     hv_loc=zeros(size(hv_v));
     hv_loc(position(:,1)==2 |position(:,1)==4)=1;
     lv_loc=zeros(size(hv_v));
     lv_loc(position(:,2)==2 |position(:,2)==4)=1;
     di_loc=zeros(size(hv_v));
     di_loc(position(:,3)==2 |position(:,3)==4)=1;
     dloc_for_glm=di_loc;
     
     %value percentage on contra
     clear perc
     for row=1:length(all_v)
         perc(row,:)=sum(all_v(row,position(row,:)==2 |position(row,:)==4))/sum(all_v(row,:));
     end
    
    ACCURACY(:,partic)=(accuracy_bolton);
    HV_V(:,partic)=(hv_v);
    D_V(:,partic)=d_v;
    LV_V(:,partic)=lv_v;
    POSITION(:,:,partic)=position;
    SESSION(:,partic)=session;
    TMS(:,partic)=tms;
    RT(:,partic)=rt;

    
     %% EXCLUDE TRIALS 
    rt_keep=rt;
    accuracy_deletion=(isnan(accuracy_bolton));%delete all trials where distractor or empty was chosen
    tms_deletion=zeros(length(accuracy_bolton),1);

   
    session(accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max |tms_deletion)=[];        
    trial(accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
    rt(accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
    tms( accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
    accuracy_bolton( accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
    accuracy_abs(accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];            
    d_v( accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
    lv_v( accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
    hv_v( accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
    all_v( accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion,:)=[];        
    all_prob( accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion,:)=[];        
    all_mag( accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion,:)=[];        
    position(accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion,:)=[];   
    perc(accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion,:)=[];  
    dloc_for_glm(accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion,:)=[];  
      
      
    DELETE(:,partic)=(tms_deletion | accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max);
    %make congruency matrix (if HV is better than LV in both prob and
    %mag==1)
    Congr=zeros(size(trial));
    Congr(all_prob(:,1)>all_prob(:,2) & all_mag(:,1)>all_mag(:,2))=1;
    Congr_cont=(all_prob(:,1)-all_prob(:,2)).*(all_mag(:,1)-all_mag(:,2));
end

ACCall=[];
RTall=[];
for partic=1:length(Partic)
    ACC_MT(partic)=mean(ACCURACY(SESSION(:,partic)==0 & DELETE(:,partic)==0 & TMS(:,partic)==0,partic));
    ACC_MIP(partic)=mean(ACCURACY(SESSION(:,partic)==1 & DELETE(:,partic)==0 & TMS(:,partic)==0,partic));
   
    RT_MT(partic)=mean(RT(SESSION(:,partic)==0 & DELETE(:,partic)==0 & TMS(:,partic)==0,partic));
    RT_MIP(partic)=mean(RT(SESSION(:,partic)==1 & DELETE(:,partic)==0 & TMS(:,partic)==0,partic));
end

%% Report descriptives
fprintf('Descriptives:\n')
fprintf('Non-tms accuracy (across sessions): %2.2f (SD=%2.2f) \n',mean([ACC_MT,ACC_MIP]).*100, std([ACC_MT,ACC_MIP]).*100)
fprintf('Non-tms RT (across sessions): %2.2f (SD=%2.2f) \n',mean([RT_MT,RT_MIP]), std([RT_MT,RT_MIP]))

fprintf('Differences between sessions in non-tms trials:\n')
[binary,p,ci,tstat]=ttest(ACC_MT,ACC_MIP);
Mean1= mean(ACC_MT);
Mean2= mean(ACC_MIP);
SD1 = std(ACC_MT);
SD2 = std(ACC_MIP);
D= (Mean2-Mean1)/(sqrt(((SD1)^2 +(SD2)^2)/2));
fprintf('t(%i) = %2.2f, p = %1.3f, D = %2.2f\n', tstat.df, tstat.tstat, p,D)
[binary,p,ci,tstat]=ttest(RT_MT,RT_MIP);
Mean1= mean(RT_MT);
Mean2= mean(RT_MIP);
SD1 = std(RT_MT);
SD2 = std(RT_MIP);
D= (Mean2-Mean1)/(sqrt(((SD1)^2 +(SD2)^2)/2));
fprintf('t(%i) = %2.2f, p = %1.3f, D = %2.2f\n', tstat.df, tstat.tstat, p,D)



%% Figure 2 a, b
x=(ACC_MT);
SEM = std(x)/sqrt(length(x));               % Standard Error
ts = tinv([0.025  0.975],length(x)-1);      % T-Score
CI = mean(x) + ts*SEM;                      % Confidence Intervals
error_ACC_MT=[SEM];
x=(ACC_MIP);
SEM = std(x)/sqrt(length(x));               % Standard Error
ts = tinv([0.025  0.975],length(x)-1);      % T-Score
CI = mean(x) + ts*SEM;                      % Confidence Intervals
error_ACC_MIP=[SEM];
x=(RT_MT);
SEM = std(x)/sqrt(length(x));               % Standard Error
ts = tinv([0.025  0.975],length(x)-1);      % T-Score
CI = mean(x) + ts*SEM;                      % Confidence Intervals
error_RT_MT=[SEM];
x=(RT_MIP);
SEM = std(x)/sqrt(length(x));               % Standard Error
ts = tinv([0.025  0.975],length(x)-1);      % T-Score
CI = mean(x) + ts*SEM;                      % Confidence Intervals
error_RT_MIP=[SEM];

clf
subplot(1,2,1)
hold on
bar([mean(ACC_MT), mean(ACC_MIP)])
errorbar([1], mean(ACC_MT), error_ACC_MT, 'k', 'linestyle', 'none');
   
errorbar([2], mean(ACC_MIP), error_ACC_MIP, 'k', 'linestyle', 'none');
ylim([0 1])
set(gca,'YTick', [0 .5 1])
set(gca,'XTick', [1 2])

subplot(1,2,2)
hold on
bar([mean(RT_MT), mean(RT_MIP)])
errorbar([1], mean(RT_MT), error_RT_MT, 'k', 'linestyle', 'none');
   
errorbar([2], mean(RT_MIP), error_RT_MIP, 'k', 'linestyle', 'none');
ylim([0 1000])
set(gca,'YTick', [0 500 1000])
set(gca,'XTick', [1 2])

%save figure
cd('J:\PolyU\TMS\Paper')

print -depsc basic_beh






