%% Behavioural Analysis for TMS7 Publication
% This is just for Figure 2 a/b
% only looks at average ACC and RT for non-tms trials.

% first part very similar to Beh.m
clearvars 



% Print to ppt (1) or not (0)?    
ppt=1;
%% --- ----- --- %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beh_data_folder='D:\PolyU\TMS\Data\Beh';

Partic=1:31;
exclude_bad_tms=0;
tms_dev_cutoff=5;
stepwise=0;
rt_min=0;%100;
rt_max=1500;%1499;
accs=[];
exclude=[];
magprob_weight=0;


ACCURACY=[];
HV_V=[];
D_V=[];
LV_V=[];

for partic=1:length(Partic)
% partic=partic+1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% GET DATA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    partic_char=num2str(Partic(partic));
    if Partic(partic)<10
        partic_char=strcat('0',num2str(Partic(partic)));
    end

    MIP=load(strcat(beh_data_folder,'\7',partic_char,'1\transformed\data_behavior_tms'));
    V5=load(strcat(beh_data_folder,'\7',partic_char,'0\transformed\data_behavior_tms'));     
%     MIP=load(strcat(beh_data_folder,'\7',partic_char,'1\transformed\data_behavior_tms_with_practice'));
%     V5=load(strcat(beh_data_folder,'\7',partic_char,'0\transformed\data_behavior_tms_with_practice'));   
%     
    
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
%     ACCURACY(partic)=nanmean(accuracy_bolton);
    
    
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
%        
%    ACCURACY_BOLTON(:,partic )=accuracy_bolton;
%    ACC(partic,:)=[partic,mean(accuracy_bolton(~isnan(accuracy_bolton)))];
%    
    
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
    
    ACC_MT1(partic)=mean(ACCURACY(SESSION(:,partic)==0 & DELETE(:,partic)==0 & TMS(:,partic)==1,partic));
    ACC_MIP1(partic)=mean(ACCURACY(SESSION(:,partic)==1 & DELETE(:,partic)==0 & TMS(:,partic)==1,partic));
   
    RT_MT(partic)=mean(RT(SESSION(:,partic)==0 & DELETE(:,partic)==0 & TMS(:,partic)==0,partic));
    RT_MIP(partic)=mean(RT(SESSION(:,partic)==1 & DELETE(:,partic)==0 & TMS(:,partic)==0,partic));
    
    RT_MT1(partic)=mean(RT(SESSION(:,partic)==0 & DELETE(:,partic)==0 & TMS(:,partic)==1,partic));
    RT_MIP1(partic)=mean(RT(SESSION(:,partic)==1 & DELETE(:,partic)==0 & TMS(:,partic)==1,partic));
end

x=(ACC_MT);
SEM = std(x)/sqrt(length(x));               
ts = tinv([0.025  0.975],length(x)-1);     
error_ACC_MT=[SEM];
x=(ACC_MIP);
SEM = std(x)/sqrt(length(x));               
ts = tinv([0.025  0.975],length(x)-1);   
error_ACC_MIP=[SEM];
x=(RT_MT);
SEM = std(x)/sqrt(length(x));            
ts = tinv([0.025  0.975],length(x)-1);   
error_RT_MT=[SEM];
x=(RT_MIP);
SEM = std(x)/sqrt(length(x));               
ts = tinv([0.025  0.975],length(x)-1);      
error_RT_MIP=[SEM];

x=(ACC_MT1);
SEM = std(x)/sqrt(length(x));               
ts = tinv([0.025  0.975],length(x)-1);     
error_ACC_MT1=[SEM];
x=(ACC_MIP1);
SEM = std(x)/sqrt(length(x));               
ts = tinv([0.025  0.975],length(x)-1);   
error_ACC_MIP1=[SEM];
x=(RT_MT1);
SEM = std(x)/sqrt(length(x));            
ts = tinv([0.025  0.975],length(x)-1);   
error_RT_MT1=[SEM];
x=(RT_MIP1);
SEM = std(x)/sqrt(length(x));               
ts = tinv([0.025  0.975],length(x)-1);      
error_RT_MIP1=[SEM];


figure
clf
subplot(1,2,1)
hold on
b = bar([mean(ACC_MT),mean(ACC_MT1); mean(ACC_MIP),mean(ACC_MIP1)]);

x = [];
for i = 1:2
    x(i,:) = b(i).XEndPoints;
end
errorbar([x(1,1)], mean(ACC_MT), error_ACC_MT, 'k', 'linestyle', 'none');
errorbar([x(2,1)], mean(ACC_MT1), error_ACC_MT1, 'k', 'linestyle', 'none');
errorbar([x(1,2)], mean(ACC_MIP), error_ACC_MIP, 'k', 'linestyle', 'none');
errorbar([x(2,2)], mean(ACC_MIP1), error_ACC_MIP1, 'k', 'linestyle', 'none');
ylim([0.4 1])
set(gca,'YTick', [0 .5 1])
set(gca,'XTick', [1 2])
title('Accuracy')
set(gca,'XTicklabels', {'MT','MIP'})


subplot(1,2,2)
hold on
b = bar([mean(RT_MT),mean(RT_MT1); mean(RT_MIP),mean(RT_MIP1)]);
x = [];
for i = 1:2
    x(i,:) = b(i).XEndPoints;
end
errorbar([x(1,1)], mean(RT_MT), error_RT_MT, 'k', 'linestyle', 'none');
errorbar([x(2,1)], mean(RT_MT1), error_RT_MT1, 'k', 'linestyle', 'none');
errorbar([x(1,2)], mean(RT_MIP), error_RT_MIP, 'k', 'linestyle', 'none');
errorbar([x(2,2)], mean(RT_MIP1), error_RT_MIP1, 'k', 'linestyle', 'none');
ylim([400 1000])
set(gca,'YTick', [0 500 1000])
set(gca,'XTick', [1 2])
title('RT')
set(gca,'XTicklabels', {'MT','MIP'})

% print -depsc descr




    disp('MT vs MIP - Accuracy - NonTMS')
    [h p ci stats]=ttest(ACC_MT,ACC_MIP);
    t=stats.tstat;
    Mean1= mean(ACC_MT);
    Mean2= mean(ACC_MIP);
    SD1 = std(ACC_MT);
    SD2 = std(ACC_MIP);
    D= (Mean2-Mean1)/(sqrt(((SD1)^2 +(SD2)^2)/2));
    fprintf('Acc: t(%i) = %2.2f, p = %1.3f, D = %2.2f \n', stats.df, t,p,D)

    disp('MT vs MIP - RT - NonTMS')
    [h p ci stats]=ttest(RT_MT,RT_MIP);
    t=stats.tstat;
    Mean1= mean(RT_MT);
    Mean2= mean(RT_MIP);
    SD1 = std(RT_MT);
    SD2 = std(RT_MIP);
    D= (Mean2-Mean1)/(sqrt(((SD1)^2 +(SD2)^2)/2));
    fprintf('RT: t(%i) = %2.2f, p = %1.3f, D = %2.2f \n', stats.df, t,p,D)

    disp('MT vs MIP - Accuracy - TMS')
    [h p ci stats]=ttest(ACC_MT1,ACC_MIP1);
    t=stats.tstat;
    Mean1= mean(ACC_MT1);
    Mean2= mean(ACC_MIP1);
    SD1 = std(ACC_MT1);
    SD2 = std(ACC_MIP1);
    D= (Mean2-Mean1)/(sqrt(((SD1)^2 +(SD2)^2)/2));
    fprintf('Acc: t(%i) = %2.2f, p = %1.3f, D = %2.2f \n', stats.df, t,p,D)

    disp('MT vs MIP - RT - TMS')
    [h p ci stats]=ttest(RT_MT1,RT_MIP1);
    t=stats.tstat;
    Mean1= mean(RT_MT1);
    Mean2= mean(RT_MIP1);
    SD1 = std(RT_MT1);
    SD2 = std(RT_MIP1);
    D= (Mean2-Mean1)/(sqrt(((SD1)^2 +(SD2)^2)/2));
    fprintf('RT: t(%i) = %2.2f, p = %1.3f, D = %2.2f \n', stats.df, t,p,D)
    
    

    disp('MIP - RT - TMS vs Non')
    [h p ci stats]=ttest(RT_MIP,RT_MIP1);
    t=stats.tstat;
    Mean1= mean(RT_MIP);
    Mean2= mean(RT_MIP1);
    SD1 = std(RT_MIP);
    SD2 = std(RT_MIP1);
    D= (Mean2-Mean1)/(sqrt(((SD1)^2 +(SD2)^2)/2));
    fprintf('RT: t(%i) = %2.2f, p = %1.3f, D = %2.2f \n', stats.df, t,p,D)
    
    disp('MT - RT - TMS vs Non')
    [h p ci stats]=ttest(RT_MT,RT_MT1);
    t=stats.tstat;
    Mean1= mean(RT_MT);
    Mean2= mean(RT_MT);
    SD1 = std(RT_MT);
    SD2 = std(RT_MT1);
    D= (Mean2-Mean1)/(sqrt(((SD1)^2 +(SD2)^2)/2));
    fprintf('RT: t(%i) = %2.2f, p = %1.3f, D = %2.2f \n', stats.df, t,p,D)

    

    
titles={'MIP0','MIP1','MT0','MT1'};
for i = 1:2
    if i==1
        fprintf('Accuracy\n')
        tab=array2table([ACC_MIP', ACC_MIP1', ACC_MT', ACC_MT1'],'variablenames',titles);
    else
        fprintf('RT\n')
        tab=array2table([RT_MIP', RT_MIP1', RT_MT', RT_MT1'],'variablenames',titles);
    end
    %set up ANOVA
    within1 = categorical([1 1 0 0 ])';
    within2 = categorical([0 1 0 1  ])';
    within = table(within1,within2,'variablenames',{'Session','Stim'});
    factors='Session*Stim';

    allvars=strcat(titles{1},'-',titles{end},' ~1');
    rm = fitrm(tab,allvars,'WithinDesign',within);
    ranovatbl = ranova(rm,'withinmodel',factors);

    ranovatbl = table2array(ranovatbl);
    p_spots=[3 5 7];
     p_spot_names={'Session' 'Stim' 'Session*Stim'};

    sig_idx=[p_spots];
    for i=1:length(sig_idx)
       fprintf('%s : F(%d,%d) = %4.2f, p = %1.3f, n = %4.3f \n',p_spot_names{find(p_spots==sig_idx(i))},  ranovatbl(1,2),ranovatbl(2,2), ranovatbl ((sig_idx(i)),4),ranovatbl((sig_idx(i)),5), (ranovatbl((sig_idx(i)),1)/(ranovatbl((sig_idx(i)),1) + ranovatbl((sig_idx(i))+1,1))))
    end   
end
    
    

%follow up1
    disp('MIP - Acc - TMS vs Non')
    [h p ci stats]=ttest(ACC_MIP,ACC_MIP1);
    t=stats.tstat;
    Mean1= mean(ACC_MIP);
    Mean2= mean(ACC_MIP1);
    SD1 = std(ACC_MIP);
    SD2 = std(ACC_MIP1);
    D= (Mean2-Mean1)/(sqrt(((SD1)^2 +(SD2)^2)/2));
    fprintf('RT: t(%i) = %2.2f, p = %1.3f, D = %2.2f \n', stats.df, t,p,D)
    
    disp('MT - Acc - TMS vs Non')
    [h p ci stats]=ttest(ACC_MT,ACC_MT1);
    t=stats.tstat;
    Mean1= mean(ACC_MT);
    Mean2= mean(ACC_MT);
    SD1 = std(ACC_MT);
    SD2 = std(ACC_MT1);
    D= (Mean2-Mean1)/(sqrt(((SD1)^2 +(SD2)^2)/2));
    fprintf('RT: t(%i) = %2.2f, p = %1.3f, D = %2.2f \n', stats.df, t,p,D)

    
    %%folow up 2
    
     disp('MIP - RT - TMS vs Non')
    [h p ci stats]=ttest(ACC_MIP,ACC_MIP1);
    t=stats.tstat;
    Mean1= mean(ACC_MIP);
    Mean2= mean(ACC_MIP1);
    SD1 = std(ACC_MIP);
    SD2 = std(ACC_MIP1);
    D= (Mean2-Mean1)/(sqrt(((SD1)^2 +(SD2)^2)/2));
    fprintf('RT: t(%i) = %2.2f, p = %1.3f, D = %2.2f \n', stats.df, t,p,D)
    
    disp('MT - Acc - TMS vs Non')
    [h p ci stats]=ttest(ACC_MT,ACC_MT1);
    t=stats.tstat;
    Mean1= mean(ACC_MT);
    Mean2= mean(ACC_MT);
    SD1 = std(ACC_MT);
    SD2 = std(ACC_MT1);
    D= (Mean2-Mean1)/(sqrt(((SD1)^2 +(SD2)^2)/2));
    fprintf('RT: t(%i) = %2.2f, p = %1.3f, D = %2.2f \n', stats.df, t,p,D)



% cd('J:\PolyU\TMS\Paper')
% 
% print -depsc basic_beh






