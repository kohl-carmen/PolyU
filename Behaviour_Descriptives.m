%% Descriptives: Behaviour
% based on: BasicBeh_for_Publication.m 

% - reports descriptive stats (Accuracy and RT) for non-tms behaviour (MT vs MIP)
% - reports difference between MIP and MT non-tms behaviour (Accuracy and RT)
% - plots Figure 2 a/b

clearvars 
beh_data_folder='D:\PolyU\TMS\Data\Beh';
plot_dir = 'C:\Users\ckohl\Desktop';

%% load and sort data
Partic=1:31;
rt_min=0;
rt_max=1500;

ACCURACY = [];
SESSION = [];
TMS = [];
RT = [];
ACC_MT = [];
ACC_MIP = [];
RT_MT = [];
RT_MIP = [];
for partic=1:length(Partic)
    partic_char = num2str(Partic(partic));
    if Partic(partic)<10
        partic_char=strcat('0',num2str(Partic(partic)));
    end

    MIP = load(strcat(beh_data_folder,'\7',partic_char,'1\transformed\data_behavior_tms'));
    V5 = load(strcat(beh_data_folder,'\7',partic_char,'0\transformed\data_behavior_tms'));     
   
    for file={'MIP', 'V5'}
        eval([file{1} '=' file{1} '.data.behavior;'])
    end
   
    %     session         1=MIP   0=V5
    %     trial           trial number
    %     rt              RT in ms
    %     tms             0=no 1=stim
    %     tms_cond        0=no 1=MIP 2=V5
    %     accuracy        1=corr 0=err nan=empty/distractor

    % variables of interest
    session=[ones(length(MIP.trial{:}),1); zeros(length(V5.trial{:}),1)];%1=MIP, 0=V5;
    tms=[MIP.tms{:};V5.tms{:}];  
    rt=[MIP.RT{:};V5.RT{:}];
    accuracy=[MIP.accuracy{:};V5.accuracy{:}]; 
   
    % collect per partic
    ACCURACY(:,partic)=(accuracy);
    SESSION(:,partic)=session;
    TMS(:,partic)=tms;
    RT(:,partic)=rt;

    % exclude trials
    rt_keep=rt;
    accuracy_deletion=(isnan(accuracy));%delete all trials where distractor or empty was chosen
    DELETE(:,partic)=(accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max);
    
    ACC_MT(partic)=mean(ACCURACY(SESSION(:,partic)==0 & DELETE(:,partic)==0 & TMS(:,partic)==0,partic));
    ACC_MIP(partic)=mean(ACCURACY(SESSION(:,partic)==1 & DELETE(:,partic)==0 & TMS(:,partic)==0,partic));
   
    RT_MT(partic)=mean(RT(SESSION(:,partic)==0 & DELETE(:,partic)==0 & TMS(:,partic)==0,partic));
    RT_MIP(partic)=mean(RT(SESSION(:,partic)==1 & DELETE(:,partic)==0 & TMS(:,partic)==0,partic));
end


%% Report descriptives
fprintf('\nDescriptives:\n')
fprintf('Non-tms accuracy (across sessions): %2.2f (SD=%2.2f) \n',mean([ACC_MT,ACC_MIP]).*100, std([ACC_MT,ACC_MIP]).*100)
fprintf('Non-tms RT (across sessions): %2.2f (SD=%2.2f) \n',mean([RT_MT,RT_MIP]), std([RT_MT,RT_MIP]))

%% Report t-tests between MT and MIP (non-tms)
fprintf('\nDifferences between sessions in non-tms trials:\n')
[binary,p,ci,tstat]=ttest(ACC_MT,ACC_MIP);
Mean1= mean(ACC_MT);
Mean2= mean(ACC_MIP);
SD1 = std(ACC_MT);
SD2 = std(ACC_MIP);
D= (Mean2-Mean1)/(sqrt(((SD1)^2 +(SD2)^2)/2));
fprintf('t(%i) = %2.2f, p = %1.3f, d = %2.2f\n', tstat.df, tstat.tstat, p,D)
[binary,p,ci,tstat]=ttest(RT_MT,RT_MIP);
Mean1= mean(RT_MT);
Mean2= mean(RT_MIP);
SD1 = std(RT_MT);
SD2 = std(RT_MIP);
D= (Mean2-Mean1)/(sqrt(((SD1)^2 +(SD2)^2)/2));
fprintf('t(%i) = %2.2f, p = %1.3f, d = %2.2f\n', tstat.df, tstat.tstat, p,D)



%% Plot Figure 2 a, b
x=(ACC_MT);
SEM = std(x)/sqrt(length(x));               % Standard Error
error_ACC_MT=[SEM];

x=(ACC_MIP);
SEM = std(x)/sqrt(length(x));               % Standard Error
error_ACC_MIP=[SEM];

x=(RT_MT);
SEM = std(x)/sqrt(length(x));               % Standard Error
error_RT_MT=[SEM];

x=(RT_MIP);
SEM = std(x)/sqrt(length(x));               % Standard Error
error_RT_MIP=[SEM];

clf
subplot(1,2,1)
hold on
bar([mean(ACC_MT), mean(ACC_MIP)],'r')
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

%% save figure
cd(plot_dir)
print -depsc basic_beh



