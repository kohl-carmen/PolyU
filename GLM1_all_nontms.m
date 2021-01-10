%% Basic Behavioural Analysis
% Runs GLM1 on non-tms data (across conditions) to get base behaviour
%       GLM1:	β0 + β1 z(HV-LV) + β2 z(HV+LV) + β3 z(D-HV) + β4 z(HV-LV) z(D-HV) + ε

% Reports ttests of each predictor's beta weights against zero
% Plots Figure 2d (each beta weight)

clearvars 

beh_data_folder='D:\PolyU\TMS\Data\Beh';
plot_dir = 'C:\Users\ckohl\Desktop';
output_dir='C:\Users\ckohl\Desktop\Current\Other\Bolton\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Partic=1:31;
rt_min=0;
rt_max=1500;

Conds={'all_nontms'};
 

% GLM1
Regs={'HV+LV' 'HV-LV' 'D-HV' '(HV-LV)*(D-HV)'};
model='[hv_v(idx)+lv_v(idx), hv_v(idx)-lv_v(idx),d_v(idx)-hv_v(idx),normalise(hv_v(idx)-lv_v(idx)).*normalise(d_v(idx)-hv_v(idx))]';
 % We always want to normalise the individual terms of an interaction (and
% then everything below at normalise(regressors)

for partic=1:length(Partic)
    partic_char = sprintf('%02d', Partic(partic));
    MIP=load(strcat(beh_data_folder,'\7',partic_char,'1\transformed\data_behavior_tms'));
    V5=load(strcat(beh_data_folder,'\7',partic_char,'0\transformed\data_behavior_tms'));      
    for file={'MIP', 'V5'}
        eval([file{1} '=' file{1} '.data.behavior;'])
    end
    
    %     session         1=MIP   0=V5
    %     trial           trial number
    %     rt              RT in ms
    %     tms             0=no 1=stim
    %     tms_cond        0=no 1=MIP 2=V5
    %     accuracy        1=corr 0=err nan=empty/distractor
    %     d_v             distractor value(prob*rew)
    %     lv_v            low value option value(prob*rew)
    %     hv_v            high value option value(prob*rew)
    
    % select variables of interest
    session=[ones(length(MIP.trial{:}),1); zeros(length(V5.trial{:}),1)];%1=MIP, 0=V5;
    tms=[MIP.tms{:};V5.tms{:}];  
    rt=[MIP.RT{:};V5.RT{:}];
    accuracy=[MIP.accuracy{:};V5.accuracy{:}]; 
    d_v=[MIP.D{:};V5.D{:}]; 
    lv_v=[MIP.LV{:};V5.LV{:}]; 
    hv_v=[MIP.HV{:};V5.HV{:}];
   
    % exclude 
    rt_keep=rt;
    accuracy_deletion=(isnan(accuracy));%delete all trials where distractor or empty was chosen
    tms_deletion=zeros(length(accuracy),1);  
   
    session(accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max |tms_deletion)=[];        
    rt(accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
    tms( accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
    accuracy( accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];                
    d_v( accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
    lv_v( accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
    hv_v( accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
    
    clear acc data file MIP V5 temp
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% GLM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        idx=( tms==0 );

        regressors=[eval(model)];
        regressors=[normalise(regressors)];
        criterion=[accuracy(idx)];% Y is a vector of response values.  If DISTR is 'binomial' Y may a binary vector indicating success/failure, and the total number of trials is taken to be 1 for all observations.  If DISTR is 'binomial', Y may also be a two column matrix, the first column containing the number of successes for each observation,and the second containing the total number of trials.

        [betas(partic,:),dev,stats]=glmfit(regressors,criterion,'binomial');
      
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% test each regressor against 0
%only non-tms
for regressor=2:size(betas,2)
    [h p ci stats]=ttest(betas(:,regressor));
    D=(mean(betas(:,regressor))-0)/std(betas(:,regressor));
    fprintf('\nTesting %s against 0:\n',Regs{regressor-1})
    fprintf('- all nontms: t(%d) = %2.2f, p = %2.3f, d = %2.2f\n',stats.df,stats.tstat, p,D)
    
end


%% Plot Figure 2 c
% GLM1 on non-tms trials
errors=[];
for r=1:length(Regs)+1
    x=(betas(~isnan(betas(:,r )),r ));
    SEM = std(x)/sqrt(length(x));% Standard Error
    errors(r)=[SEM];
end
%plot
figure
bar(mean(betas(:,2:end)));
hold on
err=errorbar(1:4,mean(betas(:,2:end)),errors(2:end));
err.Color = [0 0 0];                            
err.LineStyle = 'none';  
set(gca,'xticklabel',Regs)

%% save
cd(plot_dir)
print -depsc GLM1bar
