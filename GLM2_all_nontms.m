%% Basic Behavioural Analysis
% Runs GLM2 on median split (HV-LV) non-tms data (across conditions) to get base behaviour
%       GLM2:	Step 1, β0 + β1 z(HV-LV) + β2 z(HV+LV) + ε1
%               Step 2, β3 + β4 z(D-HV) + ε2

% Reports ttests of each predictor's beta weights against zero
% Plots Figure 2d (each D-HV beta weight)

clearvars 

beh_data_folder='D:\PolyU\TMS\Data\Beh';
plot_dir = 'J:\PolyU\TMS\Paper';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Partic=1:31;
rt_min=0;
rt_max=1500;

Conds={'all_nontms'};

Partic=1:31;
stepwise=1;
rt_min=0;%100;
rt_max=1500;%1499;


Conds={'all_nontms'};
% We always want to normalise the individual terms of an interaction (and
% then everything below at normalise(regressors)


%now:
Regs={'HV+LV' 'HV-LV' 'D-HV' '(HV-LV)*(D-HV)'};
model='[hv_v(idx)+lv_v(idx), hv_v(idx)-lv_v(idx),d_v(idx)-hv_v(idx),normalise(hv_v(idx)-lv_v(idx)).*normalise(d_v(idx)-hv_v(idx))]';
   

if stepwise
    Regs={'HV+LV' 'HV-LV' };
    model='[hv_v(idx)+lv_v(idx), hv_v(idx)-lv_v(idx)]';
    Regs2={'D-HV'};
    model2='[d_v(idx)-hv_v(idx)]';
end


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

     for file={'MIP', 'V5'};
         eval([file{1} '=' file{1} '.data.behavior;'])
     end
    
     
    %combine and make a data matrix that I'm more comfi with
    %     session         1=MIP   0=V5
    %     trial           trial number
    %     rt              RT in ms
    %     tms             0=no 1=stim
    %     tms_cond        0=no 1=MIP 2=V5
    %     accuracy	      1=corr 0=err nan=empty/distractor
    %     d_v             distractor value(prob*rew)
    %     lv_v            low value option value(prob*rew)
    %     hv_v            high value option value(prob*rew)

    %% caps ones for all partics

    session=[ones(length(MIP.trial{:}),1); zeros(length(V5.trial{:}),1)];%1=MIP, 0=V5;
    trial=[MIP.trial{:};V5.trial{:}];  
    tms=[MIP.tms{:};V5.tms{:}];  
    rt=[MIP.RT{:};V5.RT{:}];
    accuracy=[MIP.accuracy{:};V5.accuracy{:}]; 
    d_v=[MIP.D{:};V5.D{:}]; 
    lv_v=[MIP.LV{:};V5.LV{:}]; 
    hv_v=[MIP.HV{:};V5.HV{:}];

    
     %% EXCLUDE TRIALS 
    
    rt_keep=rt;
    
    accuracy_deletion=(isnan(accuracy));%delete all trials where distractor or empty was chosen
    
    tms_deletion=zeros(length(accuracy),1);

   
    session(accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max |tms_deletion)=[];        
    trial(accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
    rt(accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
    tms( accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
    accuracy( accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];                
    d_v( accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
    lv_v( accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
    hv_v( accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
           
    %% median split(HV-LV)
    mediansplit=zeros(size(hv_v));
    mediansplit(hv_v-lv_v>=median(hv_v-lv_v))=1;
    
    clear acc data file MIP V5 temp
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SETUP ANALYSIS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for median_split=[0 1]

        idx=( tms==0 & mediansplit==median_split);

        regressors=[eval(model)];
        regressors=[normalise(regressors)];
        criterion=[accuracy(idx)];% Y is a vector of response values.  If DISTR is 'binomial' Y may a binary vector indicating success/failure, and the total number of trials is taken to be 1 for all observations.  If DISTR is 'binomial', Y may also be a two column matrix, the first column containing the number of successes for each observation,and the second containing the total number of trials.

        if median_split==0
            [betas_low(partic,:),dev,stats]=glmfit(regressors,criterion,'binomial');
        else
            [betas_high(partic,:),dev,stats]=glmfit(regressors,criterion,'binomial');
        end
        if stepwise
            regressors=[eval(model2)];
            regressors=[normalise(regressors)];
            criterion=stats.resid;

            if median_split==0
                [betas2_low(partic,:),dev,stats]=glmfit(regressors,criterion);
            else
                [betas2_high(partic,:),dev,stats]=glmfit(regressors,criterion);
            end
        end
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT BETAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if stepwise==1 
    betas_low(:,end+1,:)=betas2_low(:,end,:);
    betas_high(:,end+1,:)=betas2_high(:,end,:);
    Regs={Regs{:},Regs2{:}};
end


%% test beta weights (D-HV) against zero (once for easy/hard)
regressor=4;
% easy
fprintf('\nD-HV: EASY')
betas=betas_high(:,:);
[h p ci stats]=ttest(betas(:,regressor));
D=(mean(betas(:,regressor))-0)/std(betas(:,regressor));
fprintf('\nTesting %s against 0:\n',Regs{regressor-1})
fprintf('- all_nontms: t(%d) = %2.2f, p = %2.3f, d = %2.2f\n',stats.df,stats.tstat, p,D)


% easy
fprintf('\nD-HV: HARD')
betas=betas_low(:,:);
[h p ci stats]=ttest(betas(:,regressor));
D=(mean(betas(:,regressor))-0)/std(betas(:,regressor));
fprintf('\nTesting %s against 0:\n',Regs{regressor-1})
fprintf('- all-nontms: t(%d) = %2.2f, p = %2.3f, d = %2.2f\n',stats.df,stats.tstat, p,D)


%% Figure 2 d
%get errors:
regressor=4;
x=betas_high(:,regressor);
SEM = std(x)/sqrt(length(x));% Standard Error
errors_high=[SEM];

x=betas_low(:,regressor);
SEM = std(x)/sqrt(length(x));     
errors_low=[SEM];

%plot
figure
bar([mean(betas_low(:,regressor));mean(betas_high(:,regressor))]);
hold on
err=errorbar(1:2,[mean(betas_low(:,regressor)),mean(betas_high(:,regressor))], [mean(betas_low(:,regressor)),mean(betas_high(:,regressor))]-[errors_low,errors_high],[mean(betas_low(:,regressor)),mean(betas_high(:,regressor))]+[errors_low,errors_high]);
err.Color = [0 0 0];                            
err.LineStyle = 'none';  
set(gca,'xticklabel',{'hard','easy'})

%% save
cd(plot_dir)
rint -depsc GLM2bar







