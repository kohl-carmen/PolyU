%% this is GLM1 for Behaviour, but collpased over all non-tms trials
%% we just run the GLM and test the betas against 0
%% this is now in the paper, due to Bolton's changes
%% see paper for definition of GLM1

%% Behavioural Analysis for TMS7 Publication
% Makes Figure 2.
% copied from, 'BehGLM_forsimplestacceffect.m' (just a bit cut down and
% with different figure part)

% GLMs: Predict Accuracy using the regrresors:
% HV+LV
% HV-LV
% D-HV

% run that GLM per person, per session(MIP/MT) and by stim cond(tms/non)
% Then run Session(MIP/MT) x (TMS(stim/non) ANOVA on each regressor

% NOTE: Posthoc tests are done in SPSS. include covariates.

% first part very similar to Beh.m
clearvars 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% --- DEFINE --- %%
% Testing effects on ACC(1) or RT (2)?
Testing=1; %1=ACC; 2=RT

% Criterion: Accuracy relative (1) or absolute (2)? Or (3): porbability of choosing each option. As defiend by Gluth et al (2018).(relative == bolton accuracy; abolsute==mine)
Criterion_Type=1; 
    %% NOTE:(1) or (2)makes no difference when Testing==2;
    %% NOTE:(3) does not work for all options yet (doesn't work for RT).
    %% NOTE:(3) only works for very few cases (we often don't have people choosing all options in all conditions)
    %% NOTE:(3) doesn't currently have a working ANOVA attached to it (I guess just put a loop around it though)

% Which ANOVA: Stim*Sess (1) or Stim*Sess*Dlocation (2)?
ANOVA=1; % which ANOVA 1:Stim*Sess  2:Stim*Sess*Dlocation 
    ANOVA2_val=1; % (1) D location, (2) HV location? if ANOVA==2, the 3rd factoris contra vs ipsi. ANOVA2_val defines which values we're looking at.

% RT cutoff defined per person (1) or not (0)?
define_RT_cutoff_per_partic=0	; % per participant or not
    
% Exclusions based on task performance (1) or not (0)? Define cut off
task_exclusion=0;
    accuracy_cutoff=.7;
    %% NOTE: if task_exclusion is on, match_exclusion does not matter

% Print to ppt (1) or not (0)?    
ppt=0;
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
%% MODELS

if ANOVA==1
    Conds={'MIP0' 'MIP1' 'MT0' 'MT1'};
elseif ANOVA==2 
    Conds={'MIP0ipsi' 'MIP1ipsi' 'MT0ipsi' 'MT1ipsi' 'MIP0contra' 'MIP1contra' 'MT0contra' 'MT1contra'};
end
Conds={'all_nontms'};
 
% We always want to normalise the individual terms of an interaction (and
% then everything below at normalise(regressors)


Regs={'HV+LV' 'HV-LV' 'D-HV' };
model='[hv_v(idx)+lv_v(idx), hv_v(idx)-lv_v(idx),d_v(idx)-hv_v(idx)]';

%now:
Regs={'HV+LV' 'HV-LV' 'D-HV' '(HV-LV)*(D-HV)'};
model='[hv_v(idx)+lv_v(idx), hv_v(idx)-lv_v(idx),d_v(idx)-hv_v(idx),normalise(hv_v(idx)-lv_v(idx)).*normalise(d_v(idx)-hv_v(idx))]';
   

if stepwise
    Regs={'HV+LV' 'HV-LV' };
    model='[hv_v(idx)+lv_v(idx), hv_v(idx)-lv_v(idx)]';
    Regs2={'D'};
    model2='[d_v(idx)]';
end


%exclude bad tms
if exclude_bad_tms
    load('J:\PolyU\TMS\Data\Navigator\TMS_DEV')
    count_tms_exclusions=0;
end

if magprob_weight
    load('J:\PolyU\TMS\Data\Beh\magprob_betas_diff')
end

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
    ACCURACY(partic)=nanmean(accuracy_bolton);
    
    if magprob_weight
        weight=magprob_betas_diff(partic,2:end);
        weight=[1 1 1 1];
        hv_v=(all_mag(:,1).*weight(1))+(all_prob(:,1).*weight(2));
        lv_v=(all_mag(:,2).*weight(1))+(all_prob(:,2).*weight(2));
        d_v=(all_mag(:,3).*weight(1))+(all_prob(:,3).*weight(2));
    end
    
    
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
       
   ACCURACY_BOLTON(:,partic )=accuracy_bolton;
   ACC(partic,:)=[partic,mean(accuracy_bolton(~isnan(accuracy_bolton)))];
   
       ACCURACY=[];
    ACCURACY(:,partic)=(accuracy_bolton);
    HV_V(:,partic)=(hv_v);
    D_V(:,partic)=d_v;
    LV_V(:,partic)=lv_v;
    POSITION(:,:,partic)=position;
    SESSION(:,partic)=session;
    TMS(:,partic)=tms;

    
     %% EXCLUDE TRIALS 
    
    if define_RT_cutoff_per_partic
        rt_max=median(rt);
    end
    rt_keep=rt;
    
    %% ?!?!?!?!!? %%
    if Criterion_Type==1
        accuracy_deletion=(isnan(accuracy_bolton));%delete all trials where distractor or empty was chosen
    elseif Criterion_Type>1
        accuracy_deletion=zeros(length(accuracy_bolton),1);%don't
    end
    
    if exclude_bad_tms
%         tms_dev=[TMS_DEV(:,partic,2);TMS_DEV(:,partic,1)];
        %these have practice stims in them, so delete those
        tms_dev=[TMS_DEV(11:end,Partic(partic),2);TMS_DEV(11:end,Partic(partic),1)];
        if sum(tms)~=length(tms_dev)
            error
        end
        tms_dev_temp=zeros(size(tms));
        tms_dev_temp(tms~=0)=tms_dev;
        tms_dev=tms_dev_temp;
        tms_deletion=tms_dev>tms_dev_cutoff;% & session==1;
        count_tms_exclusions=count_tms_exclusions+sum(tms_deletion);
    else
        tms_deletion=zeros(length(accuracy_bolton),1);
    end
     
   
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


    if Criterion_Type==1
        accuracy_used=accuracy_bolton;
    elseif Criterion_Type>1
        accuracy_used=accuracy_abs;
    end
    
    clear acc data file MIP V5 temp
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SETUP ANALYSIS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
        %four models per person
        for run=1:length(Conds)
                warning('') 
                idx=( tms==0 );
                
                if sum(idx)==0
                    exclude=[exclude,partic];
                else
            
                regressors=[eval(model)];

                regressors=[normalise(regressors)];

                criterion=[accuracy_used(idx)];% Y is a vector of response values.  If DISTR is 'binomial' Y may a binary vector indicating success/failure, and the total number of trials is taken to be 1 for all observations.  If DISTR is 'binomial', Y may also be a two column matrix, the first column containing the number of successes for each observation,and the second containing the total number of trials.
                if magprob_weight
                    criterion(hv_v(idx)<lv_v(idx) & accuracy_used(idx)==0)=1;
                    criterion(hv_v(idx)<lv_v(idx) & accuracy_used(idx)==1)=0;
                end
                        
                [betas(partic,:,run),dev,stats]=glmfit(regressors,criterion,'binomial');
                if stepwise
                    regressors=[eval(model2)];

                    regressors=[normalise(regressors)];

                    criterion=stats.resid;
                    [betas2(partic,:,run),dev,stats]=glmfit(regressors,criterion);
                end
                [warnMsg, warnId] = lastwarn;
                if ~isempty(warnMsg)
                    exclude=[exclude,partic];
                end
                end

            
        end

end
exclude=unique(exclude);
fprintf('%i participants excluded \n',length(exclude))

discarded_acc=0;
all=0;
for partic=1:length(Partic)
    discarded_acc=discarded_acc+sum(isnan(ACCURACY_BOLTON(:,partic)));
    all=all+length(ACCURACY_BOLTON(:,partic));
end
fprintf(' %2.2f %% of trials were discarded because they picked nonsense \n',discarded_acc./all.*100)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT BETAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf

if stepwise==1 
    betas(:,end+1,:)=betas2(:,end,:);
    Regs={Regs{:},Regs2{:}};
end

betas(exclude,:,:)=[];

for c=1:length(Conds)
    for r=1:length(Regs)+1
        x=(betas(~isnan(betas(:,r ,c)),r ,c));
        SEM = std(x)/sqrt(length(x));               % Standard Error
        ts = tinv([0.025  0.975],length(x)-1);      % T-Score
        CI = mean(x) + ts*SEM;                      % Confidence Intervals
        errors.(strcat('Reg',num2str(r))).(Conds{c})=[CI(2)-mean(x)];
        errors.(strcat('Reg',num2str(r))).(Conds{c})=[SEM];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  FIGURE (if you want)
plott=0;
if plott
    %Conds: 'MIP0ipsi'    'MIP1ipsi'    'MT0ipsi'    'MT1ipsi'    'MIP0contra'    'MIP1contra'    'MT0contra'    'MT1contra'
    close all
    for rfig=2: length(Regs)+1 %one figure per regressor (except intercept)

        figure(rfig-1)
        to_plot=[];

        %% if split into MIP/MT instead of contra/ipsi (Bolton said contra/ipsi but MIP/MT looks more impressive)
        to_plot=[];
        error_to_plot=[];
        for MTMIP=1:2
            if MTMIP==1
                these_conds=[3  4 ];
            else
                these_conds=[1  2 ];
            end
            for cond=1:length(Conds)/2% ipsi0 contra0 ipsi1 contra1
                to_plot(MTMIP,cond)=mean(betas(:,rfig,these_conds(cond)));
                error_to_plot(MTMIP,cond)=([errors.(strcat('Reg',num2str(rfig-1))).(Conds{these_conds(cond)})]);

            end
        end
        h = bar(to_plot,.9);
        colours={ [.5 .75 1]   [ 0 .438 .875] [ 1 .75 .5]  [.875 .438 0]};
        for i=1:length(Conds)/2
             h(i).FaceColor=colours{i};
        end
        hold on

        [numgroups,numbars ] = size(to_plot); 
        groupwidth = min(0.8, numbars/(numbars+1.5));
        for i = 1:numbars
          % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange: https://uk.mathworks.com/matlabcentral/answers/102220-how-do-i-place-errorbars-on-my-grouped-bar-graph-using-function-errorbar-in-matlab-7-13-r2011b
          x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
          errorbar(x, to_plot(:,i), error_to_plot(:,i), 'k', 'linestyle', 'none');
        end

    %     if rfig==2 | rfig==4
    %         temp1=sprintf('Non-TMS \n D: ipsilateral');
    %         temp2=sprintf('Non-TMS \n D: contralateral');
    %         temp3=sprintf('TMS \n D: ipsilateral');
    %         temp4=sprintf('TMS \n D: contralateral');
    % %         legend({'Non-TMS  D-ipsilateral'  'Non-TMS  D-contralateral' 'TMS  D-ipsilateral' 'TMS D-contralateral'})
    %          legend(temp1, temp2, temp3, temp4)
    %     end

        if rfig <4
            ylim([-.6 1.2])
            set(gca,'YTick',[-.5 0 .5 1 1.5])
        else
            set(gca,'YTick',[-.2 0 .2 .4])
        end
        set(gca,'xticklabel',{'MT' 'MIP'})
        ticks=get(gca,'YTick');
        title(Regs{rfig-1})



           if rfig==2
               print -depsc HV+LV_bar
           elseif rfig==3
               print -depsc HV-LV_bar
           elseif rfig==4
               print -depsc D-HV_bar
           end

    end
end

%% NEW T TEST
%% test each regressor against 0
%only non-tms
for regressor=2:size(betas,2)

    [h p ci stats]=ttest(betas(:,regressor));
    fprintf('\nTesting %s against 0:\n',Regs{regressor-1})
    fprintf('- all nontms: t(%d) = %2.2f, p = %2.3f\n',stats.df,stats.tstat, p)
    
end

%% same but with effect size
for regressor=2:size(betas,2)
    [h p ci stats]=ttest(betas(:,regressor));
    D=(mean(betas(:,regressor))-0)/std(betas(:,regressor));
    fprintf('\nTesting %s against 0:\n',Regs{regressor-1})
    fprintf('- all nontms: t(%d) = %2.2f, p = %2.3f, d = %2.2f\n',stats.df,stats.tstat, p,D)
    
end



%bar chart
% GLM1 on non-tms trials
%get errors:
for r=1:length(Regs)+1
    x=(betas(~isnan(betas(:,r ,c)),r ,c));
    SEM = std(x)/sqrt(length(x));               % Standard Error
    ts = tinv([0.025  0.975],length(x)-1);      % T-Score
    CI = mean(x) + ts*SEM;                      % Confidence Intervals
    errors(r)=[CI(2)-mean(x)];
    errors(r)=[SEM];
end
%plot
clf
a=bar(mean(betas(:,2:end)));
hold on
err=errorbar(1:4,mean(betas(:,2:end)),mean(betas(:,2:end))-errors(2:end),mean(betas(:,2:end))+errors(2:end));
err.Color = [0 0 0];                            
err.LineStyle = 'none';  
set(gca,'xticklabel',Regs)

print -depsc GLM1bar