%% Behavioural Analysis for TMS7 Publication
%== 'Beh_for_Publication.m' on drive
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
ANOVA=2; % which ANOVA 1:Stim*Sess  2:Stim*Sess*Dlocation 
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
 
% We always want to normalise the individual terms of an interaction (and
% then everything below at normalise(regressors)


Regs={'HV+LV' 'HV-LV' 'D-HV' };
model='[hv_v(idx)+lv_v(idx), hv_v(idx)-lv_v(idx),d_v(idx)-hv_v(idx)]';




if stepwise
    Regs={'HV+LV' 'HV-LV' };
    model='[hv_v(idx)+lv_v(idx), hv_v(idx)-lv_v(idx)]';
    Regs2={'D'};
    model2='[d_v(idx)]';
end


%exclude bad tms
if exclude_bad_tms
    load('D:\PolyU\TMS\Data\Navigator\TMS_DEV')
    count_tms_exclusions=0;
end

if magprob_weight
    load('D:\PolyU\TMS\Data\Beh\magprob_betas_diff')
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
            if Conds{run}(1:2)=='MI'
                this_session=1;
                tms_index=4;
            else
                this_session=0;
                tms_index=3;
            end
            
            if Conds{run}(tms_index)=='1'
                this_tms=1;
            else
                this_tms=0;
            end
            
            if Conds{run}(end-2:end)=='psi'
                this_loc=0; %ipsi
            else
                this_loc=1; %contra
            end
             
            pos_oi=3;
           

            if Testing==1
                warning('') 
                idx=(session==this_session & tms==this_tms );
                if length(Conds)>4 %considers distractor location
                    idx=(session==this_session & tms==this_tms  & (position(:,pos_oi)==locations{this_loc+1}(1) | position(:,pos_oi)==locations{this_loc+1}(2) ));
                end
                
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

            elseif Testing==2
                idx=(session==this_session & tms==this_tms & accuracy_used==1);

                
                if length(Conds)>4 %considers distractor location
                    idx=( accuracy_used==1 & session==this_session & tms==this_tms  & (position(:,pos_oi)==locations{loc_oi+1}(1) | position(:,pos_oi)==locations{loc_oi+1}(2) ));
                end


                regressors=eval(model);

                regressors=normalise(regressors);

                criterion=([rt(idx)]);% Y is a vector of response values.  If DISTR is 'binomial' Y may a binary vector indicating success/failure, and the total number of trials is taken to be 1 for all observations.  If DISTR is 'binomial', Y may also be a two column matrix, the first column containing the number of successes for each observation,and the second containing the total number of trials.

                [betas(partic,:,run)]=glmfit(regressors,criterion);
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
%% PUBLICATION FIGURE

%Conds: 'MIP0ipsi'    'MIP1ipsi'    'MT0ipsi'    'MT1ipsi'    'MIP0contra'    'MIP1contra'    'MT0contra'    'MT1contra'
close all
% for rfig=2: length(Regs)+1 %one figure per regressor (except intercept)
%     
%     figure(rfig-1)
%     to_plot=[];
%     
%     %% if split into MIP/MT instead of contra/ipsi (Bolton said contra/ipsi but MIP/MT looks more impressive)
%     to_plot=[];
%     error_to_plot=[];
%     for MTMIP=1:2
%         if MTMIP==1
%             these_conds=[3 7 4 8];
%         else
%             these_conds=[1 5 2 6];
%         end
%         for cond=1:length(Conds)/2% ipsi0 contra0 ipsi1 contra1
%             to_plot(MTMIP,cond)=mean(betas(:,rfig,these_conds(cond)))
%             error_to_plot(MTMIP,cond)=([errors.(strcat('Reg',num2str(rfig-1))).(Conds{these_conds(cond)})]);
%     
%         end
%     end
%     h = bar(to_plot,.9);
%     colours={ [.5 .75 1]   [ 0 .438 .875] [ 1 .75 .5]  [.875 .438 0]};
%     for i=1:length(Conds)/2
%          h(i).FaceColor=colours{i};
%     end
%     hold on
%   
%     [numgroups,numbars ] = size(to_plot); 
%     groupwidth = min(0.8, numbars/(numbars+1.5));
%     for i = 1:numbars
%       % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange: https://uk.mathworks.com/matlabcentral/answers/102220-how-do-i-place-errorbars-on-my-grouped-bar-graph-using-function-errorbar-in-matlab-7-13-r2011b
%       x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
%       errorbar(x, to_plot(:,i), error_to_plot(:,i), 'k', 'linestyle', 'none');
%     end
% 
%     if rfig==2 | rfig==4
%         temp1=sprintf('Non-TMS \n D: ipsilateral');
%         temp2=sprintf('Non-TMS \n D: contralateral');
%         temp3=sprintf('TMS \n D: ipsilateral');
%         temp4=sprintf('TMS \n D: contralateral');
% %         legend({'Non-TMS  D-ipsilateral'  'Non-TMS  D-contralateral' 'TMS  D-ipsilateral' 'TMS D-contralateral'})
%          legend(temp1, temp2, temp3, temp4)
%     end
%     
%     if rfig <4
%         ylim([-.6 1.2])
%         set(gca,'YTick',[-.5 0 .5 1 1.5])
%     else
%         set(gca,'YTick',[-.2 0 .2 .4])
%     end
%     set(gca,'xticklabel',{'MT' 'MIP'})
%     ticks=get(gca,'YTick');
%     title(Regs{rfig-1})
%     
%     
% 
%        if rfig==2
%            print -depsc HV+LV_bar
%        elseif rfig==3
%            print -depsc HV-LV_bar
%        elseif rfig==4
%            print -depsc D-HV_bar
%        end
% 
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try line plot

%Conds: 'MIP0ipsi'    'MIP1ipsi'    'MT0ipsi'    'MT1ipsi'    'MIP0contra'    'MIP1contra'    'MT0contra'    'MT1contra'
close all
for rfig=2: length(Regs)+1 %one figure per regressor (except intercept)
    
    figure(rfig-1)
    to_plot=[];

     %plot interaction
    clf
    %1 MIP0ips 2 MIP1ips 3 MT0ips 4 MT1ips
    %5 MIP0con 6 MIP1con 7 MT1con 8 MT0con
    hold on
    xlim([0 5])
%     ylim([-.05 .05])
     
    for lines=1:4
        if lines==1 %MT0ipsi MT0con
            these_conds=[3 8]; 
            these_x=[1 2];
            col=[.5 .5 .5];
            linetype='--';
        elseif lines==2 %MT1ipsi MT1con
            these_conds=[4 7];
            these_x=[1 2];
            col=['k'];
            linetype='-';
        else
            these_x=[3 4];
            if lines==3 %MIp0ipsi MIP0con
                these_conds=[1 5];
                col=[ 1 .75 .5]; 
                linetype='--';
            else
                these_conds=[2 6];
                col=[.875 .438 0];
                linetype='-';
            end
        end
        to_plot=squeeze(mean(betas(:,rfig,these_conds)));
        error_to_plot=[([errors.(strcat('Reg',num2str(rfig-1))).(Conds{these_conds(1)})]); ([errors.(strcat('Reg',num2str(rfig-1))).(Conds{these_conds(2)})])];
        
        plot([0:5], zeros(1,6),':', 'Color','k');%[.7 .7 .7])
    
        leg(lines)=plot(these_x,[to_plot],linetype, 'Color', col, 'Linewidth',2);
        errorbar(these_x(1), to_plot(1), error_to_plot(1), 'Color', col, 'Linewidth',2);
        errorbar(these_x(2), to_plot(2), error_to_plot(2), 'Color', col, 'Linewidth',2);

        
    end
        
      if rfig==2 | rfig==4
        temp1=sprintf('MT Non-TMS');
        temp2=sprintf('MT TMS');
        temp3=sprintf('MIP Non-TMS');
        temp4=sprintf('MIP TMS');
%         legend({'Non-TMS  D-ipsilateral'  'Non-TMS  D-contralateral' 'TMS  D-ipsilateral' 'TMS D-contralateral'})
        legend(leg,temp1, temp2, temp3, temp4)
      end
    
    if rfig <4
      ylim([-.6 1.2])
        set(gca,'YTick',[-.5 0 .5 1 1.5])
    else
        ylim([-.3 .6])
        set(gca,'YTick',[-.3 0 .3 .6])
    end
    set(gca,'XTick',[1:4])
    set(gca,'xticklabel',{'Contra' 'Ipsi' 'Contra' 'Ipsi'})
    ticks=get(gca,'YTick');
    title(Regs{rfig-1})
    
    
       if rfig==2
           print -depsc HV+LV_line
       elseif rfig==3
           print -depsc HV-LV_line
       elseif rfig==4
           print -depsc D-HV_line
       end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANOVA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for regressor=2:length(Regs)+1

        titles=Conds;%{'MIP0','MIP1','MT0','MT1'};

        tab=array2table(squeeze(betas(:,regressor,:)),'variablenames',titles);
        if ANOVA==1
            within1 = categorical([1 1 0 0])';
            within2 = categorical([0 1 0 1])';
            within = table(within1,within2,'variablenames',{'Session','Stim'});
            factors='Session*Stim';
        elseif ANOVA==2
            within1 = categorical([1 1 0 0 1 1 0 0 ])';
            within2 = categorical([0 1 0 1 0 1 0 1])';
            within3=  categorical([0 0 0 0 1 1 1 1])';
            within = table(within1,within2,within3,'variablenames',{'Session','Stim','DLoc'});
            factors='Session*Stim*DLoc';
        end
        allvars=strcat(titles{1},'-',titles{end},' ~1');
        rm = fitrm(tab,allvars,'WithinDesign',within);
        ranovatbl = ranova(rm,'withinmodel',factors)
        if ANOVA==2
            p_spots=[3 5 7 9 11 13 15];
            p_spot_names={'Session' 'Stim' 'DLoc' 'Session*Stim' 'Session*DLoc*' 'Stim*DLoc' 'Session*Stim*DLoc'};
        elseif ANOVA==1
            p_spots=[3 5 7];
            p_spot_names={'Session' 'Stim' 'Session*Stim*'};
        end
        fprintf('\n- - %s - -\n',Regs{regressor-1})
        ranova(rm,'withinmodel',factors);
        ranovatbl=table2array(ranovatbl);
        if sum(ranovatbl([p_spots],5)<.1)>0
           sig_idx=find(ranovatbl([p_spots],5)<.1);
           sig_idx=p_spots(sig_idx);
           for i=1:length(sig_idx)
               fprintf('%s : F(%d,%d) = %4.2f, p = %1.3f \n',p_spot_names{find(p_spots==sig_idx(i))},  ranovatbl(1,2),ranovatbl(2,2), ranovatbl ((sig_idx(i)),4),ranovatbl((sig_idx(i)),5))
           end   
        else
            disp('nothing <.1')
        end
    end
    fprintf('\n For posthoc: SPSS \n')
                    


%% test each regressor against 0
mean_betas=squeeze(mean(betas,3));
p=[];
t=[];
D=[];
fprintf('Ttest of mean betas against 0:\n')
for i=2:size(mean_betas,2)
    [h p(i) ci stats]=ttest(mean_betas(:,i));
    t(i)=stats.tstat;
    p(i);
    Mean1= mean(mean_betas(:,i));
    Mean2= 0;
    SD1 = std(mean_betas(:,i));
    SD2 = 0;
%     D(i)= (Mean2-Mean1)/(sqrt(((SD1)^2 +(SD2)^2)/2))
    D(i)= (Mean2-Mean1)/(SD1)';   
    fprintf(Regs{i-1})
    fprintf('\nt(%i) = %2.3f, p = %2.3f, d = %2.3f\n',stats.df, t(i), p(i), D(i))
end
                       
                        
                        
%% Explore Interaction



clf
%%-------------------------------------------------------------------------
%% ANOVA:       2
%% CRITERION:   ACC
%% REGRESSORS:  (HV-LV)*D_E
%%-------------------------------------------------------------------------
ANOVA=2;
Comparison={ 'D_V'};  
bin_x=[0:.125:0.5];
avg_window=0.5;
% bin_x=[0 .5];
% bin_y=[0 .5];
% avg_window=0.5;
% 
% bin_x=[0:1/3:2/3];
% avg_window=1/3;
%  bin_x=[bin_x;bin_x+avg_window];
%  
%  bin_x=[0:1/4:3/4];
% avg_window=1/4;
 bin_x=[bin_x;bin_x+avg_window];
 
if ANOVA==2
    position_oi=3; %ipsi/contra loc refers to distractor
elseif ANOVA==3
    position_oi=1; %ipsi/contra loc refers to HV
end

% tic
%   h = actxserver('PowerPoint.Application');
%     Presentation = h.Presentation.Add;

    Conds={'MIP0ipsi' 'MIP1ipsi' 'MIP0contra' 'MIP1contra'};
Runs=length(Conds);
ACC_PER=[];
ERR_PER=[];
TRIALS_PER=[];
DV_PER=[];
keep_acc_per=[];
 for run=1:Runs
        if any(run==[1 2 5 6]) % MIP
            this_session=1;
        elseif any(run==[3 4 7 8]) % MT
            this_session=0;
        end
        if any(run==[1 3 5 7]) % no stim
            this_tms=0;%no stim
        elseif any(run==[2 4 6 8]) %  stim
            this_tms=1;
        end
        if any(run==[1 2 3 4]) % d ipsilateral
            d_loc=[1 3];
        elseif any(run==[5 6 7 8]) %  d contralateral
            d_loc=[2 4];
        end

        one_dimension=0;
        acc_per=[];
        x_idx=[];
        y_idx=[];
        rt_per=[];
        trials_per_rt=[];
        trials_per_acc=[];
         X=eval(Comparison{1});
  
        if length(Comparison)>1
            Y=eval(Comparison{2});
        else 
            Y=X;
            one_dimension=1;
        end
        for bin_x_count=1:length(bin_x) 
                for partic=1:length(Partic)
                      %pick trials according to cond etc
                      if ANOVA==1
                          idx=(DELETE(:,partic)==0 & SESSION(:,partic)==this_session & TMS(:,partic)==this_tms & ~isnan(ACCURACY_BOLTON(:,partic)));
                      else
                        idx=(DELETE(:,partic)==0 & (POSITION(:,position_oi,partic)==d_loc(1) | POSITION(:,position_oi,partic)==d_loc(2)) & SESSION(:,partic)==this_session & TMS(:,partic)==this_tms & ~isnan(ACCURACY_BOLTON(:,partic)));
                      end
                      
                      
                      %% quantiles based on trials
%                       accuracy=ACCURACY_BOLTON(idx,partic);
%                       
%                       bin_idx= round(quantile(1:length((X(idx,partic))),bin_x(1,bin_x_count))):round(quantile(1:length((X(idx,partic))),bin_x(2,bin_x_count)));
%                       
%                       [sorted_X, sort_i]=sort(X(idx,partic));
%                       
%                       dv_per(bin_x_count,partic)=mean(sorted_X);
%                       
%                       accuracy= accuracy(sort_i);
%                       acc_per(bin_x_count,partic)=mean(accuracy(bin_idx));
%                       trials_per(bin_x_count,partic)=length(accuracy(bin_idx));
%                       
                      %% quantiles based on D value
%                       %select bin
                      idx= idx & X(:,partic)>= quantile(X(idx,partic),bin_x(1,bin_x_count)) &  X(:,partic)<= quantile(X(idx,partic),bin_x(2,bin_x_count));     
                                               
                      acc_per(bin_x_count,partic)=mean(ACCURACY_BOLTON(idx,partic));
                            
                      trials_per(bin_x_count,partic)=sum(sum(idx,partic));
                      
                      dv_per(bin_x_count,partic)=mean(D_V(idx,partic));
                     
%                      % this is just to label x/y axes appropriately later
%                      x_idx(bin_x_count)=mean(quantile(X,bin_x(2,bin_x_count)));
                         
                end
                ACC_PER(bin_x_count)=mean(acc_per(bin_x_count,:),2);
                DV_PER(bin_x_count)=mean(dv_per(bin_x_count,:),2);
%                 TRIALS_PER(bin_x_count)=mean(trials_per(bin_x_count,:),2);
                errorx=acc_per(bin_x_count,:);
                SEM = std(errorx)/sqrt(length(errorx));               % Standard Error
                ts = tinv([0.025  0.975],length(errorx)-1);      % T-Score
                CI = mean(errorx) + ts*SEM;                      % Confidence Intervals
                ERR_PER(bin_x_count) =[SEM];
               
                
        end
        keep_acc_per(:,:,run)=acc_per;
        
%         acc_per=nanmean(acc_per,3);
%         rt_per=nanmean(rt_per,3);
%         trials_per_rt=mean(trials_per_rt,3);
%         trials_per_acc=mean(trials_per_acc,3);
             
        figure(1)
%         subplot(2,length(Conds)/2,run)
        hold on
        if Conds{run}(2)=='I' & Conds{run}(4)=='1' & Conds{run}(end)=='a'
            linewidth=4;
            linestyle='-';
            color=[223/255 112/255 14/255];
        elseif Conds{run}(2)=='I' & Conds{run}(4)=='1' & Conds{run}(end)=='i'
            linewidth=1;
            color='k';
            linestyle='-';
        elseif Conds{run}(2)=='I' & Conds{run}(4)=='0' & Conds{run}(end)=='i'
            linewidth=1;
            color=[128/255 128/255 128/255];
            linestyle='--';
        else
            linewidth=1;
            color=[250/255 190/255 129/255];
            linestyle='--';
        end
            
%         plot(bin_x(2,:),ACC_PER,'Linewidth',linewidth, 'Color', color)
%         errorbar([bin_x(2,:)], (ACC_PER), ERR_PER, 'linestyle', 'none','Linewidth',linewidth, 'Color', color);
        plot(DV_PER,ACC_PER,'Linewidth',linewidth, 'Color', color,'Linestyle',linestyle)
        errorbar(DV_PER, (ACC_PER), ERR_PER, 'linestyle', 'none','Linewidth',linewidth, 'Color', color,'Linestyle',linestyle);
        
%         plot(bin_x(2,:),TRIALS_PER/100,'Linewidth',1, 'Color', [.5 .5 .5])
     
        title(strcat('Acc -', Conds{run}))
%         Conds_title=Conds([1 2 5 6]);
%         title(strcat('Acc -', Conds_title{run}))
        xlabel(Comparison{1})
        if one_dimension==0
            ylabel(Comparison{2})
        end
        

        set(gca,'YTick',[.68 .7 .72 .74]);
        set(gca,'XTick',[100 200 300]);


 end
 
cd('D:\PolyU\TMS\Paper')

print -depsc intplot


           
%% NOTE: the one just aboveplots D not D-HV!!
%% If you want D-HV: this is the best I can do
% 
% clf
% ANOVA=2;
% Comparison={ 'D_V-HV_V'};  
% bin_x=[0:.125:0.5];
% avg_window=0.5;
% bin_x=[0:.1:0.6];
% avg_window=0.4;
% bin_x=[0 .5];
% bin_y=[0 .5];
% avg_window=0.5;
% % 
% bin_x=[0:1/3:2/3];
% avg_window=1/3;
% %  bin_x=[bin_x;bin_x+avg_window];
% %  
% %  bin_x=[0:1/4:3/4];
% % avg_window=1/4;
%  bin_x=[bin_x;bin_x+avg_window];
%  
% if ANOVA==2
%     position_oi=3; %ipsi/contra loc refers to distractor
% elseif ANOVA==3
%     position_oi=1; %ipsi/contra loc refers to HV
% end
% 
% % tic
% %   h = actxserver('PowerPoint.Application');
% %     Presentation = h.Presentation.Add;
% 
%     Conds={'MIP0ipsi' 'MIP1ipsi' 'MIP0contra' 'MIP1contra'};
% Runs=length(Conds);
% ACC_PER=[];
% ERR_PER=[];
% TRIALS_PER=[];
% DV_PER=[];
% keep_acc_per=[];
%  for run=1:Runs
%         if any(run==[1 2 5 6]) % MIP
%             this_session=1;
%         elseif any(run==[3 4 7 8]) % MT
%             this_session=0;
%         end
%         if any(run==[1 3 5 7]) % no stim
%             this_tms=0;%no stim
%         elseif any(run==[2 4 6 8]) %  stim
%             this_tms=1;
%         end
%         if any(run==[1 2 3 4]) % d ipsilateral
%             d_loc=[1 3];
%         elseif any(run==[5 6 7 8]) %  d contralateral
%             d_loc=[2 4];
%         end
% 
%         one_dimension=0;
%         acc_per=[];
%         x_idx=[];
%         y_idx=[];
%         rt_per=[];
%         trials_per_rt=[];
%         trials_per_acc=[];
%          X=eval(Comparison{1});
%   
%         if length(Comparison)>1
%             Y=eval(Comparison{2});
%         else 
%             Y=X;
%             one_dimension=1;
%         end
%         for bin_x_count=1:length(bin_x) 
%                 for partic=1:length(Partic)
%                       %pick trials according to cond etc
%                       if ANOVA==1
%                           idx=(DELETE(:,partic)==0 & SESSION(:,partic)==this_session & TMS(:,partic)==this_tms & ~isnan(ACCURACY_BOLTON(:,partic)));
%                       else
%                         idx=(DELETE(:,partic)==0 & (POSITION(:,position_oi,partic)==d_loc(1) | POSITION(:,position_oi,partic)==d_loc(2)) & SESSION(:,partic)==this_session & TMS(:,partic)==this_tms & ~isnan(ACCURACY_BOLTON(:,partic)));
%                       end
%                       
%                       
%                       %% quantiles based on trials
%                       accuracy=ACCURACY_BOLTON(idx,partic);
%                       
%                       bin_idx= round(quantile(1:length((X(idx,partic))),bin_x(1,bin_x_count))):round(quantile(1:length((X(idx,partic))),bin_x(2,bin_x_count)));
%                       
%                       [sorted_X, sort_i]=sort(X(idx,partic));
%                       
%                       dv_per(bin_x_count,partic)=mean(sorted_X);
%                       
%                       accuracy= accuracy(sort_i);
%                       acc_per(bin_x_count,partic)=mean(accuracy(bin_idx));
%                       trials_per(bin_x_count,partic)=length(accuracy(bin_idx));
%                       
%                       % quantiles based on D value
% %                       %select bin
% %                       idx= idx & X(:,partic)>= quantile(X(idx,partic),bin_x(1,bin_x_count)) &  X(:,partic)<= quantile(X(idx,partic),bin_x(2,bin_x_count));     
% %                                                
% %                       acc_per(bin_x_count,partic)=mean(ACCURACY_BOLTON(idx,partic));
% %                             
% %                       trials_per(bin_x_count,partic)=sum(sum(idx,partic));
% %                       
% %                       dv_per(bin_x_count,partic)=mean(D_V(idx,partic));
% %                      
% % %                      % this is just to label x/y axes appropriately later
% %                      x_idx(bin_x_count)=mean(quantile(X,bin_x(2,bin_x_count)));
%                          
%                 end
%                 ACC_PER(bin_x_count)=mean(acc_per(bin_x_count,:),2);
%                 DV_PER(bin_x_count)=mean(dv_per(bin_x_count,:),2);
% %                 TRIALS_PER(bin_x_count)=mean(trials_per(bin_x_count,:),2);
%                 errorx=acc_per(bin_x_count,:);
%                 SEM = std(errorx)/sqrt(length(errorx));               % Standard Error
%                 ts = tinv([0.025  0.975],length(errorx)-1);      % T-Score
%                 CI = mean(errorx) + ts*SEM;                      % Confidence Intervals
%                 ERR_PER(bin_x_count) =[SEM];
%                
%                 
%         end
%         keep_acc_per(:,:,run)=acc_per;
%         
% %         acc_per=nanmean(acc_per,3);
% %         rt_per=nanmean(rt_per,3);
% %         trials_per_rt=mean(trials_per_rt,3);
% %         trials_per_acc=mean(trials_per_acc,3);
%              
%         figure(1)
% %         subplot(2,length(Conds)/2,run)
%         hold on
%         if Conds{run}(2)=='I' & Conds{run}(4)=='1' & Conds{run}(end)=='a'
%             linewidth=4;
%             linestyle='-';
%             color=[223/255 112/255 14/255];
%         elseif Conds{run}(2)=='I' & Conds{run}(4)=='1' & Conds{run}(end)=='i'
%             linewidth=1;
%             color='k';
%             linestyle='-';
%         elseif Conds{run}(2)=='I' & Conds{run}(4)=='0' & Conds{run}(end)=='i'
%             linewidth=1;
%             color=[128/255 128/255 128/255];
%             linestyle='--';
%         else
%             linewidth=1;
%             color=[250/255 190/255 129/255];
%             linestyle='--';
%         end
%             
%         plot(bin_x(2,:),ACC_PER,'Linewidth',linewidth, 'Color', color)
%         errorbar([bin_x(2,:)], (ACC_PER), ERR_PER, 'linestyle', 'none','Linewidth',linewidth, 'Color', color);
% %         plot(DV_PER,ACC_PER,'Linewidth',linewidth, 'Color', color,'Linestyle',linestyle)
% %         errorbar(DV_PER, (ACC_PER), ERR_PER, 'linestyle', 'none','Linewidth',linewidth, 'Color', color,'Linestyle',linestyle);
%         
% %         plot(bin_x(2,:),TRIALS_PER/100,'Linewidth',1, 'Color', [.5 .5 .5])
%      
%         title(strcat('Acc -', Conds{run}))
% %         Conds_title=Conds([1 2 5 6]);
% %         title(strcat('Acc -', Conds_title{run}))
%         xlabel(Comparison{1})
%         if one_dimension==0
%             ylabel(Comparison{2})
%         end
%         
% % 
% %         set(gca,'YTick',[.68 .7 .72 .74]);
% %         set(gca,'XTick',[100 200 300]);
% 
% 
%  end
% cd('J:\PolyU\TMS\Paper')
% 
% print -depsc d-HV2f
