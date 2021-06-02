%% Behavioural Analysis fo Supplement S3 (RT and Acc)
% 	- applies GLM3 ({'HV+LV' 'HV-LV' 'D-HV' }) to behavioural data (tms x
% 	session x side) once for RT and once for Acc
% 	- reports GLM effects (tms x session x side x AccvsRT ANOVA) 
% 	- plots supplementary figure S3



clearvars 

beh_data_folder='D:\PolyU\TMS\Data\Beh';
plot_dir = 'C:\Users\ckohl\Desktop';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Partic=1:31;
rt_min=0;%100;
rt_max=1500;%1499;
exclude=[];

%% GLM
Conds={'MIP0ipsi' 'MIP1ipsi' 'MT0ipsi' 'MT1ipsi' 'MIP0contra' 'MIP1contra' 'MT0contra' 'MT1contra'};

Regs={'HV+LV' 'HV-LV' 'D-HV' };
model='[hv_v(idx)+lv_v(idx), hv_v(idx)-lv_v(idx),d_v(idx)-hv_v(idx)]';

% We always want to normalise the individual terms of an interaction (and
% then everything below at normalise(regressors)

ACCURACY = [];

for partic=1:length(Partic)
    partic_char = sprintf('%02d', Partic(partic));
    MIP=load(strcat(beh_data_folder,'\7',partic_char,'1\transformed\data_behavior_tms'));
    V5=load(strcat(beh_data_folder,'\7',partic_char,'0\transformed\data_behavior_tms'));     

    for file={'MIP', 'V5'}
        eval([file{1} '=' file{1} '.data.behavior;'])
    end
      
    %combine and make a data matrix that I'm more comfi with
    %     session         1=MIP   0=V5
    %     trial           trial number
    %     rt              RT in ms
    %     tms             0=no 1=stim
    %     tms_cond        0=no 1=MIP 2=V5
    %     accuracy        1=corr 0=err nan=empty/distractor
    %     d_v             distractor value(prob*rew)
    %     lv_v            low value option value(prob*rew)
    %     hv_v            high value option value(prob*rew)
    %     position        positions of stims on screen - 1st column: HV, 2nd column: LV, 3rd column: D - 1=top left, 2=top right, 3=bottom left, 4=bottom right
    
    %% caps ones for all partics

    session=[ones(length(MIP.trial{:}),1); zeros(length(V5.trial{:}),1)];%1=MIP, 0=V5;
    trial=[MIP.trial{:};V5.trial{:}];  
    tms=[MIP.tms{:};V5.tms{:}];  
    rt=[MIP.RT{:};V5.RT{:}];
    accuracy=[MIP.accuracy{:};V5.accuracy{:}]; 
    d_v=[MIP.D{:};V5.D{:}]; 
    lv_v=[MIP.LV{:};V5.LV{:}]; 
    hv_v=[MIP.HV{:};V5.HV{:}];
    position=[MIP.pos_rews{:};V5.pos_rews{:}];   
    
    %for distractor location
    locations={[1 3], [2 4]}; %ipsi=[1 3]=left
    di_loc=zeros(size(hv_v));
    di_loc(position(:,3)==2 |position(:,3)==4)=1;
    dloc_for_glm=di_loc;    
   
    %collect partics
    ACCURACY(:,partic)=(accuracy);
   
    % exclude   
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
    position(accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion,:)=[];   
    dloc_for_glm(accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion,:)=[];       
         
    clear acc data file MIP V5 temp
    
     rt = rt/1000;
    %     for i = 1:length(rt)
    %     rt(i)=1/rt(i);
    %     end
    %     rt = (rt-min(rt))/(max(rt)-min(rt));
    %% GLM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

        warning('') 
        idx=(session==this_session & tms==this_tms  & (position(:,pos_oi)==locations{this_loc+1}(1) | position(:,pos_oi)==locations{this_loc+1}(2) ));

        if sum(idx)==0
            exclude=[exclude,partic];
        else
            regressors=[eval(model)];
            regressors=[normalise(regressors)];
            criterion=[accuracy(idx)];% Y is a vector of response values.  If DISTR is 'binomial' Y may a binary vector indicating success/failure, and the total number of trials is taken to be 1 for all observations.  If DISTR is 'binomial', Y may also be a two column matrix, the first column containing the number of successes for each observation,and the second containing the total number of trials.

            [betas_acc(partic,:,run),dev,stats]=glmfit(regressors,criterion,'binomial');
            


            idx=(session==this_session & accuracy == 1 & tms==this_tms  & (position(:,pos_oi)==locations{this_loc+1}(1) | position(:,pos_oi)==locations{this_loc+1}(2) ));
            regressors=[eval(model)];
            regressors=[normalise(regressors)];criterion=[rt(idx)]; 
            
            [betas_rt(partic,:,run),dev,stats]=glmfit(regressors,criterion);

            [warnMsg, warnId] = lastwarn;
            if ~isempty(warnMsg)
                exclude=[exclude,partic];
            end
        end

    end

end
% exclude
exclude=unique(exclude);
fprintf('%i participants excluded \n',length(exclude))

discarded_acc=0;
all=0;
for partic=1:length(Partic)
    discarded_acc=discarded_acc+sum(isnan(ACCURACY(:,partic)));
    all=all+length(ACCURACY(:,partic));
end
fprintf(' %2.2f %% of trials were discarded because they picked nonsense \n',discarded_acc./all.*100)

betas_acc(exclude,:,:)=[];
betas_rt(exclude,:,:)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANOVA ON RT AND ACC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for regressor=4%:length(Regs)+1
    titles ={};
    for i = 1:length(Conds)
        titles{i} = strcat('Acc_',Conds{i});
        titles{i+length(Conds)} = strcat('RT_',Conds{i});
    end
    tab=array2table([squeeze(betas_acc(:,regressor,:)),squeeze(betas_rt(:,regressor,:))],'variablenames',titles);
    
    %set up ANOVA
    within1 = categorical([1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0])';
    within2 = categorical([0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1])';
    within3=  categorical([0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1])';
    within4 = categorical([0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1])';
    within = table(within1,within2,within3,within4,'variablenames',{'Session','Stim','DLoc','AccRT'});
    factors='Session*Stim*DLoc*AccRT';

    allvars=strcat(titles{1},'-',titles{end},' ~1');
    rm = fitrm(tab,allvars,'WithinDesign',within);
    ranovatbl = ranova(rm,'withinmodel',factors)
    
%     %report
%     p_spots=[3 5 7 9 11 13 15];
%     p_spot_names={'Session' 'Stim' 'DLoc' 'Session*Stim' 'Session*DLoc*' 'Stim*DLoc' 'Session*Stim*DLoc'};
% 
%     fprintf('\n- - %s - -\n',Regs{regressor-1})
%     ranova(rm,'withinmodel',factors);
%     ranovatbl=table2array(ranovatbl);
%     %%this is to only report significant ones
%     % if sum(ranovatbl([p_spots],5)<.1)>0
%     %    sig_idx=find(ranovatbl([p_spots],5)<.1);
%     %    sig_idx=p_spots(sig_idx);
%     %    for i=1:length(sig_idx)
%     %        fprintf('%s : F(%d,%d) = %4.2f, p = %1.3f, n = %4.3f \n',p_spot_names{find(p_spots==sig_idx(i))},  ranovatbl(1,2),ranovatbl(2,2), ranovatbl ((sig_idx(i)),4),ranovatbl((sig_idx(i)),5), (ranovatbl((sig_idx(i)),1)/(ranovatbl((sig_idx(i)),1) + ranovatbl((sig_idx(i))+1,1))))
%     %    end   
%     % else
%     %     disp('nothing <.1')
%     % end
%     sig_idx=[p_spots];
%     for i=1:length(sig_idx)
%        fprintf('%s : F(%d,%d) = %4.2f, p = %1.3f, n = %4.3f \n',p_spot_names{find(p_spots==sig_idx(i))},  ranovatbl(1,2),ranovatbl(2,2), ranovatbl ((sig_idx(i)),4),ranovatbl((sig_idx(i)),5), (ranovatbl((sig_idx(i)),1)/(ranovatbl((sig_idx(i)),1) + ranovatbl((sig_idx(i))+1,1))))
%     end   

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT BETAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 
% errors
    figure
    subplot(1,2,1)
    hold on
betas = betas_acc;

for c=1:length(Conds)
    for r=1:length(Regs)+1
        x=(betas(~isnan(betas(:,r ,c)),r ,c));
        SEM = std(x)/sqrt(length(x));               % Standard Error
        errors.(strcat('Reg',num2str(r))).(Conds{c})=[SEM];
    end
end

rfig=4%: length(Regs)+1 %one figure per regressor (except intercept)   

    to_plot=[];


    %1 MIP0ips 2 MIP1ips 3 MT0ips 4 MT1ips
    %5 MIP0con 6 MIP1con 7 MT0con 8 MT1con
    hold on
    xlim([0 5])
     
    for lines=1:4
        if lines==1 %MT0ipsi MT0con
            these_conds=[3 7]; 
            these_x=[1 2];
            col=[.5 .5 .5];
            linetype='--';
        elseif lines==2 %MT1ipsi MT1con
            these_conds=[4 8];
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
        error_to_plot=[([errors.(strcat('Reg',num2str(rfig))).(Conds{these_conds(1)})]); ([errors.(strcat('Reg',num2str(rfig))).(Conds{these_conds(2)})])];
        
        plot([0:5], zeros(1,6),':', 'Color','k');%[.7 .7 .7])
    
        leg(lines)=plot(these_x,[to_plot],linetype, 'Color', col, 'Linewidth',2);
        errorbar(these_x(1), to_plot(1), error_to_plot(1), 'Color', col, 'Linewidth',2);
        errorbar(these_x(2), to_plot(2), error_to_plot(2), 'Color', col, 'Linewidth',2);      
    end
        
    temp={'MT Non-TMS','MT TMS','MIP Non-TMS','MIP TMS'};
    legend(leg,temp)
   
    if rfig <4
        ylim([-.6 1.2])
        set(gca,'YTick',[-.5 0 .5 1 1.5])
    else
        ylim([-.3 .6])
        set(gca,'YTick',[-.3 0 .3 .6])
    end
    set(gca,'XTick',[1:4])
    set(gca,'xticklabel',{'Ipsi' 'Contra' 'Ipsi' 'Contra'})
    ticks=get(gca,'YTick');
    title(strcat('Accuracy: ',Regs{rfig-1}))
   
    
betas = betas_rt;

for c=1:length(Conds)
    for r=1:length(Regs)+1
        x=(betas(~isnan(betas(:,r ,c)),r ,c));
        SEM = std(x)/sqrt(length(x));               % Standard Error
        errors.(strcat('Reg',num2str(r))).(Conds{c})=[SEM];
    end
end

rfig=4%: length(Regs)+1 %one figure per regressor (except intercept)   
subplot(1,2,2)
hold on
    to_plot=[];

    %1 MIP0ips 2 MIP1ips 3 MT0ips 4 MT1ips
    %5 MIP0con 6 MIP1con 7 MT0con 8 MT1con
    hold on
    xlim([0 5])
     
    for lines=1:4
        if lines==1 %MT0ipsi MT0con
            these_conds=[3 7]; 
            these_x=[1 2];
            col=[.5 .5 .5];
            linetype='--';
        elseif lines==2 %MT1ipsi MT1con
            these_conds=[4 8];
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
        error_to_plot=[([errors.(strcat('Reg',num2str(rfig))).(Conds{these_conds(1)})]); ([errors.(strcat('Reg',num2str(rfig))).(Conds{these_conds(2)})])];
        
        plot([0:5], zeros(1,6),':', 'Color','k');%[.7 .7 .7])
    
        leg(lines)=plot(these_x,[to_plot],linetype, 'Color', col, 'Linewidth',2);
        errorbar(these_x(1), to_plot(1), error_to_plot(1), 'Color', col, 'Linewidth',2);
        errorbar(these_x(2), to_plot(2), error_to_plot(2), 'Color', col, 'Linewidth',2);      
    end
         set(gca,'XTick',[1:4])
    set(gca,'xticklabel',{'Ipsi' 'Contra' 'Ipsi' 'Contra'})   
    temp={'MT Non-TMS','MT TMS','MIP Non-TMS','MIP TMS'};
%     legend(leg,temp)
   
    title(strcat('RT: ',Regs{rfig-1}))
    
   
% cd('C:\Users\ckohl\Desktop\')
% print -depsc HV+LV_line



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANOVA JUST ON RT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


titles=Conds;
tab=array2table(squeeze(betas_rt(:,regressor,:)),'variablenames',titles);

%set up ANOVA
within1 = categorical([1 1 0 0 1 1 0 0 ])';
within2 = categorical([0 1 0 1 0 1 0 1])';
within3=  categorical([0 0 0 0 1 1 1 1])';
within = table(within1,within2,within3,'variablenames',{'Session','Stim','DLoc'});
factors='Session*Stim*DLoc';

allvars=strcat(titles{1},'-',titles{end},' ~1');
rm = fitrm(tab,allvars,'WithinDesign',within);
ranovatbl = ranova(rm,'withinmodel',factors);

%report
p_spots=[3 5 7 9 11 13 15];
p_spot_names={'Session' 'Stim' 'DLoc' 'Session*Stim' 'Session*DLoc*' 'Stim*DLoc' 'Session*Stim*DLoc'};

fprintf('\n- - %s - -\n',Regs{regressor-1})
ranova(rm,'withinmodel',factors);
ranovatbl=table2array(ranovatbl);

sig_idx=[p_spots];
for i=1:length(sig_idx)
fprintf('%s : F(%d,%d) = %4.2f, p = %1.3f, n = %4.3f \n',p_spot_names{find(p_spots==sig_idx(i))},  ranovatbl(1,2),ranovatbl(2,2), ranovatbl ((sig_idx(i)),4),ranovatbl((sig_idx(i)),5), (ranovatbl((sig_idx(i)),1)/(ranovatbl((sig_idx(i)),1) + ranovatbl((sig_idx(i))+1,1))))
end   
    
    
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANOVA JUST ON ACC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
titles=Conds;
tab=array2table(squeeze(betas_acc(:,regressor,:)),'variablenames',titles);

%set up ANOVA
within1 = categorical([1 1 0 0 1 1 0 0 ])';
within2 = categorical([0 1 0 1 0 1 0 1])';
within3=  categorical([0 0 0 0 1 1 1 1])';
within = table(within1,within2,within3,'variablenames',{'Session','Stim','DLoc'});
factors='Session*Stim*DLoc';

allvars=strcat(titles{1},'-',titles{end},' ~1');
rm = fitrm(tab,allvars,'WithinDesign',within);
ranovatbl = ranova(rm,'withinmodel',factors);

%report
p_spots=[3 5 7 9 11 13 15];
p_spot_names={'Session' 'Stim' 'DLoc' 'Session*Stim' 'Session*DLoc*' 'Stim*DLoc' 'Session*Stim*DLoc'};

fprintf('\n- - %s - -\n',Regs{regressor-1})
ranova(rm,'withinmodel',factors);
ranovatbl=table2array(ranovatbl);

sig_idx=[p_spots];
for i=1:length(sig_idx)
fprintf('%s : F(%d,%d) = %4.2f, p = %1.3f, n = %4.3f \n',p_spot_names{find(p_spots==sig_idx(i))},  ranovatbl(1,2),ranovatbl(2,2), ranovatbl ((sig_idx(i)),4),ranovatbl((sig_idx(i)),5), (ranovatbl((sig_idx(i)),1)/(ranovatbl((sig_idx(i)),1) + ranovatbl((sig_idx(i))+1,1))))
end   




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GITHUB VERSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - same code but runs with github data

% %% Behavioural Analysis fo Supplement S3 (RT and Acc)
% % 	- applies GLM3 ({'HV+LV' 'HV-LV' 'D-HV' }) to behavioural data (tms x
% % 	session x side) once for RT and once for Acc
% % 	- reports GLM effects (tms x session x side x AccvsRT ANOVA) 
% % 	- plots supplementary figure S3
% 
% 
% 
% clearvars 
% 
% beh_data_folder='D:\PolyU\TMS\Data\Beh';
% plot_dir = 'C:\Users\ckohl\Desktop';
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% Partic=1:31;
% rt_min=0;%100;
% rt_max=1500;%1499;
% exclude=[];
% 
% %% GLM
% Conds={'MIP0ipsi' 'MIP1ipsi' 'MT0ipsi' 'MT1ipsi' 'MIP0contra' 'MIP1contra' 'MT0contra' 'MT1contra'};
% 
% Regs={'HV+LV' 'HV-LV' 'D-HV' };
% model='[hv_v(idx)+lv_v(idx), hv_v(idx)-lv_v(idx),d_v(idx)-hv_v(idx)]';
% 
% % We always want to normalise the individual terms of an interaction (and
% % then everything below at normalise(regressors)
% 
% ACCURACY = [];
% 
% for partic=1:length(Partic)
%     partic_char = sprintf('%02d', Partic(partic));
%     MIP=load(strcat(beh_data_folder,'\7',partic_char,'1\transformed\data_behavior_tms'));
%     V5=load(strcat(beh_data_folder,'\7',partic_char,'0\transformed\data_behavior_tms'));     
% 
%     for file={'MIP', 'V5'}
%         eval([file{1} '=' file{1} '.data.behavior;'])
%     end
%       
%     %combine and make a data matrix that I'm more comfi with
%     %     session         1=MIP   0=V5
%     %     trial           trial number
%     %     rt              RT in ms
%     %     tms             0=no 1=stim
%     %     tms_cond        0=no 1=MIP 2=V5
%     %     accuracy        1=corr 0=err nan=empty/distractor
%     %     d_v             distractor value(prob*rew)
%     %     lv_v            low value option value(prob*rew)
%     %     hv_v            high value option value(prob*rew)
%     %     position        positions of stims on screen - 1st column: HV, 2nd column: LV, 3rd column: D - 1=top left, 2=top right, 3=bottom left, 4=bottom right
%     
%     %% caps ones for all partics
% 
%     session=[ones(length(MIP.trial{:}),1); zeros(length(V5.trial{:}),1)];%1=MIP, 0=V5;
%     trial=[MIP.trial{:};V5.trial{:}];  
%     tms=[MIP.tms{:};V5.tms{:}];  
%     rt=[MIP.RT{:};V5.RT{:}];
%     accuracy=[MIP.accuracy{:};V5.accuracy{:}]; 
%     d_v=[MIP.D{:};V5.D{:}]; 
%     lv_v=[MIP.LV{:};V5.LV{:}]; 
%     hv_v=[MIP.HV{:};V5.HV{:}];
%     position=[MIP.pos_rews{:};V5.pos_rews{:}];   
%     
%     %for distractor location
%     locations={[1 3], [2 4]}; %ipsi=[1 3]=left
%     di_loc=zeros(size(hv_v));
%     di_loc(position(:,3)==2 |position(:,3)==4)=1;
%     dloc_for_glm=di_loc;    
%    
%     %collect partics
%     ACCURACY(:,partic)=(accuracy);
%    
%     % exclude   
%     rt_keep=rt;   
%     accuracy_deletion=(isnan(accuracy));%delete all trials where distractor or empty was chosen   
%     tms_deletion=zeros(length(accuracy),1);
%         
%     session(accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max |tms_deletion)=[];        
%     trial(accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
%     rt(accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
%     tms( accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
%     accuracy( accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
%     d_v( accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
%     lv_v( accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
%     hv_v( accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion)=[];        
%     position(accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion,:)=[];   
%     dloc_for_glm(accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max|tms_deletion,:)=[];       
%          
%     clear acc data file MIP V5 temp
%     
%      rt = rt/1000;
%     %     for i = 1:length(rt)
%     %     rt(i)=1/rt(i);
%     %     end
%     %     rt = (rt-min(rt))/(max(rt)-min(rt));
%     %% GLM
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for run=1:length(Conds)
%         if Conds{run}(1:2)=='MI'
%             this_session=1;
%             tms_index=4;
%         else
%             this_session=0;
%             tms_index=3;
%         end
% 
%         if Conds{run}(tms_index)=='1'
%             this_tms=1;
%         else
%             this_tms=0;
%         end
% 
%         if Conds{run}(end-2:end)=='psi'
%             this_loc=0; %ipsi
%         else
%             this_loc=1; %contra
%         end
%         pos_oi=3;
% 
%         warning('') 
%         idx=(session==this_session & tms==this_tms  & (position(:,pos_oi)==locations{this_loc+1}(1) | position(:,pos_oi)==locations{this_loc+1}(2) ));
% 
%         if sum(idx)==0
%             exclude=[exclude,partic];
%         else
%             regressors=[eval(model)];
%             regressors=[normalise(regressors)];
%             criterion=[accuracy(idx)];% Y is a vector of response values.  If DISTR is 'binomial' Y may a binary vector indicating success/failure, and the total number of trials is taken to be 1 for all observations.  If DISTR is 'binomial', Y may also be a two column matrix, the first column containing the number of successes for each observation,and the second containing the total number of trials.
% 
%             [betas_acc(partic,:,run),dev,stats]=glmfit(regressors,criterion,'binomial');
%             
% 
% 
%             idx=(session==this_session & accuracy == 1 & tms==this_tms  & (position(:,pos_oi)==locations{this_loc+1}(1) | position(:,pos_oi)==locations{this_loc+1}(2) ));
%             regressors=[eval(model)];
%             regressors=[normalise(regressors)];criterion=[rt(idx)]; 
%             
%             [betas_rt(partic,:,run),dev,stats]=glmfit(regressors,criterion);
% 
%             [warnMsg, warnId] = lastwarn;
%             if ~isempty(warnMsg)
%                 exclude=[exclude,partic];
%             end
%         end
% 
%     end
% 
% end
% % exclude
% exclude=unique(exclude);
% fprintf('%i participants excluded \n',length(exclude))
% 
% discarded_acc=0;
% all=0;
% for partic=1:length(Partic)
%     discarded_acc=discarded_acc+sum(isnan(ACCURACY(:,partic)));
%     all=all+length(ACCURACY(:,partic));
% end
% fprintf(' %2.2f %% of trials were discarded because they picked nonsense \n',discarded_acc./all.*100)
% 
% betas_acc(exclude,:,:)=[];
% betas_rt(exclude,:,:)=[];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% ANOVA ON RT AND ACC
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% for regressor=4%:length(Regs)+1
%     titles ={};
%     for i = 1:length(Conds)
%         titles{i} = strcat('Acc_',Conds{i});
%         titles{i+length(Conds)} = strcat('RT_',Conds{i});
%     end
%     tab=array2table([squeeze(betas_acc(:,regressor,:)),squeeze(betas_rt(:,regressor,:))],'variablenames',titles);
%     
%     %set up ANOVA
%     within1 = categorical([1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0])';
%     within2 = categorical([0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1])';
%     within3=  categorical([0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1])';
%     within4 = categorical([0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1])';
%     within = table(within1,within2,within3,within4,'variablenames',{'Session','Stim','DLoc','AccRT'});
%     factors='Session*Stim*DLoc*AccRT';
% 
%     allvars=strcat(titles{1},'-',titles{end},' ~1');
%     rm = fitrm(tab,allvars,'WithinDesign',within);
%     ranovatbl = ranova(rm,'withinmodel',factors)
%     
% %     %report
% %     p_spots=[3 5 7 9 11 13 15];
% %     p_spot_names={'Session' 'Stim' 'DLoc' 'Session*Stim' 'Session*DLoc*' 'Stim*DLoc' 'Session*Stim*DLoc'};
% % 
% %     fprintf('\n- - %s - -\n',Regs{regressor-1})
% %     ranova(rm,'withinmodel',factors);
% %     ranovatbl=table2array(ranovatbl);
% %     %%this is to only report significant ones
% %     % if sum(ranovatbl([p_spots],5)<.1)>0
% %     %    sig_idx=find(ranovatbl([p_spots],5)<.1);
% %     %    sig_idx=p_spots(sig_idx);
% %     %    for i=1:length(sig_idx)
% %     %        fprintf('%s : F(%d,%d) = %4.2f, p = %1.3f, n = %4.3f \n',p_spot_names{find(p_spots==sig_idx(i))},  ranovatbl(1,2),ranovatbl(2,2), ranovatbl ((sig_idx(i)),4),ranovatbl((sig_idx(i)),5), (ranovatbl((sig_idx(i)),1)/(ranovatbl((sig_idx(i)),1) + ranovatbl((sig_idx(i))+1,1))))
% %     %    end   
% %     % else
% %     %     disp('nothing <.1')
% %     % end
% %     sig_idx=[p_spots];
% %     for i=1:length(sig_idx)
% %        fprintf('%s : F(%d,%d) = %4.2f, p = %1.3f, n = %4.3f \n',p_spot_names{find(p_spots==sig_idx(i))},  ranovatbl(1,2),ranovatbl(2,2), ranovatbl ((sig_idx(i)),4),ranovatbl((sig_idx(i)),5), (ranovatbl((sig_idx(i)),1)/(ranovatbl((sig_idx(i)),1) + ranovatbl((sig_idx(i))+1,1))))
% %     end   
% 
% end
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% PLOT BETAS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Figure 
% % errors
%     figure
%     subplot(1,2,1)
%     hold on
% betas = betas_acc;
% 
% for c=1:length(Conds)
%     for r=1:length(Regs)+1
%         x=(betas(~isnan(betas(:,r ,c)),r ,c));
%         SEM = std(x)/sqrt(length(x));               % Standard Error
%         errors.(strcat('Reg',num2str(r))).(Conds{c})=[SEM];
%     end
% end
% 
% rfig=4%: length(Regs)+1 %one figure per regressor (except intercept)   
% 
%     to_plot=[];
% 
% 
%     %1 MIP0ips 2 MIP1ips 3 MT0ips 4 MT1ips
%     %5 MIP0con 6 MIP1con 7 MT0con 8 MT1con
%     hold on
%     xlim([0 5])
%      
%     for lines=1:4
%         if lines==1 %MT0ipsi MT0con
%             these_conds=[3 7]; 
%             these_x=[1 2];
%             col=[.5 .5 .5];
%             linetype='--';
%         elseif lines==2 %MT1ipsi MT1con
%             these_conds=[4 8];
%             these_x=[1 2];
%             col=['k'];
%             linetype='-';
%         else
%             these_x=[3 4];
%             if lines==3 %MIp0ipsi MIP0con
%                 these_conds=[1 5];
%                 col=[ 1 .75 .5]; 
%                 linetype='--';
%             else
%                 these_conds=[2 6];
%                 col=[.875 .438 0];
%                 linetype='-';
%             end
%         end
%         to_plot=squeeze(mean(betas(:,rfig,these_conds)));
%         error_to_plot=[([errors.(strcat('Reg',num2str(rfig))).(Conds{these_conds(1)})]); ([errors.(strcat('Reg',num2str(rfig))).(Conds{these_conds(2)})])];
%         
%         plot([0:5], zeros(1,6),':', 'Color','k');%[.7 .7 .7])
%     
%         leg(lines)=plot(these_x,[to_plot],linetype, 'Color', col, 'Linewidth',2);
%         errorbar(these_x(1), to_plot(1), error_to_plot(1), 'Color', col, 'Linewidth',2);
%         errorbar(these_x(2), to_plot(2), error_to_plot(2), 'Color', col, 'Linewidth',2);      
%     end
%         
%     temp={'MT Non-TMS','MT TMS','MIP Non-TMS','MIP TMS'};
%     legend(leg,temp)
%    
%     if rfig <4
%         ylim([-.6 1.2])
%         set(gca,'YTick',[-.5 0 .5 1 1.5])
%     else
%         ylim([-.3 .6])
%         set(gca,'YTick',[-.3 0 .3 .6])
%     end
%     set(gca,'XTick',[1:4])
%     set(gca,'xticklabel',{'Ipsi' 'Contra' 'Ipsi' 'Contra'})
%     ticks=get(gca,'YTick');
%     title(strcat('Accuracy: ',Regs{rfig-1}))
%    
%     
% betas = betas_rt;
% 
% for c=1:length(Conds)
%     for r=1:length(Regs)+1
%         x=(betas(~isnan(betas(:,r ,c)),r ,c));
%         SEM = std(x)/sqrt(length(x));               % Standard Error
%         errors.(strcat('Reg',num2str(r))).(Conds{c})=[SEM];
%     end
% end
% 
% rfig=4%: length(Regs)+1 %one figure per regressor (except intercept)   
% subplot(1,2,2)
% hold on
%     to_plot=[];
% 
%     %1 MIP0ips 2 MIP1ips 3 MT0ips 4 MT1ips
%     %5 MIP0con 6 MIP1con 7 MT0con 8 MT1con
%     hold on
%     xlim([0 5])
%      
%     for lines=1:4
%         if lines==1 %MT0ipsi MT0con
%             these_conds=[3 7]; 
%             these_x=[1 2];
%             col=[.5 .5 .5];
%             linetype='--';
%         elseif lines==2 %MT1ipsi MT1con
%             these_conds=[4 8];
%             these_x=[1 2];
%             col=['k'];
%             linetype='-';
%         else
%             these_x=[3 4];
%             if lines==3 %MIp0ipsi MIP0con
%                 these_conds=[1 5];
%                 col=[ 1 .75 .5]; 
%                 linetype='--';
%             else
%                 these_conds=[2 6];
%                 col=[.875 .438 0];
%                 linetype='-';
%             end
%         end
%         to_plot=squeeze(mean(betas(:,rfig,these_conds)));
%         error_to_plot=[([errors.(strcat('Reg',num2str(rfig))).(Conds{these_conds(1)})]); ([errors.(strcat('Reg',num2str(rfig))).(Conds{these_conds(2)})])];
%         
%         plot([0:5], zeros(1,6),':', 'Color','k');%[.7 .7 .7])
%     
%         leg(lines)=plot(these_x,[to_plot],linetype, 'Color', col, 'Linewidth',2);
%         errorbar(these_x(1), to_plot(1), error_to_plot(1), 'Color', col, 'Linewidth',2);
%         errorbar(these_x(2), to_plot(2), error_to_plot(2), 'Color', col, 'Linewidth',2);      
%     end
%          set(gca,'XTick',[1:4])
%     set(gca,'xticklabel',{'Ipsi' 'Contra' 'Ipsi' 'Contra'})   
%     temp={'MT Non-TMS','MT TMS','MIP Non-TMS','MIP TMS'};
% %     legend(leg,temp)
%    
%     title(strcat('RT: ',Regs{rfig-1}))
%     
%    
% % cd('C:\Users\ckohl\Desktop\')
% % print -depsc HV+LV_line
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% ANOVA JUST ON RT
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% titles=Conds;
% tab=array2table(squeeze(betas_rt(:,regressor,:)),'variablenames',titles);
% 
% %set up ANOVA
% within1 = categorical([1 1 0 0 1 1 0 0 ])';
% within2 = categorical([0 1 0 1 0 1 0 1])';
% within3=  categorical([0 0 0 0 1 1 1 1])';
% within = table(within1,within2,within3,'variablenames',{'Session','Stim','DLoc'});
% factors='Session*Stim*DLoc';
% 
% allvars=strcat(titles{1},'-',titles{end},' ~1');
% rm = fitrm(tab,allvars,'WithinDesign',within);
% ranovatbl = ranova(rm,'withinmodel',factors);
% 
% %report
% p_spots=[3 5 7 9 11 13 15];
% p_spot_names={'Session' 'Stim' 'DLoc' 'Session*Stim' 'Session*DLoc*' 'Stim*DLoc' 'Session*Stim*DLoc'};
% 
% fprintf('\n- - %s - -\n',Regs{regressor-1})
% ranova(rm,'withinmodel',factors);
% ranovatbl=table2array(ranovatbl);
% 
% sig_idx=[p_spots];
% for i=1:length(sig_idx)
% fprintf('%s : F(%d,%d) = %4.2f, p = %1.3f, n = %4.3f \n',p_spot_names{find(p_spots==sig_idx(i))},  ranovatbl(1,2),ranovatbl(2,2), ranovatbl ((sig_idx(i)),4),ranovatbl((sig_idx(i)),5), (ranovatbl((sig_idx(i)),1)/(ranovatbl((sig_idx(i)),1) + ranovatbl((sig_idx(i))+1,1))))
% end   
%     
%     
%  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% ANOVA JUST ON ACC
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
% titles=Conds;
% tab=array2table(squeeze(betas_acc(:,regressor,:)),'variablenames',titles);
% 
% %set up ANOVA
% within1 = categorical([1 1 0 0 1 1 0 0 ])';
% within2 = categorical([0 1 0 1 0 1 0 1])';
% within3=  categorical([0 0 0 0 1 1 1 1])';
% within = table(within1,within2,within3,'variablenames',{'Session','Stim','DLoc'});
% factors='Session*Stim*DLoc';
% 
% allvars=strcat(titles{1},'-',titles{end},' ~1');
% rm = fitrm(tab,allvars,'WithinDesign',within);
% ranovatbl = ranova(rm,'withinmodel',factors);
% 
% %report
% p_spots=[3 5 7 9 11 13 15];
% p_spot_names={'Session' 'Stim' 'DLoc' 'Session*Stim' 'Session*DLoc*' 'Stim*DLoc' 'Session*Stim*DLoc'};
% 
% fprintf('\n- - %s - -\n',Regs{regressor-1})
% ranova(rm,'withinmodel',factors);
% ranovatbl=table2array(ranovatbl);
% 
% sig_idx=[p_spots];
% for i=1:length(sig_idx)
% fprintf('%s : F(%d,%d) = %4.2f, p = %1.3f, n = %4.3f \n',p_spot_names{find(p_spots==sig_idx(i))},  ranovatbl(1,2),ranovatbl(2,2), ranovatbl ((sig_idx(i)),4),ranovatbl((sig_idx(i)),5), (ranovatbl((sig_idx(i)),1)/(ranovatbl((sig_idx(i)),1) + ranovatbl((sig_idx(i))+1,1))))
% end   