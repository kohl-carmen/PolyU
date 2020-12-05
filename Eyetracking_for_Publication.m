%% Based on EyeTracking.m
%from EyeTrackingIV.m
clear
% runs GLM to see the impact of different predictors on different gaze
% shifts & creates figure 3



Which_GLM=1;    
                % 1: HV+LV, HV-LV, D-HV  => these values to predict the number of gaze shifts in different conds (TMSxSess)
                % 2: HV, LV, D  => these values to predict the number of gaze shifts in different conds (TMSxSess)
                % 3: hv_lv, lv_hv, lv_dv, dv_lv, hv_dv, dv_hv  => these gaze shifts to predict accuracy across all non-tms trials

Which_criterion =  1 ;
                % If GLM is 1 or 2 (predicting gaze shifts), which gaze shifts do you wana predict?
                % 1: bidirectional D to HV and HV to D
                % 2: D to HV
                % 3: HV to D


if Which_GLM==1
    Regs={'HV+LV' 'HV-LV' 'D-HV' };
    model='[hv_v(idx)+lv_v(idx), hv_v(idx)-lv_v(idx),d_v(idx)-hv_v(idx)]';
    Testing=2;
    ANOVA=1;
elseif Which_GLM==2   
    Regs={'HV' 'LV' 'DV'};
    model='[hv_v(idx), lv_v(idx), d_v(idx)]';
    Testing=2;
    ANOVA=1;
elseif Which_GLM==3   
    Regs={'hv_lv' 'lv_hv' 'lv_dv' 'dv_lv' 'hv_dv' 'dv_hv'};
    model='[hv_lv(idx),lv_hv(idx),lv_dv(idx), dv_lv(idx), hv_dv(idx), dv_hv(idx)]';
    Testing=1;
    ANOVA=10;
end


if Which_criterion==1
    criterion_name=['bi_dv_hv'];
elseif Which_criterion==2
    criterion_name=['dv_hv'];
elseif Which_criterion==3
	criterion_name=['hv_dv'];
end


Partic=1:31;

Sessions={'MT' 'MIP'};
ppt=0;
file_dir='D:\PolyU\TMS\Data\Tobii\';
StimType={'hv_e' 'lv_e' 'd_e'};


cd('D:\PolyU\TMS\Data\Beh\')
load('magprob_betas_sumdiff')
load('magprob_betas_diff')


contra_only=0;
stepwise=0;
ind=0;
for all_these_loops1 =1%:2%[1:3] %ANOVA
    Testing=Testing;  %1=acc, 2=rt
    tms_dev_cutoff=5;
%     ANOVA=2%all_these_loops1;
    for all_these_loops2 =1% [1:3] %DV
        DV=all_these_loops2;
        for all_these_loops3=0;%[0 1]
            exclude_hard=all_these_loops3;
            exclude_easy=0;
            for all_these_loops4 =0;%[0 1]
                exclude_bad_tms =all_these_loops4;
                for all_these_loops5 =0;% [0 1]
                    delete_trials_without_Dfix=all_these_loops5;
                    for all_these_loops6 =0;% [0 1]
                        late_fixations_only= all_these_loops6;
                        for all_these_loops7=0;%[0 1]
                            delete_extreme_fix= all_these_loops7;
                            for all_these_loops8=0;%[0 1]
                                only_highdv_dfix= all_these_loops8;
                        
                        
  
                        
ind=ind+1;
exclude=[]; 
%exclude bad tms
if exclude_bad_tms
    load('D:\PolyU\TMS\Data\Navigator\TMS_DEV')
    count_tms_exclusions=0;
end
if ANOVA==1
    Conds={'MIP0' 'MIP1' 'MT0' 'MT1'};
elseif ANOVA==2 
    Conds={'MIP0ipsi' 'MIP1ipsi' 'MT0ipsi' 'MT1ipsi' 'MIP0contra' 'MIP1contra' 'MT0contra' 'MT1contra'};
elseif ANOVA==3
    Conds={'MIP0_Dipsi_HVipsi' 'MIP1_Dipsi_HVipsi' 'MT0_Dipsi_HVipsi' 'MT1_Dipsi_HVipsi' ...
           'MIP0_Dcontra_HVipsi' 'MIP1_Dcontra_HVipsi' 'MT0_Dcontra_HVipsi' 'MT1_Dcontra_HVipsi'...
           'MIP0_Dipsi_HVcontra' 'MIP1_Dipsi_HVcontra' 'MT0_Dipsi_HVcontra' 'MT1_Dipsi_HVcontra'...
           'MIP0_Dcontra_HVcontra' 'MIP1_Dcontra_HVcontra' 'MT0_Dcontra_HVcontra' 'MT1_Dcontra_HVcontra'};
elseif ANOVA==10
    Conds={'nontms'};
end
 
betas=[];
count_bad_trials=[];
for partic=1:length(Partic)
%   partic=partic+1
    sess_i=[];
    nr_fix=[];
    dur_fix=[];
    rel_dur_fix=[];
    tms=[];
    dloc=[];
    hvloc=[];
    hv_v=[];
    lv_v=[];
    d_v=[];
    acc=[];
    rt=[];
    pos=[];
    weightediff=[];
    bad_trials=[];
    for sess=1:length(Sessions)
        subj_code=7000+Partic(partic)*10+sess-1;
        % load fix structure
        loaddir=strcat(file_dir,num2str(subj_code),'\Preproc\');
        cd(loaddir)
        load('FIXATIONS')        
        % first, let's see which stimuli are fixated the most

                
        temp=length(nr_fix);
        for i=1:length(FIXATIONS.nr_fix)
            nr_fix(temp+i,:)=FIXATIONS.nr_fix{i};
            dur_fix(temp+i,:)=FIXATIONS.total_dur{i};
            rel_dur_fix(temp+i,:)=FIXATIONS.total_dur{i}./sum(FIXATIONS.total_dur{i});
            
            
%             if FIXATIONS.delete{i} ==1
%              bad_trials=[bad_trials,temp+i];
%             end
            
            %gaze
            gaze{temp+i}=[];
            if sum(FIXATIONS.nr_fix{i})>0
                cont=FIXATIONS.timeseries{i};
                latest_fix=0;
                for j=2:length(cont)
                    for k=1:4
                        if cont(j,k)==1 & (cont(j-1,k)==0 | j==2)
                           gaze{temp+i}=[gaze{temp+i};latest_fix, k];
                           latest_fix=k;
                        end
                    end
                end
            end
      
            
        
        end
        
   
    
        
       
        if late_fixations_only
            for trial=1:length(FIXATIONS.timeseries)
                mid=round(length(FIXATIONS.timeseries{trial})/2);
                dur_fix(trial,:)=sum(FIXATIONS.timeseries{trial}(mid:end,:));
                rel_dur_fix(trial,:)=dur_fix(trial,:)./sum(dur_fix(trial,:));
            end
        end
        %match trials to beh data
        load(strcat(file_dir(1:18),'Beh\',num2str(subj_code),'\transformed\data_behavior_tms.mat'))       
        tms=[tms;data.behavior.tms{:}];
        dloc_temp=data.behavior.pos_rews{:}(:,3);
        hvloc_temp=data.behavior.pos_rews{:}(:,1);
        pos=[pos;data.behavior.pos_rews{:}];
        weightedhv=data.behavior.vals{:}(:,1).*magprob_betas_sumdiff(partic,2) + data.behavior.probs{:}(:,1).*magprob_betas_sumdiff(partic,3);
        weightedlv=data.behavior.vals{:}(:,2).*magprob_betas_sumdiff(partic,2) + data.behavior.probs{:}(:,2).*magprob_betas_sumdiff(partic,3);
        weightediff=[weightediff;weightedhv-weightedlv];
        
        
        if isempty(dloc)
            dloc(dloc_temp==1 |dloc_temp==3)=0;%ipsi
            dloc(dloc_temp==2 |dloc_temp==4)=1;%contra
        else
            dloc2(dloc_temp==1 |dloc_temp==3)=0;%ipsi
            dloc2(dloc_temp==2 |dloc_temp==4)=1;%contra
            dloc=[dloc;dloc2'];
        end
        
        if isempty(hvloc)
            hvloc(hvloc_temp==1 |hvloc_temp==3)=0;%ipsi
            hvloc(hvloc_temp==2 |hvloc_temp==4)=1;%contra
        else
            hvloc2(hvloc_temp==1 |hvloc_temp==3)=0;%ipsi
            hvloc2(hvloc_temp==2 |hvloc_temp==4)=1;%contra
            hvloc=[hvloc;hvloc2'];
        end
        
        if size(dloc,1)==1
            dloc=dloc';
        end
        
        if size(hvloc,1)==1
            hvloc=hvloc';
        end
        hv_v=[hv_v;data.behavior.HV{:}];
        lv_v=[lv_v;data.behavior.LV{:}];
        d_v=[d_v;data.behavior.D{:}];
        acc=[acc; data.behavior.accuracy{:}];
        rt=[rt; data.behavior.RT{:}];
        newdiff=(magprob_betas_sumdiff(partic,2).*hv_v+lv_v) + (magprob_betas_sumdiff(partic,3).*hv_v-lv_v);
        
%         if isempty(sess_i)
%             sess_i=ones(size(dloc_temp));
%         else
%             sess_i=[sess_i;ones(size(dloc_temp)).*2];
%         end
        if sess==1
            sess_i=ones(size(hv_v));
        else
            sess_i=[sess_i;ones(size(sess_i)).*2];
        end
    end
    
     hv_lv=zeros(540,1);
    hv_dv=zeros(540,1);
    hv_all=zeros(540,1);

    lv_hv=zeros(540,1);
    lv_dv=zeros(540,1);
    lv_all=zeros(540,1);

    dv_hv=zeros(540,1);
    dv_lv=zeros(540,1);
    dv_all=zeros(540,1);
    
    all_dv=zeros(540,1);
    all_hv=zeros(540,1);
    all_lv=zeros(540,1);
    
    %bidirectional
    every_hv=zeros(540,1);
    every_lv=zeros(540,1);
    every_dv=zeros(540,1);
    
    bi_dv_hv=zeros(540,1);
    bi_dv_lv=zeros(540,1);
    bi_hv_lv=zeros(540,1);

    for trial=1:540
        if ~isempty(gaze{trial})
            for i=2:size(gaze{trial},1)
                if gaze{trial}(i,1)==1
                    hv_all(trial)=hv_all(trial)+1;
                    every_hv(trial)=every_hv(trial)+1;

                    if  gaze{trial}(i,2)==2
                        hv_lv(trial)=hv_lv(trial)+1;
                        bi_hv_lv(trial)=bi_hv_lv(trial)+1;
                    elseif gaze{trial}(i,2)==3
                        hv_dv(trial)=hv_dv(trial)+1;
                        bi_dv_hv(trial)=bi_dv_hv(trial)+1;
                    end

                elseif gaze{trial}(i,1)==2
                    lv_all(trial)=lv_all(trial)+1;
                    every_lv(trial)=every_lv(trial)+1;

                     if  gaze{trial}(i,2)==1
                        lv_hv(trial)=lv_hv(trial)+1;
                        bi_hv_lv(trial)=bi_hv_lv(trial)+1;
                    elseif gaze{trial}(i,2)==3
                        lv_dv(trial)=lv_dv(trial)+1;
                        bi_dv_lv(trial)=bi_dv_lv(trial)+1;
                    end


                elseif gaze{trial}(i,1)==3
                    dv_all(trial)=dv_all(trial)+1;
                    every_dv(trial)=every_dv(trial)+1;

                     if  gaze{trial}(i,2)==1
                        dv_hv(trial)=dv_hv(trial)+1;
                         bi_dv_hv(trial)=bi_dv_hv(trial)+1;
                    elseif gaze{trial}(i,2)==2
                        dv_lv(trial)=dv_lv(trial)+1;
                        bi_dv_lv(trial)=bi_dv_lv(trial)+1;
                    end
                end
                
                if gaze{trial}(i,2)==1
                    all_hv(trial)=all_hv(trial)+1;
                    every_hv(trial)=every_hv(trial)+1;
                elseif gaze{trial}(i,2)==3
                    all_dv(trial)=all_dv(trial)+1;
                    every_dv(trial)=every_dv(trial)+1;
                elseif gaze{trial}(i,2)==2
                    all_lv(trial)=all_lv(trial)+1;
                    every_lv(trial)=every_lv(trial)+1;
                end
                    
            end
        end
    end
    
    LV_V(:,partic)=lv_v; 
    HV_V(:,partic)=hv_v; 
    if DV==1
        D_E(:,partic)=nr_fix(:,3);
    elseif DV==2
        D_E(:,partic)=dur_fix(:,3);
    elseif DV==3
        D_E(:,partic)=rel_dur_fix(:,3);
    end
    ACCURACY_BOLTON(:,partic) =acc;
    SESSION(:,partic)=sess_i-1;
    TMS(:,partic)=tms;
    POSITION(:,:,partic)=pos;
    WEIGHTEDIFF(:,partic)=weightediff;
    
    dtemp=zeros(size(tms));
    
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
        tms_deletion=tms_dev>tms_dev_cutoff;% & sess_i==2;
        count_tms_exclusions=count_tms_exclusions+sum(tms_deletion);
    else
        tms_deletion=zeros(length(tms),1);
    end
        
        %find trials with no fixations detected in roi
%         delete=find(sum(nr_fix,2)==0);
% delete=[];
%         delete=find(sum(nr_fix(:,1:3),2)==0 | sum(nr_fix(:,1:3),2)> mean(sum(nr_fix(:,1:3),2))+s(sum(nr_fix(:,1:3),2)).*2);
        delete=find( tms_deletion==1);
%         delete=find(sum(dur_fix(:,1:3),2)<5 | tms_deletion==1);
        if delete_trials_without_Dfix
            delete=find(dur_fix(:,3)<5 | sum(dur_fix(:,1:3),2)<5 | tms_deletion==1);
        end
        
        if delete_extreme_fix
%             delete=[delete; find( sum(dur_fix(:,1:3),2)> mean(sum(dur_fix(:,1:3),2))+std(sum(dur_fix(:,1:3),2)).*2)];
             delete=[delete; find(sum(nr_fix(:,1:3),2)>= 10 | sum(dur_fix(:,1:3),2)> mean(sum(dur_fix(:,1:3),2))+std(sum(dur_fix(:,1:3),2)).*2)];
       
        end
%         delete=find(sum(nr_fix,2)==0 | tms_deletion==1);
%         delete=find(sum(dur_fix(:,1:3),2)<5 | dur_fix(:,3)<5 |tms_deletion==1);

        if only_highdv_dfix
            
%             delete=[delete; find(d_v < quantile(d_v,.6))];
            delete=[delete; find(d_v < lv_v)];
% delete=[delete; find(d_v < median(d_v))];
        end
        
        delete=[delete; bad_trials'];
        count_bad_trials=[count_bad_trials,length(bad_trials)];
        
        nr_fix(delete,:)=[];
        dur_fix(delete,:)=[];
        rel_dur_fix(delete,:)=[];
        tms(delete,:)=[];
        dloc(delete,:)=[];
        hvloc(delete,:)=[];
        hv_v(delete,:)=[];
        lv_v(delete,:)=[];
        d_v(delete,:)=[];
        sess_i(delete,:)=[];
        acc(delete,:)=[];
        rt(delete,:)=[];
        weightediff(delete,:)=[];
        newdiff(delete,:)=[];
        
        dtemp(delete)=1;
        
        
             
    hv_lv(delete,:)=[];
    hv_dv(delete,:)=[];
    hv_all(delete,:)=[];

    lv_hv(delete,:)=[];
    lv_dv(delete,:)=[];
    lv_all(delete,:)=[];

    dv_hv(delete,:)=[];
    dv_lv(delete,:)=[];
    dv_all(delete,:)=[];
    
    all_dv(delete,:)=[];
    all_hv(delete,:)=[];
    all_lv(delete,:)=[];
    
    every_hv(delete,:)=[];
    every_lv(delete,:)=[];
    every_dv(delete,:)=[];
    bi_dv_hv(delete,:)=[];
    bi_dv_lv(delete,:)=[];
    bi_hv_lv(delete,:)=[];
    

       
    DELETE(:,partic)=dtemp;
    
    
    if DV==1
        hv_e=nr_fix(:,1);
        lv_e=nr_fix(:,2);
        d_e=nr_fix(:,3);
        x_e=nr_fix(:,4);
        all_e=sum(nr_fix,2);
    elseif DV==2
        hv_e=dur_fix(:,1);
        lv_e=dur_fix(:,2);
        d_e=dur_fix(:,3);
        x_e=dur_fix(:,4);
        all_e=sum(dur_fix,2);
    elseif DV==3
        hv_e=rel_dur_fix(:,1);
        lv_e=rel_dur_fix(:,2);
        d_e=rel_dur_fix(:,3);
        x_e=rel_dur_fix(:,4);
        all_e=sum(dur_fix,2);
    end
    
    
    
       
    
    for run=1:length(Conds)
        if length(Conds)==1
            idx=tms==0;
        else
        if any(run==[1 2 5 6 9 10 13 14]) % MIP
            this_session=2;
        elseif any(run==[3 4 7 8 11 12 15 16]) % MT
            this_session=1;
        end
        if any(run==[1 3 5 7 9 11 13 15]) % no stim
            this_tms=0;%no stim
        elseif any(run==[2 4 6 8 10 12 14 16]) %  stim
            this_tms=1;
        end
        if any(run==[1 2 3 4]) % d (or whatever, seeANOVA2_val) ipsilateral
            d_loc_oi=0;
            hv_loc_oi=0;
        elseif any(run==[5 6 7 8]) %  d (or whatever, seeANOVA2_val)contralateral
            d_loc_oi=1;
            hv_lov_oi=0;
        elseif any(run==[9 10 11 12])
            d_loc_oi=0;
            hv_lov_oi=1;
        elseif any(run==[13 14 15 16])
            d_loc_oi=1;
            hv_lov_oi=1;
        end

        if ANOVA==1
            idx= sess_i==this_session & tms==this_tms & ~isnan(acc);
        elseif ANOVA==2 
            idx= sess_i==this_session & tms==this_tms & dloc==d_loc_oi & ~isnan(acc) ;
        elseif ANOVA==3
            idx= sess_i==this_session & tms==this_tms & dloc==d_loc_oi & hvloc==hv_loc_oi & ~isnan(acc) ;
        end

        if exclude_hard
            idx= idx & hv_v-lv_v>median(hv_v-lv_v);
        end
        if contra_only
            idx=idx & dloc==1;
        end
        end
        
        if sum(idx)==0
            exclude=[exclude,partic];
        else
            regressors=[eval(model)];

            regressors=[normalise(regressors)];


            if Testing==1
                warning('') 
                criterion=[acc(idx)];
                
                [betas(partic,:,run),dev,stats]=glmfit(regressors,criterion,'binomial');
                [warnMsg, warnId] = lastwarn;
                if ~isempty(warnMsg)
                    exclude=[exclude,partic];
                end
                 if stepwise
                     warning('') 
                    regressors=[eval(model2)];

                    regressors=[normalise(regressors)];

                    criterion=stats.resid;
                    [betas2(partic,:,run),dev,stats]=glmfit(regressors,criterion);
                    [warnMsg, warnId] = lastwarn;
                    if ~isempty(warnMsg)
                        exclude=[exclude,partic];
                    end
                    if stepwise3
                        warning('') 
                        regressors=[eval(model3)];
                        regressors=[normalise(regressors)];
                        criterion=stats.resid;
                        [betas3(partic,:,run),dev,stats]=glmfit(regressors,criterion);
                        [warnMsg, warnId] = lastwarn;
                        if ~isempty(warnMsg)
                            exclude=[exclude,partic];
                        end
                    end
                 end
            else
                criterion=eval(criterion_name);
                criterion=[criterion(idx)];
%                 criterion=[rt(idx)];% Y is a vector of response values.  If DISTR is 'binomial' Y may a binary vector indicating success/failure, and the total number of trials is taken to be 1 for all observations.  If DISTR is 'binomial', Y may also be a two column matrix, the first column containing the number of successes for each observation,and the second containing the total number of trials.
                [betas(partic,:,run),dev,stats]=glmfit(regressors,criterion);
                if stepwise
                    regressors=[eval(model2)];

                    regressors=[normalise(regressors)];

                    criterion=stats.resid;
                    [betas2(partic,:,run),dev,stats]=glmfit(regressors,criterion);
                end
            end
%             if any(betas(partic,:,run) >10)
%              exclude=[exclude,partic];
%             end
        end  
    end
end
exclude=unique(exclude);
fprintf('%i participants excluded \n',length(exclude))

if length(exclude)<15
%% plot betas
clf
if stepwise==1 
    betas(:,end+1,:)=betas2(:,end,:);
    Regs={Regs{:},Regs2{:}};
    if stepwise3
        betas(:,end+1,:)=betas3(:,end,:);
        Regs={Regs{:},Regs3{:}};
    end
end

betas(exclude,:,:)=[];


%% test each regressor against 0
%only non-tms
for regressor=2:size(betas,2)
    %MIP
    [h p ci stats]=ttest(betas(:,regressor));
    fprintf('\nTesting %s against 0:\n',Regs{regressor-1})
    fprintf(' t(%d) = %2.2f, p = %2.3f\n',stats.df,stats.tstat, p)

end

for c=1:length(Conds)
    for r=1:length(Regs)+1
        x=(betas(~isnan(betas(:,r ,c)),r ,c));
        SEM = std(x)/sqrt(length(x));               % Standard Error
        ts = tinv([0.025  0.975],length(x)-1);      % T-Score
        CI = mean(x) + ts*SEM;                      % Confidence Intervals
%         errors.(strcat('Reg',num2str(r))).(Conds{c})=[CI(2)-mean(x)];
        errors.(strcat('Reg',num2str(r))).(Conds{c})=SEM;
    end
end

to_plot=[];
for r=1:length(Regs)+1
    for c=1:length(Conds)
        to_plot(r,c)=([nanmean(betas(:,r ,c))]);
    end
end

if ANOVA==10
    to_plot=to_plot(2:end);
end
h = bar(to_plot,.9);
% colours={ [0 0 1] [.5 1 0]  [ 0 .188 .375] [0 .375 .375]};
% for i=1:2%length(titles)
%      h(i).FaceColor=colours{i};
% end
hold on
%errors
for r=1:length(Regs)+1
    for c=1:length(Conds)
        error_to_plot(r,c)=([errors.(strcat('Reg',num2str(r))).(Conds{c})]);
    end
end
if ANOVA==10
    error_to_plot=error_to_plot(2:end);
end
[numgroups,numbars ] = size(to_plot); 
groupwidth = min(0.8, numbars/(numbars+1.5));
for i = 1:numbars
  % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange: https://uk.mathworks.com/matlabcentral/answers/102220-how-do-i-place-errorbars-on-my-grouped-bar-graph-using-function-errorbar-in-matlab-7-13-r2011b
  x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
  errorbar(x, to_plot(:,i), error_to_plot(:,i), 'k', 'linestyle', 'none');
end
        
hold on
set(gca,'XTick',1:length(Regs)+1)
set(gca,'xticklabel',{'INTC',Regs{:}})
xlabel('Regressors')
ylabel('Beta Value')
legend(Conds, 'Location','southwest')
if ANOVA==10
    ylim([-0.25 0.3])
    set(gca,'YTick',[-.25  0 .25])
        
    cd('D:\PolyU\TMS\Paper')
    print -depsc allgazes
else
    if DV==1
        title(strcat('Nr Fixations (mean/trial)'))
    elseif DV==2
        title(strcat('Fixation Duration (mean/trial)'))
    elseif DV==3
        title(strcat('Relative Fixation Duration (mean/trial)'))
    end
end

if ANOVA<10

%% anovas
if length(Conds)>1
result_text=[];
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
        elseif ANOVA==3
            within1 = categorical([1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0])';
            within2 = categorical([0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1])';
            within3=  categorical([0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1])'; 
            within4 = categorical([0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1])';
            within = table(within1,within2,within3,within4,'variablenames',{'Session','Stim','DLoc', 'HVLoc'});
            factors='Session*Stim*DLoc*HVLoc';
        end
        allvars=strcat(titles{1},'-',titles{end},' ~1');
        rm = fitrm(tab,allvars,'WithinDesign',within);
        ranovatbl = ranova(rm,'withinmodel',factors);
        if ANOVA==2
            p_spots=[3 5 7 9 11 13 15];
            p_spot_names={'Session' 'Stim' 'DLoc' 'Session*Stim' 'Session*DLoc*' 'Stim*DLoc' 'Session*Stim*DLoc'};
        elseif ANOVA==1
            p_spots=[3 5 7];
            p_spot_names={'Session' 'Stim' 'Session*Stim*'};
        elseif ANOVA==3 
            p_spots=[3 5 7 9 11 13 15 17 19 21 23 25 27 29 31];
            p_spot_names={'Session' 'Stim' 'DLoc' 'HVLoc' 'Session*Stim' 'Session*DLoc*' 'Stim*DLoc' 'Session*HVLoc*' 'Stim*HVLoc' 'DLoc*HVLoc' 'Session*Stim*DLoc' 'Session*Stim*HVLoc' 'Session*HVLoc*DLoc' 'Stim*HVLoc*DLoc' 'Session*Stim*HVLoc*DLoc'};
        end
        fprintf('\n- - %s - -\n',Regs{regressor-1})
        result_text=sprintf('%s \n- - %s - -\n'  , result_text,Regs{regressor-1});
        ranova(rm,'withinmodel',factors);
        ranovatbl=table2array(ranovatbl);
        if sum(ranovatbl([p_spots],5)<.1)>0
           sig_idx=find(ranovatbl([p_spots],5)<.1);
           sig_idx=p_spots(sig_idx);
           for i=1:length(sig_idx)
               fprintf('%s : F(%d,%d) = %4.2f, p = %1.3f \n',p_spot_names{find(p_spots==sig_idx(i))},  ranovatbl(1,2),ranovatbl(2,2), ranovatbl ((sig_idx(i)),4),ranovatbl((sig_idx(i)),5))
               result_text=sprintf('%s %s : F(%d,%d) = %4.2f, p = %1.3f \n'  , result_text,p_spot_names{find(p_spots==sig_idx(i))},  ranovatbl(1,2),ranovatbl(2,2), ranovatbl ((sig_idx(i)),4),ranovatbl((sig_idx(i)),5));
           end   
        else
            disp('nothing <.1')
            result_text=sprintf('%s \n nothing <.1 \n',result_text );
        
        end
        

end
                
fprintf('\nANOVA:%i',ANOVA)
fprintf('\nDV:%i',DV)
fprintf('\nTesting:%i',Testing)
fprintf('\nStepwise:%i',stepwise)
fprintf('\nExclude hard:%i',exclude_hard)
if exclude_bad_tms
    fprintf('Exclude based on bad tms: cutoff %i -> lost %i trials \n', tms_dev_cutoff,count_tms_exclusions)
end
fprintf('\nExclude Trials based on Dfix:%i',delete_trials_without_Dfix)
fprintf('\nExclude early fixs:%i',late_fixations_only),
fprintf('\n')

setting_text=sprintf('\nANOVA:%i \nDV:%i  \nTesting:%i \nStepwise:%i \nExclude Easy:%i \nExclude Trials based on Dfix:%i \nExclude early fixs:%i \nExclude based on bad tms:%i \n Exclude extreme fixations:%i \n Only highDV:%i \n',ANOVA ,DV ,Testing,stepwise,exclude_hard ,delete_trials_without_Dfix ,late_fixations_only,exclude_bad_tms, delete_extreme_fix, only_highdv_dfix );

display_text=sprintf('%s \n\ %s',result_text,setting_text);

%If I wana add just one slide

end
end
end

if ANOVA<10
close all
for rfig=2: length(Regs)+1 %one figure per regressor (except intercept)
    figure(rfig-1)
    to_plot=[]; 
    error_to_plot=[];
    for MTMIP=1:2
        if MTMIP==1
            these_conds=[3 4];
        else
            these_conds=[1 2];
        end
        for cond=1:length(Conds)/2% ipsi0 contra0 ipsi1 contra1
            to_plot(MTMIP,cond)=mean(betas(:,rfig,these_conds(cond)));
            error_to_plot(MTMIP,cond)=([errors.(strcat('Reg',num2str(rfig))).(Conds{these_conds(cond)})]);
    
        end
    end
    h = bar(to_plot,.9);
    colours={[ 1 .75 .5]  [.875 .438 0]};
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

    if rfig==2 %| rfig==4
        temp1=sprintf('MT Non-TMS');
        temp2=sprintf('MT TMS');
        temp3=sprintf('MIP Non-TMS');
        temp4=sprintf('MIP TMS');
%         legend({'Non-TMS  D-ipsilateral'  'Non-TMS  D-contralateral' 'TMS  D-ipsilateral' 'TMS D-contralateral'})
         legend(temp1, temp2, temp3, temp4)
    end
    

if criterion_name(1)=='b'
%         ylim([-.11 .07])
        ylim([-.04 .1])
        set(gca,'YTick',[-.1 -.05 0 .05 .1])
elseif criterion_name(1)=='d'
%         ylim([-.07 .05])
        ylim([-.04 .07])
        set(gca,'YTick',[-.04 -.02 0 .02 .04 .06])
elseif criterion_name(1)=='h'
%         ylim([-.05 .025])
        ylim([-.02 .04])
        set(gca,'YTick',[-.02  0 .02 .04])
end
        

    set(gca,'xticklabel',{'MT' 'MIP'})
    ticks=get(gca,'YTick');
    title(Regs{rfig-1})
    
    
    cd('D:\PolyU\TMS\Paper')
if criterion_name(1)=='b'
       if rfig==2
           print -depsc eyeHV+LV_bar
       elseif rfig==3
           print -depsc eyeLV-LV_bar
       elseif rfig==4
           print -depsc eyeDV-HV_bar
       end
elseif criterion_name(1)=='d'
    if rfig==2
           print -depsc eyeHV+LV_bar_dvhv
       elseif rfig==3
           print -depsc eyeLV-LV_bar_dvhv
       elseif rfig==4
           print -depsc eyeDV-HV_bar_dvhv
       end
elseif criterion_name(1)=='h'
     if rfig==2
           print -depsc eyeHV+LV_bar_hvdv
       elseif rfig==3
           print -depsc eyeLV-LV_bar_hvdv
       elseif rfig==4
           print -depsc eyeDV-HV_bar_hvdv
     end
end

end
end

%% line plots
% close all
% for rfig=2: length(Regs)+1 %one figure per regressor (except intercept)
%     
%     figure(rfig-1)
%     to_plot=[];
% 
%      %plot interaction
%     clf
%     %1 MIP0ips 2 MIP1ips 3 MT0ips 4 MT1ips
%     %5 MIP0con 6 MIP1con 7 MT1con 8 MT0con
%     hold on
%     xlim([0 3])
% %     ylim([-.05 .05])
%      
%     for lines=1:2
%         if lines==1 %MT0 MIP1
%             these_conds=[3 4]; 
%             these_x=[1 2];
%             col=['k'];
%             linetype='-';
%         elseif lines==2 %MIP0 MIP1
%             these_conds=[1 2];
%             these_x=[1 2];
%             col=[.875 .438 0];
%             linetype='-';
%         end
%         to_plot=squeeze(mean(betas(:,rfig,these_conds)));
%         error_to_plot=[([errors.(strcat('Reg',num2str(rfig-1))).(Conds{these_conds(1)})]); ([errors.(strcat('Reg',num2str(rfig-1))).(Conds{these_conds(2)})])];
%         
%         plot([0:5], zeros(1,6),':', 'Color','k');%[.7 .7 .7])
%     
%         leg(lines)=plot(these_x,[to_plot],linetype, 'Color', col, 'Linewidth',2)
%         errorbar(these_x(1), to_plot(1), error_to_plot(1), 'Color', col, 'Linewidth',2);
%         errorbar(these_x(2), to_plot(2), error_to_plot(2), 'Color', col, 'Linewidth',2);
% 
%         
%     end
%         
%       if rfig==2 | rfig==4
%         temp1=sprintf('MT');
%         temp2=sprintf('MIP');
% %         legend({'Non-TMS  D-ipsilateral'  'Non-TMS  D-contralateral' 'TMS  D-ipsilateral' 'TMS D-contralateral'})
%         legend(leg,temp1, temp2)
%       end
%     
% %     if rfig <4
% %       ylim([-.6 1.2])
% %         set(gca,'YTick',[-.5 0 .5 1 1.5])
% %     else
%         ylim([-.11 .07])
% %         set(gca,'YTick',[-.3 0 .3 .6])
% %     end
%     set(gca,'XTick',[1:2])
%     set(gca,'xticklabel',{'Non-TMS' 'TMS'})
%     ticks=get(gca,'YTick');
%     title(Regs{rfig-1})
%     
%     
% %        if rfig==2
% %            print -depsc HV+LV_line
% %        elseif rfig==3
% %            print -depsc HV-LV_line
% %        elseif rfig==4
% %            print -depsc D-HV_line
% %        end
% 
% end



                            end
                        end
                    end
                end
            end
        end
    end
end


if ANOVA==10
%against 0
for i=2:7
    [binary,p,ci,tstat]=ttest(betas(:,i));
    Mean1= mean(betas(:,i));
    Mean2= 0;
    SD1 = std(betas(:,i));
    SD2 = 0;
    D= (Mean2-Mean1)/(sqrt(((SD1)^2 +(SD2)^2)/2));
    Regs{i-1}
    fprintf('t(%i) = %2.2f, p = %1.3f, D = %2.3f', tstat.df, tstat.tstat, p,D)

    
end
end


% 
% 
% 
% 
% %post hoc
% tab
% %fix session to MT (then MIP) to see if stim matters
% 
% % take only MIP
% tabbi=table2array(tab);
% tabbi=tabbi(:,1:2);
% % ttest between tms and non-tms
% [h p ci]=ttest(tabbi(:,1),tabbi(:,2))
% 
% % take only MT
% tabbi=table2array(tab);
% tabbi=tabbi(:,3:4);
% % ttest between tms and non-tms
% [h p ci]=ttest(tabbi(:,1),tabbi(:,2))
