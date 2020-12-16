
% Eyetracking_for_Publication_new_ttests.m
%% This dos not replayce Eyeteacking_for_Publication.m, it's an addition to it
%% This is like Eyeteacking_for_Publication.m, but I cut it off and added some t-tests
%% Those are the first stats I report for the eyetracking (based on Bolton's changes)
%% to use all those t-tests, change 'Which_criterion' depending on what you need
%from EyeTrackingIV.m
clear
% runs GLM to see the impact of different predictors on different gaze
% shifts & creates figure 3



Which_GLM=1;    
                % 1: HV+LV, HV-LV, D-HV  => these values to predict the number of gaze shifts in different conds (TMSxSess)
                % 2: HV, LV, D  => these values to predict the number of gaze shifts in different conds (TMSxSess)
                % 3: hv_lv, lv_hv, lv_dv, dv_lv, hv_dv, dv_hv  => these gaze shifts to predict accuracy across all non-tms trials

Which_criterion =  3 ;
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
    load('J:\PolyU\TMS\Data\Navigator\TMS_DEV')
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
                            end
                        end
                    end
                end
            end
        end
    end
end

exclude=unique(exclude);
fprintf('%i participants excluded \n',length(exclude))



betas(exclude,:,:)=[];


%% NEW T TESTS
%% test each regressor against 0
for regressor=2:size(betas,2)
    [h p ci stats]=ttest(betas(:,regressor));
    D=(mean(betas(:,regressor))-0)/std(betas(:,regressor));
    fprintf('\nTesting %s against 0:\n',Regs{regressor-1})
    fprintf(' t(%d) = %2.2f, p = %2.3f, d = %2.2f\n',stats.df,stats.tstat, p,D)
    
end

