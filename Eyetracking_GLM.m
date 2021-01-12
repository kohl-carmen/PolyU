%% Main Gaze Analysis
% most recently, based on Eyetracking_for_Publication.m
% (EyeTracking.m/EyeTrackingIV.m)

% runs GLM to see the impact of different predictors on different gaze
% shifts & creates figure 3

% This does a number of things
% Set 'Which_GLM' and 'Which_criterion' variables depending on what
% analysis to run
%   - Which_GLM = 3 &  Which_criterion =  1
%           - Runs GLM3: {'HV+LV' 'HV-LV' 'D-HV' } to predict bidirectional 
%             (D->HV and HV->D) gaze shifts per Stimualtion (tms/non) and 
%             Sesssion (MIP/MT)
%           - Runs and reportst test of beta weights (non-tms conds) against zero
%           - Runs and reports Stimualtion (tms/non) x Sesssion (MIP/MT) 
%             ANOVA on beta weights
%           - Reports follow-up t-tests (test effect of TMS for MIP and MT
%             separately)
%           - Plots Figure 5 a (note bar order was changed in illustrator)
%  
%   - Which_GLM = 3 &  Which_criterion =  2
%           -  Same as above, but now the GLM is predicting gaze shifts
%              from D to HV (not bidirectional)
%           -  Still reports everything liek above, but we're only
%              reporting one ANOVA interaction in the manuscript
%           -  Plots Figure 5 b (note bar order was changed in illustrator)
% 
%   - Which_GLM = 3 &  Which_criterion =  2
%           -  Same as above, but now the GLM is predicting gaze shifts
%              from HV to D (not bidirectional)
%           -  Still reports everything like above, but we're only
%              reporting one ANOVA interaction in the manuscript
%           -  Plots Figure 5 c (note bar order was changed in illustrator)
% 
%   - Which_GLM = 4 (Which_criterion doesn't matter)
%           -  Runs GLM4: {'hv_lv' 'lv_hv' 'lv_dv' 'dv_lv' 'hv_dv' 'dv_hv'}
%              gaze shifts to predict accuracy across all non-tms trials
%           -  Runs and reports ttests for each beta against zero
%           -  Plots Figure 5 d
% 
% Two more things this does:
% 
%   - Which_GLM = 5 (for each criterion)
%           -  Runs a GLM: {'HV' 'LV' 'D'} to predict gaze shifts
%           -  Everything else like above
%           -  Plots Supplementary Figure S3 (a/b/c  with Which_criterion
%               = 1/2/3 respectively)
%
%   - Manually uncomment lines 115-116 (while Which_GLM=3, Which_criterion=1/2/3)
%           - Tests that results do not change qualitatively when we
%             use the original Site x Stim x Session ANOVA


clear

file_dir='D:\PolyU\TMS\Data\Tobii\';
output_dir='C:\Users\ckohl\Desktop';


Which_GLM = 4;    
                % 3: HV+LV, HV-LV, D-HV  => these values to predict the number of gaze shifts in different conds (TMSxSess)
                % 4: hv_lv, lv_hv, lv_dv, dv_lv, hv_dv, dv_hv  => these gaze shifts to predict accuracy across all non-tms trials
                % 5: % 2: HV, LV, D  => these values to predict the number of gaze shifts in different conds (TMSxSess)

Which_criterion =  1;
                % If GLM is 1 or 2 (predicting gaze shifts), which gaze shifts do you wana predict?
                % 1: bidirectional D to HV and HV to D
                % 2: D to HV
                % 3: HV to D


                
if Which_GLM==3
    Regs={'HV+LV' 'HV-LV' 'D-HV' };
    model='[hv_v(idx)+lv_v(idx), hv_v(idx)-lv_v(idx),d_v(idx)-hv_v(idx)]';
    Testing=2;
    ANOVA=1;
elseif Which_GLM==4   
    Regs={'hv_lv' 'lv_hv' 'lv_dv' 'dv_lv' 'hv_dv' 'dv_hv'};
    model='[hv_lv(idx),lv_hv(idx),lv_dv(idx), dv_lv(idx), hv_dv(idx), dv_hv(idx)]';
    Testing=1;
    ANOVA=10;
elseif Which_GLM==5
    Regs={'HV' 'LV' 'DV'};
    model='[hv_v(idx), lv_v(idx), d_v(idx)]';
    Testing=2;
    ANOVA=1;
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
StimType={'hv_e' 'lv_e' 'd_e'};

                        
exclude=[]; 

if ANOVA==1
    Conds={'MIP0' 'MIP1' 'MT0' 'MT1'};
elseif ANOVA==10
    Conds={'nontms'};
end

% only for control 
% ANOVA=2;
% Conds={'MIP0ipsi' 'MIP1ipsi' 'MT0ipsi' 'MT1ipsi' 'MIP0contra' 'MIP1contra' 'MT0contra' 'MT1contra'};

fprintf('\n ---------------------')
fprintf('\n ---------------------\n')
fprintf('  -- GLM%i --\n',Which_GLM)
if Which_GLM==3
    fprintf('  -- Direction: %s --\n',criterion_name)
end
fprintf('  -- %i Conds --',length(Conds))
fprintf('\n ---------------------')
fprintf('\n ---------------------\n')
betas=[];
for partic=1:length(Partic)
    sess_i=[];    nr_fix=[];
    tms=[];       dloc=[];
    hvloc=[];     hv_v=[];
    lv_v=[];      d_v=[];
    acc=[];       pos=[];

    for sess=1:length(Sessions)
        subj_code=7000+Partic(partic)*10+sess-1;
        
        % load fixations structure
        loaddir=strcat(file_dir,num2str(subj_code),'\Preproc\');
        cd(loaddir)
        load('FIXATIONS')        
        % get gaze data        
        temp=length(nr_fix);
        for i=1:length(FIXATIONS.nr_fix)
            nr_fix(temp+i,:)=FIXATIONS.nr_fix{i};          
            %gaze
            gaze{temp+i}=[];
            if sum(FIXATIONS.nr_fix{i})>0
                cont=FIXATIONS.timeseries{i};
                latest_fix=0;
                for j=2:length(cont)
                    for k=1:4
                        if cont(j,k)==1 && (cont(j-1,k)==0 || j==2)
                           gaze{temp+i}=[gaze{temp+i};latest_fix, k];
                           latest_fix=k;
                        end
                    end
                end
            end          
        end
  
        %match trials to beh data
        load(strcat(file_dir(1:18),'Beh\',num2str(subj_code),'\transformed\data_behavior_tms.mat'))       
        tms=[tms;data.behavior.tms{:}];
        dloc_temp=data.behavior.pos_rews{:}(:,3);
        hvloc_temp=data.behavior.pos_rews{:}(:,1);
        pos=[pos;data.behavior.pos_rews{:}];
                 
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

        if sess==1
            sess_i=ones(size(hv_v));
        else
            sess_i=[sess_i;ones(size(sess_i)).*2];
        end
    end
    
    %now we have behavioural and gaze data for one partic, both sessions
    % next : go through trials
    hv_lv=zeros(540,1);   hv_dv=zeros(540,1);   hv_all=zeros(540,1);

    lv_hv=zeros(540,1);  lv_dv=zeros(540,1);   lv_all=zeros(540,1);

    dv_hv=zeros(540,1);    dv_lv=zeros(540,1);   dv_all=zeros(540,1);
    
    all_dv=zeros(540,1);  all_hv=zeros(540,1);    all_lv=zeros(540,1);
    
    %bidirectional
    bi_dv_hv=zeros(540,1); bi_dv_lv=zeros(540,1);  bi_hv_lv=zeros(540,1);

    for trial=1:540
        if ~isempty(gaze{trial})
            for i=2:size(gaze{trial},1)
                if gaze{trial}(i,1)==1
                    hv_all(trial)=hv_all(trial)+1;

                    if  gaze{trial}(i,2)==2
                        hv_lv(trial)=hv_lv(trial)+1;
                        bi_hv_lv(trial)=bi_hv_lv(trial)+1;
                    elseif gaze{trial}(i,2)==3
                        hv_dv(trial)=hv_dv(trial)+1;
                        bi_dv_hv(trial)=bi_dv_hv(trial)+1;
                    end

                elseif gaze{trial}(i,1)==2
                    lv_all(trial)=lv_all(trial)+1;

                     if  gaze{trial}(i,2)==1
                        lv_hv(trial)=lv_hv(trial)+1;
                        bi_hv_lv(trial)=bi_hv_lv(trial)+1;
                    elseif gaze{trial}(i,2)==3
                        lv_dv(trial)=lv_dv(trial)+1;
                        bi_dv_lv(trial)=bi_dv_lv(trial)+1;
                    end


                elseif gaze{trial}(i,1)==3
                    dv_all(trial)=dv_all(trial)+1;

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
                elseif gaze{trial}(i,2)==3
                    all_dv(trial)=all_dv(trial)+1;
                elseif gaze{trial}(i,2)==2
                    all_lv(trial)=all_lv(trial)+1;
                end
                    
            end
        end
    end

    hv_e=nr_fix(:,1);
    lv_e=nr_fix(:,2);
    d_e=nr_fix(:,3);
    x_e=nr_fix(:,4);
    all_e=sum(nr_fix,2);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% GLM

    for run=1:length(Conds)
        if length(Conds)==1 %all non-tms
            idx=tms==0;
        else
           % sort into conds
            if Conds{run}(2)=='I' %MIP
                this_session=2;
                tms_ind=4;
            elseif Conds{run}(2)=='T'%MT
                this_session=1;
                tms_ind=3;
            end
            if Conds{run}(tms_ind)=='1' %TMS
                this_tms=0;
            elseif Conds{run}(tms_ind)=='0' %Non-TMS
                this_tms=1;
            end
            if Conds{run}(end-1:end)=='si' %ipsilateral
                d_loc_oi=0;
                hv_loc_oi=0;
            elseif Conds{run}(end-1:end)=='ra' % contralateral
                d_loc_oi=1;
                hv_lov_oi=0;
            end

            if ANOVA==1
                idx= sess_i==this_session & tms==this_tms & ~isnan(acc);
            elseif ANOVA==2 
                idx= sess_i==this_session & tms==this_tms & dloc==d_loc_oi & ~isnan(acc) ;
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

            else
                criterion=eval(criterion_name);
                criterion=[criterion(idx)];
                [betas(partic,:,run),dev,stats]=glmfit(regressors,criterion);

            end
        end  
    end
end
exclude=unique(exclude);
fprintf('%i participants excluded \n',length(exclude))
betas(exclude,:,:)=[];


 %% test each regressor against 0
fprintf('\n ---------------------\n')
fprintf('  t-tests')
fprintf('\n ---------------------\n')
if ANOVA < 10
    % only non-tms conds
    conds_oi=[];
    for cond=1:length(Conds)
        if any(Conds{cond}=='0');
            conds_oi=[conds_oi, cond];
        end
    end
    for regressor=2:size(betas,2)  
        betas_oi=squeeze(betas(:,regressor,[conds_oi]));
        betas_oi=mean(betas_oi,2);
        [h p ci stats]=ttest(betas_oi);
        D=(mean(betas_oi)-0)/std(betas_oi);
        ts = tinv([0.025 0.975],length(betas_oi)-1);
        CI = mean(betas_oi) + ts.*(std(betas_oi)/sqrt(length(betas_oi)));
        fprintf('\nTesting %s against 0:\n',Regs{regressor-1})
        fprintf(' t(%d) = %2.2f, p = %2.3f, d = %2.2f, CI = [%2.2f, %2.2f]\n',stats.df,stats.tstat, p,D,CI)

    end
elseif ANOVA==10
    for regressor=2:size(betas,2)  
        betas_oi = betas(:,regressor);
        [binary,p,ci,tstat]=ttest(betas_oi);
        D=(mean(betas_oi)-0)/std(betas_oi);
        ts = tinv([0.025 0.975],length(betas_oi)-1);
        CI = mean(betas_oi) + ts.*(std(betas_oi)/sqrt(length(betas_oi)));
        fprintf('\nTesting %s against 0:\n',Regs{regressor-1})
        fprintf('t(%i) = %2.2f, p = %1.3f, d = %2.2f, CI = [%2.2f, %2.2f]\n', tstat.df, tstat.tstat, p,D,CI)
    end
end


%% plot betas
for c=1:length(Conds)
    for r=1:length(Regs)+1
        x=(betas(~isnan(betas(:,r ,c)),r ,c));
        SEM = std(x)/sqrt(length(x)); % Standard Error
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
    cd(output_dir)
    print -depsc all_gazes
else
    title(strcat('Nr Fixations (mean/trial)'))
end

if ANOVA<10
    fprintf('\n ---------------------\n')
    fprintf('  ANOVAs')
    fprintf('\n ---------------------\n')
    %% anovas
    if length(Conds)>1
        result_text=[];
        for regressor=2:length(Regs)+1
            titles=Conds;
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
            ranovatbl = ranova(rm,'withinmodel',factors);
            if ANOVA==2
                p_spots=[3 5 7 9 11 13 15];
                p_spot_names={'Session' 'Stim' 'DLoc' 'Session*Stim' 'Session*DLoc*' 'Stim*DLoc' 'Session*Stim*DLoc'};
            elseif ANOVA==1
                p_spots=[3 5 7];
                p_spot_names={'Session' 'Stim' 'Session*Stim*'};
            end
            fprintf('\n- - %s - -\n',Regs{regressor-1})
            result_text=sprintf('%s \n- - %s - -\n'  , result_text,Regs{regressor-1});
            ranova(rm,'withinmodel',factors);
            ranovatbl=table2array(ranovatbl);
            %report
            for i=1:length(ranovatbl([p_spots],5))
                fprintf('%s : F(%d,%d) = %4.2f, p = %1.3f, n = %4.2f \n',p_spot_names{i},ranovatbl(1,2),ranovatbl(2,2), ranovatbl ((p_spots(i)),4),ranovatbl((p_spots(i)),5), (ranovatbl((p_spots(i)),1)/(ranovatbl((p_spots(i)),1) + ranovatbl((p_spots(i))+1,1))))
                result_text=sprintf('%s %s : F(%d,%d) = %4.2f, p = %1.3f \n'  , result_text,p_spot_names{find(p_spots(i))},  ranovatbl(1,2),ranovatbl(2,2), ranovatbl ((p_spots(i)),4),ranovatbl((p_spots(i)),5));
            end   
        end
    end

    
    %% plot Figure 5
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

        legend(Conds)
        fig_name=strcat('eye',Regs{rfig-1},'_bar');
        if criterion_name(1)=='b'
                ylim([-.04 .1])
                if Which_GLM==5
                    ylim([-.11 .07])
                end
                set(gca,'YTick',[-.1 -.05 0 .05 .1])
        elseif criterion_name(1)=='d'
                ylim([-.04 .07])
                if Which_GLM==5
                    ylim([-.07 .05])
                end
                set(gca,'YTick',[-.04 -.02 0 .02 .04 .06])
                fig_name=strcat(fig_name,'_',criterion_name(1),'_dvhv');
        elseif criterion_name(1)=='h'
                ylim([-.02 .04])
                if Which_GLM==5
                    ylim([-.05 .025])
                end
                set(gca,'YTick',[-.02  0 .02 .04])
                fig_name=strcat(fig_name,'_',criterion_name(1),'_hvdv');
        end
        set(gca,'xticklabel',{'MT' 'MIP'})
        ticks=get(gca,'YTick');
        title(Regs{rfig-1})

        %save
        cd(output_dir)
        print(fig_name, '-depsc') 
    end
    

    %% Post Hoc
    %fix session to MT (then MIP) to see if stim matters
    fprintf('\n ---------------------\n')
    fprintf('  Follow up t-test')
    fprintf('\n ---------------------\n')

    % only MIP
    tabbi=table2array(tab);
    tabbi=tabbi(:,1:2);
    % ttest between tms and non-tms
    [h p ci tstat]=ttest(tabbi(:,1),tabbi(:,2));
    Mean1= mean(tabbi(:,1));
    Mean2= mean(tabbi(:,2));
    SD1 = std(tabbi(:,1));
    SD2 = std(tabbi(:,2));
    D= (Mean2-Mean1)/(sqrt(((SD1)^2 +(SD2)^2)/2));
    ts = tinv([0.025 0.975],length(tabbi(:,1))-1);
    CI = (Mean2-Mean1) + ts.*(std(tabbi(:,1)-tabbi(:,2))/sqrt(length(tabbi(:,1))));
    fprintf('\nOnly MIP: compare tms and non-tms \n')
    fprintf('t(%i) = %2.2f, p = %1.3f, d = %2.2f, CI = [%2.2f, %2.2f]\n', tstat.df, tstat.tstat, p,D,CI)

    
    % only MT
    tabbi=table2array(tab);
    tabbi=tabbi(:,3:4);
    % ttest between tms and non-tms
    [h p ci]=ttest(tabbi(:,1),tabbi(:,2));
    % ttest between tms and non-tms
    [h p ci tstat]=ttest(tabbi(:,1),tabbi(:,2));
    Mean1= mean(tabbi(:,1));
    Mean2= mean(tabbi(:,2));
    SD1 = std(tabbi(:,1));
    SD2 = std(tabbi(:,2));
    D= (Mean2-Mean1)/(sqrt(((SD1)^2 +(SD2)^2)/2));
    ts = tinv([0.025 0.975],length(tabbi(:,1))-1);
    CI = (Mean2-Mean1) + ts.*(std(tabbi(:,1)-tabbi(:,2))/sqrt(length(tabbi(:,1))));
    fprintf('\nOnly MT: compare tms and non-tms \n')
    fprintf('t(%i) = %2.2f, p = %1.3f, d = %2.2f, CI = [%2.2f, %2.2f]\n', tstat.df, tstat.tstat, p,D,CI)
end