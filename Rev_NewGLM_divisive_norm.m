%% Supplementary Figure S4
% 	- applies GLM ({'HV+LV' 'HV-LV' 'D' '(HV+LV)(D)'}) to behavioural data (tms x session x side)
% 	- reports GLM effects (tms x session x side ANOVA)
% 	- plots Supplementary Figure S4


% run that GLM per person, per session(MIP/MT) and by stim cond(tms/non)
% run Session(MIP/MT) x (TMS(stim/non) ANOVA on each regressor ('HV+LV' 'HV-LV' 'D' 'D * HV+LV')

clearvars 

beh_data_folder='D:\PolyU\TMS\Data\Beh';
plot_dir = 'C:\Users\ckohl\Desktop';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Partic=1:31;
tms_dev_cutoff=5;
rt_min=0;%100;
rt_max=1500;%1499;
exclude=[];

%% GLM
Conds={'MIP0ipsi' 'MIP1ipsi' 'MT0ipsi' 'MT1ipsi' 'MIP0contra' 'MIP1contra' 'MT0contra' 'MT1contra'};

Regs={'HV+LV' 'HV-LV' 'D' '(HV+LV)(D)'};
model='[hv_v(idx)+lv_v(idx), hv_v(idx)-lv_v(idx), d_v(idx),normalise(hv_v(idx)+lv_v(idx)).*normalise(d_v(idx))]';
 
ACCURACY = [];
HV_V = [];
D_V = [];
LV_V = [];
POSITION = [];
SESSION = [];
TMS = [];
DELETE = [];

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
      
    DELETE(:,partic)=(tms_deletion | accuracy_deletion | rt_keep<=rt_min | rt_keep>=rt_max);
    
    clear acc data file MIP V5 temp
    
       
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

            [betas(partic,:,run),dev,stats]=glmfit(regressors,criterion,'binomial');

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

betas(exclude,:,:)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANOVA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for regressor=2:length(Regs)+1

    titles=Conds;
    tab=array2table(squeeze(betas(:,regressor,:)),'variablenames',titles);
    
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
    %this is to only report significant ones
    if sum(ranovatbl([p_spots],5)<.1)>0
       sig_idx=find(ranovatbl([p_spots],5)<.1);
       sig_idx=p_spots(sig_idx);
       for i=1:length(sig_idx)
           fprintf('%s : F(%d,%d) = %4.2f, p = %1.3f, n = %4.3f \n',p_spot_names{find(p_spots==sig_idx(i))},  ranovatbl(1,2),ranovatbl(2,2), ranovatbl ((sig_idx(i)),4),ranovatbl((sig_idx(i)),5), (ranovatbl((sig_idx(i)),1)/(ranovatbl((sig_idx(i)),1) + ranovatbl((sig_idx(i))+1,1))))
       end   
    else
        disp('nothing <.1')
    end
%     sig_idx=[p_spots];
%     for i=1:length(sig_idx)
%        fprintf('%s : F(%d,%d) = %4.2f, p = %1.3f, n = %4.3f \n',p_spot_names{find(p_spots==sig_idx(i))},  ranovatbl(1,2),ranovatbl(2,2), ranovatbl ((sig_idx(i)),4),ranovatbl((sig_idx(i)),5), (ranovatbl((sig_idx(i)),1)/(ranovatbl((sig_idx(i)),1) + ranovatbl((sig_idx(i))+1,1))))
%     end   

end
                    


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% TTESTS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % test each regressor against 0
% mean_betas=squeeze(mean(betas,3));
% p=[];
% t=[];
% D=[];
% fprintf('Ttest of mean betas against 0 (across partics and conds):\n')
% for regressor=2:size(mean_betas,2)
%     betas_oi = mean_betas(:,regressor);
%     [h p ci stats]=ttest(betas_oi);
%     D=(mean(betas_oi)-0)/std(betas_oi);
%     ts = tinv([0.025 0.975],length(betas_oi)-1);
%     CI = mean(betas_oi) + ts.*(std(betas_oi)/sqrt(length(betas_oi)));
%     fprintf('\nTesting %s against 0:\n',Regs{regressor-1})
%     fprintf('- all nontms: t(%d) = %2.2f, p = %2.3f, d = %2.2f, CI = [%2.2f, %2.2f]\n',stats.df,stats.tstat, p,D,CI)
% end
% fprintf('\nNow per cond')
% for regressor=2:size(mean_betas,2)
%     fprintf('\nTesting %s against 0:\n',Regs{regressor-1})
%     for cond = 1:length(Conds)
%         betas_oi = betas(:,regressor,cond);
%         [h p ci stats]=ttest(betas_oi);
%         D=(mean(betas_oi)-0)/std(betas_oi);
%         ts = tinv([0.025 0.975],length(betas_oi)-1);
%         CI = mean(betas_oi) + ts.*(std(betas_oi)/sqrt(length(betas_oi)));
%         fprintf('- %s: t(%d) = %2.2f, p = %2.3f, d = %2.2f, CI = [%2.2f, %2.2f]\n',Conds{cond},stats.df,stats.tstat, p,D,CI)
%     end
% end
%         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT BETAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 3 a-c
% errors
for c=1:length(Conds)
    for r=1:length(Regs)+1
        x=(betas(~isnan(betas(:,r ,c)),r ,c));
        SEM = std(x)/sqrt(length(x));               % Standard Error
        errors.(strcat('Reg',num2str(r))).(Conds{c})=[SEM];
    end
end

for rfig=2: length(Regs)+1 %one figure per regressor (except intercept)   
    figure(rfig-1)
    to_plot=[];

    clf
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
%     legend(leg,temp)
   

    set(gca,'XTick',[1:4])
    set(gca,'xticklabel',{'Ipsi' 'Contra' 'Ipsi' 'Contra'})
    ticks=get(gca,'YTick');
    title(Regs{rfig-1})
    
   % save
%    cd(plot_dir)
%    if rfig==2
%        print -depsc HV+LV_line
%    elseif rfig==3
%        print -depsc HV-LV_line
%    elseif rfig==4
%        print -depsc D_line
%    elseif rfig==5
%        print -depsc DxHV+LV_line
%    end
end   
     




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GITHUB VERSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - same code but runs with github data

% % Kohl C, Wong MXM, Rushworth MFS & Chau BKH: Intraparietal stimulation
% % disrupts negative distractor effects in hyman multi-alternative
% % decision-making
% 
% %% Behavioural Analysis: Supplementary Analysis S4
% % Script to generate results of Supplementary Figure S4
% % Written by Carmen Kohl, 2020.
% % github.com/kohl-carmen/MIP-TMS
% 
% % Applies GLM to predict accuracy in each condition 
% %       GLM: β0 + β1 z(HV-LV) + β2 z(HV+LV) + β3 z(D) + β4 z(D*(HV+LV)) + ε
% %       Conditions: MIP/MT x TMS/NonTMS x Ipsilateral/Contralateral D
% % Reports Session(MIP/MT) x Stimulation (TMS/NonTMS) x Distractor Location 
% % (Ipsilateral/Contralateral)ANOVAs for each GLM regressor, reports ttests
% % of each regressor's beta coefficients & plots Supplementary Figure S4
% 
% clearvars 
% % set directory
% dir = fileparts(which('GLM1.m'));
% cd(dir)
% 
% nPartic = 31; % nr of participants
% % define conditions 
% Conditions = {'MIP0ipsi', 'MIP1ipsi', 'MT0ipsi', 'MT1ipsi',...
%               'MIP0contra', 'MIP1contra', 'MT0contra', 'MT1contra'};
% CondSession = [1, 1, 0, 0, 1, 1, 0, 0];
% CondTMS = [0, 1, 0, 1, 0, 1, 0, 1];
% CondDLoc = [0, 0, 0, 0, 1, 1, 1, 1];
% % GLM
% regressor_str = {'HV+LV' 'HV-LV' 'D' '(HV+LV)(D)'};
% model='[hv(idx)+lv(idx), hv(idx)-lv(idx), d(idx),normalise(hv(idx)+lv(idx)).*normalise(d(idx))]';
% 
% beta = nan(nPartic,length(regressor_str) + 1,length(Conditions));
% for iPartic = 1:nPartic
%     %% Prepare data
%     partic_str = sprintf('%02d', iPartic);
%     load(strcat('Data\',partic_str))  
%     % select variables of interest (see data.Key)
%     session = [ones(size(data.MIP, 1),1); zeros(size(data.MT, 1),1)];
%     tms = [data.MIP(:,14); data.MT(:,14)]; % 1=TMS, 0=NonTMS 
%     d = [data.MIP(:, 4); data.MT(:, 4)]; % distractor value
%     lv = [data.MIP(:, 3); data.MT(:, 3)]; % low value
%     hv = [data.MIP(:, 2); data.MT(:, 2)]; % high value
%     d_loc = [data.MIP(:, 13); data.MT(:, 13)]; %distractor location
%     d_loc_binary = zeros(size(d_loc));
%     d_loc_binary(d_loc==2 | d_loc==4) = 1; % 0=contralateral, 1=ipsilateral
%     accuracy = [data.MIP(:, 18); data.MT(:, 18)]; % 1=high value chosen, 
%                                                   % 0=low value chosen, 
%                                                   % nan=distractor/empty 
%                                                   % quadrant chosen
%    
%     % exclude trials in which the distractor/empty quadrant was chosen
%     rmv = (isnan(accuracy));   
%     session(rmv) = [];
%     tms(rmv) = [];        
%     accuracy(rmv) = [];                
%     d(rmv) = [];        
%     lv(rmv) = [];        
%     hv(rmv) = [];     
%     d_loc_binary(rmv) = [];
%        
%     %% GLM
%     for iCondition = 1:length(Conditions)
%         % select all trials in the condition of interest 
%         idx = (session==CondSession(iCondition) & ...
%                tms==CondTMS(iCondition) & ...
%                d_loc_binary==CondDLoc(iCondition));
%         % define regressors and criterion
%         regressors = eval(model);
%         regressors = normalise(regressors);
%         criterion = accuracy(idx);
%         % fit GLM
%         beta(iPartic, :, iCondition) = glmfit(regressors,...
%                                              criterion, 'binomial');
%     end
% end
% 
% 
% %% ANOVA
% % Site x Stimulation x Distractor Location ANOVA for each regressor
% fprintf('\nANOVA (GLM):')
% iEffect = 3:2:15;
% % set up ANOVA
% Effect = {'Session', 'Stim', 'DLoc', 'Session*Stim', 'Session*DLoc*', ...
%         'Stim*DLoc', 'Session*Stim*DLoc'};
% within1 = categorical(CondSession)';
% within2 = categorical(CondTMS)';
% within3 = categorical(CondDLoc)';
% within = table(within1, within2, within3, 'variablenames',...
%                {'Session', 'Stim', 'DLoc'});
% factors = 'Session*Stim*DLoc';
% allvars = strcat(Conditions{1}, '-', Conditions{end}, ' ~1');
% for regressor = 2:length(regressor_str) + 1
%     % ANOVA
%     tab = array2table(squeeze(beta(:, regressor,:)),...
%                     'variablenames', Conditions);   
%     rm = fitrm(tab, allvars, 'WithinDesign', within);
%     ranovatbl = ranova(rm, 'withinmodel', factors);   
%     % print output (F-statistic, p-value, partial eta squared)
%     fprintf('\n- - %s - -\n', regressor_str{regressor-1})
%     ranovatbl = table2array(ranovatbl);
%     for i = 1:length(iEffect)
%          fprintf('%s : F(%d, %d) = %4.2f, p = %1.3f, n2 = %4.3f \n',...
%                  Effect{i},  ranovatbl(1, 2), ranovatbl(2, 2), ...
%                  ranovatbl (iEffect(i), 4), ranovatbl(iEffect(i), 5), ...
%                  ranovatbl(iEffect(i), 1) / (ranovatbl(iEffect(i), 1) + ...
%                  ranovatbl(iEffect(i) + 1, 1))) 
%     end
% end         
% 
%     
% %% Plot Figure S4
% xax = [1, 2; 3, 4];
% clr.s1 = {[ 1, 0.75, 0.5],[0.875, 0.438, 0]};
% clr.s0 = {[ 0.5, 0.5, 0.5], 'k'};
% linetype = {'--', '-'};
% for rfig = 2: length(regressor_str) + 1
%     figure(rfig-1)
%     hold on
%     count = 0;
%     leg = nan(1, length(unique(session))*length(unique(tms)));
%     for sessioni = 0:1
%         for tmsi = 0:1
%             count = count + 1;
%             conds_oi = find(CondSession==sessioni & CondTMS==tmsi);
%             to_plot = squeeze(mean(beta(:, rfig, conds_oi)));
%             leg(count) = plot(xax(sessioni + 1,:), to_plot, ...
%                          linetype{tmsi + 1}, 'Linewidth', 2, 'Color', ...
%                          clr.(strcat('s', num2str(sessioni))){tmsi + 1});
%             %standard error
%             x = beta(:, rfig, conds_oi(1));
%             errors(1) = std(x) / sqrt(length(x));       
%             x = beta(:, rfig, conds_oi(2));
%             errors(2) = std(x) / sqrt(length(x));          
%             errorbar(xax(sessioni + 1, :), to_plot, errors, ...
%                      'Linewidth',2, 'Linestyle', 'none', 'Color', ...
%                      clr.(strcat('s',num2str(sessioni))){tmsi + 1});
%         end
%     end
%     legend(leg,{'MT0', 'MT1', 'MIP0', 'MIP1'})
%     xlim([0 5])
%     set(gca,'XTick',1:4)
%     set(gca,'xticklabel', {'Contra' 'Ipsi' 'Contra' 'Ipsi'})
%     title(regressor_str{rfig-1})        
%     if rfig <4, ylim([-0.6 1.2]); else, ylim([-0.6 0.3]); end
% end
%      
% 
% 
% function x = normalise(x)
%     % removes mean from data x
%     dim = 1;
%     dims = size(x);
%     dimsize = size(x,dim);
%     dimrep = ones(1, length(dims));
%     dimrep(dim) = dimsize;
% 
%     x = (x - repmat(nanmean(x, dim), dimrep)) ./ ...
%          repmat(nanstd(x, 0, dim), dimrep);
% end
% 
% 

