%% Main Behavioural Analysis
% ( 'Beh_for_Publication.m' on drive and originally, 'BehGLM_forsimplestacceffect.m')
% 	- applies GLM3 ({'HV+LV' 'HV-LV' 'D-HV' }) to behavioural data (tms x session x side)
% 	- reports GLM effects (tms x session x side ANOVA) as well as ttests of beta weight against zero
% 	- plots Figure 3 a-d


% run that GLM per person, per session(MIP/MT) and by stim cond(tms/non)
% run Session(MIP/MT) x (TMS(stim/non) ANOVA on each regressor ('HV+LV' 'HV-LV' 'D-HV' )

% NOTE: Posthoc tests are done in SPSS. include covariates.


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

Regs={'HV+LV' 'HV-LV' 'D-HV' };
model='[hv_v(idx)+lv_v(idx), hv_v(idx)-lv_v(idx),d_v(idx)-hv_v(idx)]';

% We always want to normalise the individual terms of an interaction (and
% then everything below at normalise(regressors)

ACCURACY = [];
HV_V = [];
D_V = [];
LV_V = [];
POSITION = [];
SESSION = [];
TMS = [];
DELETE = [];

for partic=1:length(Partic)
    partic_char=num2str(Partic(partic));
    if Partic(partic)<10
        partic_char=strcat('0',num2str(Partic(partic)));
    end

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
    HV_V(:,partic)=(hv_v);
    D_V(:,partic)=d_v;
    LV_V(:,partic)=lv_v;
    POSITION(:,:,partic)=position;
    SESSION(:,partic)=session;
    TMS(:,partic)=tms;
   
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
    %%this is to only report significant ones
    % if sum(ranovatbl([p_spots],5)<.1)>0
    %    sig_idx=find(ranovatbl([p_spots],5)<.1);
    %    sig_idx=p_spots(sig_idx);
    %    for i=1:length(sig_idx)
    %        fprintf('%s : F(%d,%d) = %4.2f, p = %1.3f, n = %4.3f \n',p_spot_names{find(p_spots==sig_idx(i))},  ranovatbl(1,2),ranovatbl(2,2), ranovatbl ((sig_idx(i)),4),ranovatbl((sig_idx(i)),5), (ranovatbl((sig_idx(i)),1)/(ranovatbl((sig_idx(i)),1) + ranovatbl((sig_idx(i))+1,1))))
    %    end   
    % else
    %     disp('nothing <.1')
    % end
    sig_idx=[p_spots];
    for i=1:length(sig_idx)
       fprintf('%s : F(%d,%d) = %4.2f, p = %1.3f, n = %4.3f \n',p_spot_names{find(p_spots==sig_idx(i))},  ranovatbl(1,2),ranovatbl(2,2), ranovatbl ((sig_idx(i)),4),ranovatbl((sig_idx(i)),5), (ranovatbl((sig_idx(i)),1)/(ranovatbl((sig_idx(i)),1) + ranovatbl((sig_idx(i))+1,1))))
    end   

end
fprintf('\n For posthoc: SPSS \n')
                    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TTESTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test each regressor against 0
mean_betas=squeeze(mean(betas,3));
p=[];
t=[];
D=[];
fprintf('Ttest of mean betas against 0 (across partics and conds):\n')
for i=2:size(mean_betas,2)
    [h p(i) ci stats]=ttest(mean_betas(:,i));
    t(i)=stats.tstat;
    p(i);
    Mean1= mean(mean_betas(:,i));
    Mean2= 0;
    SD1 = std(mean_betas(:,i));
    SD2 = 0;
    D(i)= (Mean2-Mean1)/(SD1)';   
    fprintf(Regs{i-1})
    fprintf('\nt(%i) = %2.3f, p = %2.3f, d = %2.3f\n',stats.df, t(i), p(i), D(i))
end
                
% test D-HV regressor against 0 - once for MIP1contra one for MIP0contra
for cond=1:length(Conds)
    if Conds{cond}(2)=='I' & Conds{cond}(end)=='a'
        if Conds{cond}(4)=='1'
            cond_oi1=cond;
        else
            cond_oi2=cond;
        end
    end
end

p=[];
t=[];
D=[];
fprintf('Ttest of D-HV against 0 (MIP1contra vs MIP0contra):\n')
for conds = 1:2
    if conds==1
        cond=cond_oi1;
    else
        cond=cond_oi2;
    end
    [h p(i) ci stats]=ttest(betas(:,end,cond));
    t(i)=stats.tstat;
    p(i);
    Mean1= mean(betas(:,end,cond));
    Mean2= 0;
    SD1 = std(betas(:,end,cond));
    SD2 = 0;
    D(i)= (Mean2-Mean1)/(SD1)'; 
    fprintf('%s %s',Conds{cond},Regs{i-1})
    fprintf('\nt(%i) = %2.3f, p = %2.3f, d = %2.3f\n',stats.df, t(i), p(i), D(i))
end



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
    %5 MIP0con 6 MIP1con 7 MT1con 8 MT0con
    hold on
    xlim([0 5])
     
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
    set(gca,'xticklabel',{'Contra' 'Ipsi' 'Contra' 'Ipsi'})
    ticks=get(gca,'YTick');
    title(Regs{rfig-1})
    
   % save
   cd(plot_dir)
   if rfig==2
       print -depsc HV+LV_line
   elseif rfig==3
       print -depsc HV-LV_line
   elseif rfig==4
       print -depsc D-HV_line
   end
end   
     


%% Figure 3d
% Explore Interaction
figure

Comparison={'D_V-HV_V'};  
avg_window=0.25;
bin_x=[0:0.25:1-avg_window];
bin_x=[bin_x;bin_x+avg_window];
%conds to compare to
Conds={'MIP1contra' 'MIP1ipsi' 'MT1contra' 'MT1ipsi'};
position_oi=3; %ipsi/contra loc refers to distractor

Runs=length(Conds);
ACC_PER=[];
ERR_PER=[];
DV_PER=[];
DV_HV_PER = [];
keep_acc_per=[];
line=[];
 for run=1:Runs
    if Conds{run}(1:2)=='MI' %MIP
       this_session=1;
       tms_ind=4;
    else % MT
       this_session=0;
       tms_ind=3;
    end
    if Conds{run}(tms_ind)=='0' % no stim
        this_tms=0;%no stim
    else %  stim
        this_tms=1;
    end
    if Conds{run}(end-1:end) =='si'% d ipsilateral
        d_loc=[1 3];
    else%  d contralateral
        d_loc=[2 4];
    end

    acc_per=[];
    x_idx=[];
    y_idx=[];
    rt_per=[];
    X=eval(Comparison{1});

    for bin_x_count=1:length(bin_x) 
            for partic=1:length(Partic)
              %pick trials according to cond etc
              idx=(DELETE(:,partic)==0 & (POSITION(:,position_oi,partic)==d_loc(1) | POSITION(:,position_oi,partic)==d_loc(2)) & SESSION(:,partic)==this_session & TMS(:,partic)==this_tms & ~isnan(ACCURACY(:,partic)));
              % quantiles based on D value
              idx= idx & X(:,partic)>= quantile(X(idx,partic),bin_x(1,bin_x_count)) &  X(:,partic)<= quantile(X(idx,partic),bin_x(2,bin_x_count));     
              acc_per(bin_x_count,partic)=mean(ACCURACY(idx,partic));
              dv_per(bin_x_count,partic)=mean(D_V(idx,partic));
              dv_hv_per(bin_x_count,partic)=mean(D_V(idx,partic)-HV_V(idx,partic));
            end
            %collect partics
            ACC_PER(bin_x_count)=mean(acc_per(bin_x_count,:),2);
            DV_PER(bin_x_count)=mean(dv_per(bin_x_count,:),2);
            DV_HV_PER(bin_x_count)=mean(dv_hv_per(bin_x_count,:),2);

            errorx=acc_per(bin_x_count,:);
            SEM = std(errorx)/sqrt(length(errorx));             
            ERR_PER(bin_x_count) =[SEM];
    end
    keep_acc_per(:,:,run)=acc_per;

    hold on
    if Conds{run}(2)=='I' & Conds{run}(tms_ind)=='1' & Conds{run}(end)=='a'
        linewidth=4;
        linestyle='-';
        color=[223/255 112/255 14/255];
    elseif Conds{run}(2)=='I' & Conds{run}(tms_ind)=='1' & Conds{run}(end)=='i'
        linewidth=1;
        linestyle='-.';
        color=[223/255 112/255 14/255];
    elseif Conds{run}(2)=='I' & Conds{run}(tms_ind)=='0' 
        linewidth=1;
        color=[250/255 190/255 129/255];
        linestyle='--';
        if Conds{run}(end)=='i'
            linestyle=':';
        end
    elseif Conds{run}(2)=='T' & Conds{run}(tms_ind)=='1'     
        linewidth=1;
        color='k';
        linestyle='-';
        if Conds{run}(end)=='i'
            linestyle='-.';
        end
    else
        linewidth=1;
        color=[128/255 128/255 128/255];
        linestyle='--';
        if Conds{run}(end)=='i'
            linestyle=':';
        end
    end
    if length(Comparison{:})==3 % D_V
        line(run)=plot(DV_PER,ACC_PER,'Linewidth',linewidth, 'Color', color,'Linestyle',linestyle);
        errorbar(DV_PER, (ACC_PER), ERR_PER, 'linestyle', linestyle,'Linewidth',linewidth, 'Color', color,'Linestyle',linestyle);
    else % DV-HV
        line(run)=plot(DV_HV_PER,ACC_PER,'Linewidth',linewidth, 'Color', color,'Linestyle',linestyle);
        errorbar(DV_HV_PER, (ACC_PER), ERR_PER, 'linestyle', linestyle,'Linewidth',linewidth, 'Color', color,'Linestyle',linestyle);
    end

    title(strcat('Acc -', Conds{run}))
    xlabel(Comparison{1})
    if length(Comparison{:})==3 % D_V
        set(gca,'YTick',[.68 .7 .72 .74]);
        set(gca,'XTick',[100 200 300]);
    end
 end
legend(line,Conds)

% save
cd(plot_dir)
print -depsc intplot

