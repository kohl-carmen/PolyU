%% Multinomial Logistic Regression - Supplementary Figure S2
%   - runs multinomial logistic regression {'HV+LV' 'HV-LV' 'D-HV' '(HV-LV)*(D-HV)'}
% 	- reports ttests of each predictor's beta weights against zero
%   - plots Supplementary Figure S2 

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
 keep=[];
exclude=[];
% GLM1
Regs={'HV+LV' 'HV-LV' 'D-HV' '(HV-LV)*(D-HV)'};
model='[hv_v(idx)+lv_v(idx), hv_v(idx)-lv_v(idx),d_v(idx)-hv_v(idx),normalise(hv_v(idx)-lv_v(idx)).*normalise(d_v(idx)-hv_v(idx))]';
 % We always want to normalise the individual terms of an interaction (and
% then everything below at normalise(regressors)
Regs={ 'HV-LV' 'D-HV' '(HV-LV)*(D-HV)'};
model='[ hv_v(idx)-lv_v(idx),d_v(idx)-hv_v(idx),normalise(hv_v(idx)-lv_v(idx)).*normalise(d_v(idx)-hv_v(idx))]';
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
    position=[MIP.choice_pos{:};V5.choice_pos{:}];  
    dchosen=[MIP.D_chosen{:};V5.D_chosen{:}];  
   
    clear acc data file MIP V5 temp
    
    del = isnan(accuracy) & dchosen==0;
    d_v(del)=[];
    accuracy(del)=[];
    hv_v(del)=[];
    lv_v(del)=[];
    position(del)=[];
    tms(del)=[];
    dchosen(del)=[];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% GLM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    idx=( tms==0 );
    accuracy = accuracy(idx);
    dchosen = dchosen(idx);
    Choice = {};
    for i = 1:length(accuracy)
        if accuracy(i)==0
            Choice{end+1}='LV';
        elseif accuracy(i)==1
            Choice{end+1}='XHV';
        elseif isnan(accuracy(i)) & dchosen(i) ==1
            Choice{end+1}='D';
        end
    end
    Choice = categorical(Choice);       
    if sum(dchosen)<3
        partic
        exclude=[exclude,partic];
    end
    keep=[keep,sum(dchosen)];
    regressors=[eval(model)];
    regressors=[normalise(regressors)];
    criterion = Choice;
    temp=mnrfit(regressors,criterion);
    betas1(partic,:) = temp(:,2);
    betas2(partic,:) = temp(:,1);
end

cats = categories(criterion);
figure
for i = 1:2
    subplot(1,2,i)
    hold on
    if i ==1
        betas = betas1;
        title(strcat('Choosing ', cats{2}, ' in ref to choosing ', cats{3}))
    else
        betas=betas2;
        title(strcat('Choosing ', cats{1}, ' in ref to choosing ', cats{3}))
    end

betas(exclude,:)=[];
betas = betas.*-1;
if i==2
    betas(:,[2,3]) = betas(:,[2,3]).*-1;
    Regs{1} = 'LV-HV';
    Regs{2} = 'HV-D';
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% test each regressor against 0
%only non-tms
for regressor=2:size(betas,2)
    betas_oi = betas(:,regressor);
    [h p ci stats]=ttest(betas_oi);
    D=(mean(betas_oi)-0)/std(betas_oi);
    ts = tinv([0.025 0.975],length(betas_oi)-1);
    CI = mean(betas_oi) + ts.*(std(betas_oi)/sqrt(length(betas_oi)));
    fprintf('\nTesting %s against 0:\n',Regs{regressor-1})
    fprintf('- all nontms: t(%d) = %2.2f, p = %2.3f, d = %2.2f, CI = [%2.2f, %2.2f]\n',stats.df,stats.tstat, p,D,CI)
end


errors=[];
for r=1:length(Regs)+1
    x=(betas(~isnan(betas(:,r )),r ));
    SEM = std(x)/sqrt(length(x));% Standard Error
    errors(r)=[SEM];
end
%plot
% figure
bar(mean(betas(:,2:end)));
hold on
err=errorbar(1:length(Regs),mean(betas(:,2:end)),errors(2:end));
err.Color = [0 0 0];                            
err.LineStyle = 'none'; 
set(gca,'XTick',1:length(Regs))
set(gca,'xticklabel',Regs)
ylim([-0.5,1.5])
end
%% save
% cd(plot_dir)
% print -depsc Multinomial
