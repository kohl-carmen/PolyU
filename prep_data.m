%prep data for github
dir = fileparts(which('prep_data.m'));
cd(dir);
Partic = 1:31;

for partic = Partic
    partic_str=sprintf('%02d', partic);
    MT=load(strcat('Data_okd\',partic_str,'0_BehaviouralData'));
    MT=MT.data.behavior;
    new=nan(270,18);
    new(:,1)=MT.trial{:};
    new(:,2)=MT.HV{:};
    new(:,3)=MT.LV{:};
    new(:,4)=MT.D{:};
    new(:,5)=MT.rews{:}(:,1);
    new(:,6)=MT.rews{:}(:,2);
    new(:,7)=MT.rews{:}(:,3);
    new(:,8)=MT.probs{:}(:,1);
    new(:,9)=MT.probs{:}(:,2);
    new(:,10)=MT.probs{:}(:,3);
    new(:,11)=MT.pos_rews{:}(:,1);
    new(:,12)=MT.pos_rews{:}(:,2);
    new(:,13)=MT.pos_rews{:}(:,3);
    new(:,14)=MT.tms{:};
    new(:,15)=MT.choice{:};
    new(:,16)=MT.RT{:};
    new(:,17)=MT.outcome{:};
    new(:,18)=MT.accuracy{:};
    MT=new;
    
    
    MIP=load(strcat('Data_okd\',partic_str,'1_BehaviouralData'));
    MIP=MIP.data.behavior;
    new=nan(270,18);
    new(:,1)=MIP.trial{:};
    new(:,2)=MIP.HV{:};
    new(:,3)=MIP.LV{:};
    new(:,4)=MIP.D{:};
    new(:,5)=MIP.rews{:}(:,1);
    new(:,6)=MIP.rews{:}(:,2);
    new(:,7)=MIP.rews{:}(:,3);
    new(:,8)=MIP.probs{:}(:,1);
    new(:,9)=MIP.probs{:}(:,2);
    new(:,10)=MIP.probs{:}(:,3);
    new(:,11)=MIP.pos_rews{:}(:,1);
    new(:,12)=MIP.pos_rews{:}(:,2);
    new(:,13)=MIP.pos_rews{:}(:,3);
    new(:,14)=MIP.tms{:};
    new(:,15)=MIP.choice{:};
    new(:,16)=MIP.RT{:};
    new(:,17)=MIP.outcome{:};
    new(:,18)=MIP.accuracy{:};
    MIP=new;
    
    Key={'Trial_Nr','HV_value','LV_value','D_value','HV_magnitude','LV_magnitude','D_magnitude','HV_probability','LV_probability','D_probability',     'HV_position','LV_position','D_position','TMS','Decision','RT','Reward','Accuracy'};
    
    data=struct();
    data.MIP=MIP;
    data.MT=MT;
    data.Key=Key;
    
    save(partic_str,'data')
end



    