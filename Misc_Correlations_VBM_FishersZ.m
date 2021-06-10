%% VBM Correlations
% - runs correlations between  MT/MIP grey matter and MT/MIP effect
% - plots Figure 4 b, c


%   MIP/MT GM were extracted from scans (see VBM Notes Google Docs). MIP was
%   based on where we found a signifcant effect. MT was based on TMS
%   coordinates.

%   The TMS effect is taken from the behavioural analysis. The Predictor
%   D-HV, but non-tms subtracted from tms. So MIP_1_contra - MIP_0_contra
%   (same for MT). You can recreate what's loaded here as a design (e.g.
%   D_MIP) by running GLM3_Beh.m and subtracted contralateral non-tms from 
%   tms (e.g. MIP1contra-MIP0contra).


%Sources
% https://www.youtube.com/watch?v=I0SjHVOHztc&ab_channel=onlinestatbook
% http://www.stat.yale.edu/Courses/1997-98/101/confint.htm

clear
corr_or_partialcorr=2;

file_dir = 'D:\PolyU\TMS\Results\VBM\output\';
plot_dir='C:\Users\ckohl\Desktop';
glm2_output_dir ='C:\Users\ckohl\Desktop\Current\Other\Bolton\'; %output_dir in GLM2_all_nontms.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load

%load MIP GM
filename = strcat(file_dir,'GM_MIP.txt');
delimiter = ' ';
formatSpec = '%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true,  'ReturnOnError', false);
fclose(fileID);
dataArray=dataArray{1};
MIP=[];
for i=1:length(dataArray)
    MIP(i)=str2num(dataArray{i});
end
MIP=MIP';

%load MT GM
filename = strcat(file_dir,'GM_MT.txt');
delimiter = ' ';
formatSpec = '%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true,  'ReturnOnError', false);
fclose(fileID);
dataArray=dataArray{1};
MT=[];
for i=1:length(dataArray)
    MT(i)=str2num(dataArray{i});
end
MT=MT';

%load MIP design
filename = strcat(file_dir,'MIP_Test_Design.txt');
delimiter = '\t';
formatSpec = '%*s%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
MIPtemp1 = dataArray{:, 1};
MIPtemp2 = dataArray{:, 2};
D_MIP = dataArray{:, 3};

%load MT design
filename = strcat(file_dir,'MT_Test_Design.txt');
delimiter = '\t';
formatSpec = '%*s%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
MTtemp1 = dataArray{:, 1};
MTtemp2 = dataArray{:, 2};
D_MT = dataArray{:, 3};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% correlate
% Apparently, the z’ for each correlation coefficient r is defined as z = 0.5*log((1+r)/(1-r)) and the standard error of z’ is 1/sqrt(N-3).
% 	Given that we can get the probability that we’ll find a specific r given N

%standard error for fisher's z (Depends only on N)
N=length(D_MIP);
SEZ = 1/(sqrt(N-3));
% fisher's z transofmration:
% z = (0.5)*log((1+r)/(1-r))
if corr_or_partialcorr==1
    %correlations
    [r p]=corr(MIP,D_MIP);
    z(1)=(0.5)*log((1+r)/(1-r));
    MIP_MIP=sprintf('r(%i) = %2.2f, p = %2.3f',length(MIP)-2,r,p)   
    [r p]=corr(MIP,D_MT);
    z(2)=(0.5)*log((1+r)/(1-r));
    MIP_MT=sprintf('r(%i) = %2.2f, p = %2.3f',length(MIP)-2,r,p)
    [r p]=corr(MT,D_MIP);
    z(3)=(0.5)*log((1+r)/(1-r));
    MT_MIP=sprintf('r(%i) = %2.2f, p = %2.3f',length(MIP)-2,r,p)
    [r p]=corr(MT,D_MT);
    z(4)=(0.5)*log((1+r)/(1-r));
    MT_MT=sprintf('r(%i) = %2.2f, p = %2.3f',length(MIP)-2,r,p)
elseif corr_or_partialcorr==2
    %partial correlations
    [r p]=partialcorr(MIP,D_MIP,[MIPtemp1,MIPtemp2]);
    z(1)=(0.5)*log((1+r)/(1-r));
    MIP_MIP=sprintf('r(%i) = %2.2f, p = %2.3f',length(MIP)-4,r,p)
    [r p]=partialcorr(MIP,D_MT,[MTtemp1,MTtemp2]);
    z(2)=(0.5)*log((1+r)/(1-r));
    MIP_MT=sprintf('r(%i) = %2.2f, p = %2.3f',length(MIP)-4,r,p)
    [r p]=partialcorr(MT,D_MIP,[MIPtemp1,MIPtemp2]);
    z(3)=(0.5)*log((1+r)/(1-r));
    MT_MIP=sprintf('r(%i) = %2.2f, p = %2.3f',length(MIP)-4,r,p)
    [r p]=partialcorr(MT,D_MT,[MTtemp1,MTtemp2]);
    z(4)=(0.5)*log((1+r)/(1-r));
    MT_MT=sprintf('r(%i) = %2.2f, p = %2.3f',length(MIP)-4,r,p)
end

% now we should assume same sort of actual population r which has t be the
% mean of the normal distribution:https://www.youtube.com/watch?v=I0SjHVOHztc&ab_channel=onlinestatbook
% we don't have that. but whta if I assume it's zero? hoping that our
% correlation is diff from it and the rest isnt?

figure
clf
mu = 0;%pop_distr_mean
sigma = SEZ;%pop_distr_std
x = [-1:.01:1];
y = normpdf(x,mu,sigma);
plot(x,y,'k')
hold on
c = 1.96;
lim= c*sigma;
plot([lim,lim],[0,normpdf(lim,mu,sigma)],'k')
plot([-lim,-lim],[0,normpdf(-lim,mu,sigma)],'k')
%plot points
p(1) = plot(z(1),normpdf(z(1),mu,sigma),'r.','Markersize',20);
p(2) = plot(z(2),normpdf(z(2),mu,sigma),'b.','Markersize',20);
p(3) = plot(z(3),normpdf(z(3),mu,sigma),'g.','Markersize',20);
p(4) = plot(z(4),normpdf(z(4),mu,sigma),'y.','Markersize',20);


%get p-value for given z
prob = normcdf(z,mu,sigma);
labs = {'MIP-DMIP','MIP-DMT','MT-DMIP','MT-DMT'};
labels = {};
for i = 1:4
    labels{i} = strcat(labs{i},': p = ',num2str(prob(i)))
end
legend(p,labels{:})


%exclude outlier
outlier=find(D_MIP==max(D_MIP));
fprintf('Exclude outlier (partial correlation): \n')
[r_check p_check]=partialcorr(MIP([1:outlier-1,outlier+1:end]),D_MIP([1:outlier-1,outlier+1:end]),[MIPtemp1([1:outlier-1,outlier+1:end]),MIPtemp2([1:outlier-1,outlier+1:end])]);
fprintf('r(%i) = %2.2f, p = %2.3f \n',length(MIP([1:outlier-1,outlier+1:end]))-4,r_check,p_check)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot
figure
%1: MIP, MIP_Test
subplot (2,2,1)
P=polyfit(D_MIP,MIP,1);
fit=P(1)*[-2 2.5]+P(2);
plot([-2 2.5],fit,'Color',[.8 .8 .8],'Linewidth',2.5) %plot line fitted
hold on
plot(D_MIP,MIP,'.','Markersize',20,'Color',[223/255 112/255 14/255])
xlim([-2 2.5])
set(gca,'YTick', [.1 .3 .5 .7 .9])
ylim([.1 .8])
set(gca,'XTick', [-2 -1 0 1 2])
title(MIP_MIP)

%2: MIP, MT_Test
subplot (2,2,2)
P=polyfit(D_MT,MIP,1);
fit=P(1)*[-2 2.5]+P(2);
plot([-2 2.5],fit,'Color',[.8 .8 .8],'Linewidth',2.5,'Linestyle','--') %plot line fitted
hold on
plot(D_MT,MIP,'.','Markersize',20,'Color','k')
xlim([-2 2.5])
set(gca,'YTick', [.1 .3 .5 .7 .9])
ylim([.1 .8])
set(gca,'XTick', [-2 -1 0 1 2])
title(MIP_MT)

%3: MT, MIP_Test
subplot (2,2,3)
P=polyfit(D_MIP,MT,1);
fit=P(1)*[-2 2.5]+P(2);
plot([-2 2.5],fit,'Color',[.8 .8 .8],'Linewidth',2.5,'Linestyle','--') %plot line fitted
hold on
plot(D_MIP,MT,'.','Markersize',20,'Color',[223/255 112/255 14/255])
xlim([-2 2.5])
set(gca,'YTick', [.1 .3 .5 .7 .9])
ylim([.1 .7])
set(gca,'XTick', [-2 -1 0 1 2])
title(MT_MIP)

%4: MT, MT_Test
subplot (2,2,4)
P=polyfit(D_MT,MT,1);
fit=P(1)*[-2 2.5]+P(2);
plot([-2 2.5],fit,'Color',[.8 .8 .8],'Linewidth',2.5,'Linestyle','--') %plot line fitted
hold on
plot(D_MT,MT,'.','Markersize',20,'Color','k')
xlim([-2 2.5])
set(gca,'YTick', [.1 .3 .5 .7 .9])
ylim([.1 .7])
set(gca,'XTick', [-2 -1 0 1 2])
title(MT_MT)

% save
% cd(plot_dir)
% print -depsc vbm_scatter
