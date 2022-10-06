
%Note that when calculating the median number of days surveillance is
%required for different confidence levels, the outbreak simulation needs to
%be run for large times (e.g. >50), and Sdays>25.

UnifStart=0; %Set to 1 if sampling start time is uniform between 1 and 14 days
Sdays=13;
ReqProb=0.99; %The required confidence of detecting FMDV, used to calculate how long sampling lasts (Table 3)

%% Calculate prob from single sample

var1=load('Data/EnvSampProb');
samp=randi(10000,Ntrials,1);
bsamp=var1.bsamp(samp,:);
[~, ~, ~, Dprob]=Detection(E,Cs,Ni,TotInf,1,1,1,bsamp,N,dt);

for tau=1:size(Dprob,2)
    for k=1:size(Dprob,3)
        Etype=randi([2 5],1,1);
        SingleSamp(tau,k)=Dprob(Etype,tau,k);
    end
end
plot(dt:dt:Ndays,nanmedian(SingleSamp,2),'linewidth',3)


 hold on
    patch([dt:dt:Ndays, flip(dt:dt:Ndays,2)],...
    [prctile(SingleSamp',2.5), flip(prctile(SingleSamp',97.5),2)],...
    'b','linestyle','none','FaceVertexAlpha',0.5,...
    'FaceAlpha',0.3);

ax = gca;
ax.FontSize = 16;

ylabel('Probability of detection','interpreter','latex','FontSize',22);
xlabel('Day since first infection','interpreter','latex','FontSize',22);

set(1,'paperunits','centimeters');
set(1,'papersize',[20 10]);
set(1,'paperposition',[0 0 20 10]);
% fig1='SingleSampleProbability';
% print('-dtiff','-r300',[fig1,'.tif'])
% print(1,'-dpdf',[fig1,'.pdf']);
% savefig([fig1,'.fig']);

%% Calculate likely infection times
% Estimates the probability of the time between a premises being infected
% by another and the time that FMDV is detected at the previous premises

%Calculate herd generation time:
xp=0:dt:Ndays;
Tgd=trapz(xp,TotInf.*xp') ./ trapz(xp,TotInf) ;
Tge=trapz(xp,sum(E,2).*xp') ./ trapz(xp,sum(E,2)) ;
GTd=datasample(Tgd,Ntrials);
GTe=permute(datasample(Tge,Ntrials),[1 3 2]);
 
% The detection times estimated from 2001 epidemic:
 PrevDetMean=8.07; PrevDetVar=6.67;
 
 % Calculate the time sampling would start after infection (swap GTd for GTe if 
 % environmental contamination is used to estimate infectiousness rather than direct contact)
 SampleStart=gamrnd(PrevDetMean^2/PrevDetVar,PrevDetVar/PrevDetMean,1,Ntrials)-GTd;
 

figure
histogram(SampleStart(SampleStart>=0),0:11,'Normalization','probability')
% histogram(SampleStart,-28:10,'Normalization','probability')
set(1,'paperunits','centimeters');
set(1,'papersize',[18 14]);
set(1,'paperposition',[0 0 18 14]);
ax = gca;
ax.FontSize = 20;
ylabel('Probability','interpreter','latex','FontSize',28);
xlabel('Day since first infection','interpreter','latex','FontSize',28);
% fig1='SampleStartDayPos_Tge'
% print('-dtiff','-r300',[fig1,'.tif'])
% print(1,'-dpdf',[fig1,'.pdf']);
% savefig([fig1,'.fig']);

%% Probability of detection for random samples with an unknown inf. time 
%(testing starts at SampleStart days after first infection and ends before Ndays)
clear cumPr PrNeg
  
if UnifStart==1
    SampleStart=1:14;
end

Nsarr=[1 3 5 10 20];
Titles=["Daily sampling", "2 day interval", "3 days interval"];

figure

for j=1:3
    
    Sint=j;
    cumPr=zeros(ceil(Sdays/Sint),Ntrials,size(Nsarr,2));
    PrNeg=zeros(ceil(Sdays/Sint),Ntrials);
    
    SSt=ceil(datasample(SampleStart(SampleStart>=0),Ntrials)/dt);

    for s=1:size(Nsarr,2)
        Nsamp=Nsarr(s);
    for k=1:Ntrials
        Sampletimes=(0:Sint/dt:(Sdays-1)/dt)+SSt(k);
        Etypes=randi([2 5],Sdays,Nsamp);
        for i=1:ceil(Sdays/Sint)
            PrNeg(i,k)=prod(1-Dprob(Etypes(i,:),Sampletimes(i),k));
        end
        if max(cumprod(PrNeg(:,k))<1-ReqProb)>0
        DaysNeeded(k,s)=find(cumprod(PrNeg(:,k))<1-ReqProb,1)*Sint-(Sint-1);
        else DaysNeeded(k,s)=nan;
        end
        cumPr(:,k,s)=cumprod(PrNeg(:,k));

    end

    subplot(2,3,j)
    plot(1:Sint:Sdays,nanmedian(1-cumPr(:,:,s),2),'linewidth',2)
    hold on

    
    subplot(2,3,j+3)
    plot(1:Sint:Sdays,nanmedian(cumPr(:,:,s),2),'linewidth',2)
    hold on


end
plot([1 Sdays],[0.01 0.01],'--','color','black','linewidth',1.5)
plot([1 Sdays],[0.001 0.001],'--','color','black','linewidth',1.5)
plot([1 Sdays],[0.0001 0.0001],'--','color','black','linewidth',1.5)
set(gca,'yScale','log')
yticks([10^-5, 10^-4, 10^-3, 10^-2, 0.1, 1])
ylim([10^-5 1])
xlim([1 Sdays])
ax = gca;
ax.FontSize = 20;
if j==1
ylabel('Pr negative sample','interpreter','latex','FontSize',24);
end
xlabel('Day of sampling','interpreter','latex','FontSize',24);



subplot(2,3,j)
xlim([1 Sdays])
ax = gca;
ax.FontSize = 20;
title(Titles(j),'interpreter','latex','FontSize',24);
if j==1
ylabel('Pr positive sample','interpreter','latex','FontSize',24);
hl=legend({'$s=1$', '$s=3$', '$s=5$', '$s=10$', '$s=20$'},'Location','southeast');
set(hl,'interpreter','latex','FontSize',16);
end


DN=round(nanmedian(DaysNeeded),3,'significant');
DN5=round(prctile(DaysNeeded,5),3,'significant');
DN95=round(prctile(DaysNeeded,95),3,'significant');

% Display the table data so each entry can be copied into Latex
disp(['------- ',num2str(j),' day interval -------'])
for i=1:size(DN,2)
    disp([num2str(Nsarr(i)),' samples: ',num2str(DN(i)),' \newline(',num2str(DN5(i)),'-',num2str(DN95(i)),') & '])
end



end

set(1,'paperunits','centimeters');
set(1,'papersize',[32 18]);
set(1,'paperposition',[-1 0 34 18]);
% set(1,'papersize',[18 14]);
% set(1,'paperposition',[0 0 18 14]);
ax = gca;
ax.FontSize = 20;
% fig1=['PrCumED_Tge'];
% fig1=['PrNegativeCumED_',num2str(Sint),'dayint'];
% print('-dtiff','-r300',[fig1,'.tif'])
% print(1,'-dpdf',[fig1,'.pdf']);
% savefig([fig1,'.fig']);



