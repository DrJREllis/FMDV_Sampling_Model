% Matlab script to plot observed vs expected titres for the environmental
% transmission experiments. This version uses a simple phenomenological for
% viral shedding and fits to titre, not log titre.
clear all
%==========================================================================
% PREPARE THE DATA
% Load the data
varload=load('Data/EnvironmentalTransmissionData');

% Extract the data
tC=varload.tC;
tObs=varload.tObs;
tEC=varload.tEC;
chalOut=varload.chalOut;
E=varload.E; % [room, day, type, VI]
roomExpt=varload.roomExpt;
roomNum=[1 2 1 2 1 2 3 1 2 3];


% Combine the viral titres for shedding (carefully!)
VI=NaN(size(varload.N));
z=(varload.N==0 & varload.O==0);
VI(z)=0;
z=(varload.N==0 & varload.O>0);
VI(z)=varload.O(z);
z=(varload.N>0 & varload.O==0);
VI(z)=varload.N(z);
z=(varload.N>0 & varload.O>0);
VI(z)=log10(10.^varload.N(z)+10.^varload.O(z));

% Set which infected animals were in which room(s)
aR=[1 2; 3 4; 5 6; 7 8; 9 10; 11 12; 11 12; 13 14; 15 16; 15 16];

% Set when infected animals were in each room
tR=zeros(10,2);
for r=1:size(tR,1)
    if ismember(r,[1 3])==1
        tR(r,:)=[0 tEC(r)];
    elseif ismember(r,[2 4])==1
        tR(r,:)=[tEC(r-1) tEC(r)];
    elseif ismember(r,[5 8])==1
        tR(r,:)=[0 tEC(r)-1];
    elseif ismember(r,[6 9])==1
        tR(r,:)=[tEC(r-1)-1 tEC(r)];
    elseif ismember(r,[7 10])==1
        tR(r,:)=[tEC(r-1) tEC(r)];
    end
end

% Determine the number of animals and rooms
nAnim=length(tC);
nRoom=length(tEC);
%==========================================================================

%==========================================================================
% LOAD THE MCMC SAMPLES
% Load the samples
varload=load('Data/EnvTrans_SimpleModelToo4_MCMCSamples');

% Extract the shedding parameters
Vp=exp([varload.ParSamp{1}(:,1:nAnim);  varload.ParSamp{2}(:,1:nAnim)]);
tp=[varload.ParSamp{1}(:,nAnim+1:2*nAnim); ...
    varload.ParSamp{2}(:,nAnim+1:2*nAnim)];
lg=[varload.ParSamp{1}(:,2*nAnim+1:3*nAnim); ...
    varload.ParSamp{2}(:,2*nAnim+1:3*nAnim)];
ld=[varload.ParSamp{1}(:,3*nAnim+1:4*nAnim); ...
    varload.ParSamp{2}(:,3*nAnim+1:4*nAnim)];

% Extract the environmental parameters
PS=exp([varload.ParSamp{1}(:,4*nAnim+8+1:end-2);
        varload.ParSamp{2}(:,4*nAnim+8+1:end-2)]);
a=PS(:,1:4);
d=PS(:,5:8);
b=PS(:,9);
%==========================================================================

%==========================================================================
% % SET PLOTTING STUFF
% % Set the plotting weight
% bandweight=0.1;
% 
% % Specify the surface descriptors
% surfLab={'floor', 'wall', 'feed trough', 'faeces'};
% 
% % Set the plot size and locations for the environmental dynamics
% wd=0.16;
% ht=0.08;
% ax_xpos=[0.20 0.39 0.58 0.77];
% ax_ypos=0.91:-0.095:0.055;
% 
% % Set the surface descriptor locations
% xl=[20 20.5 13.2 18];


%% ========================================================================== 
% CALCULATE THE ENVIRONMENTAL DYNAMICS

% Initialise the challenge does
muEChal=zeros(length(b),nRoom);

OvP=[];

% For each room (and, hence, environmental challenge) ...
for r=1:nRoom
    display(num2str(r))

% Compute the expected amount of virus shed into the room
    t=0:max(tObs);
    VRoom=zeros(size(tp,1),length(t));
    for j=1:size(aR,2)
        tau=repmat(t-tC(aR(r,j)),size(tp,1),1);
        VpR=repmat(Vp(:,aR(r,j)),1,length(t));
        tpR=repmat(tp(:,aR(r,j)),1,length(t));
        lgR=repmat(lg(:,aR(r,j)),1,length(t));
        ldR=repmat(ld(:,aR(r,j)),1,length(t));
        V=2.*VpR./(exp(-lgR.*(tau-tpR))+exp(ldR.*(tau-tpR)));
        VRoom=VRoom+V;
    end
    VRoom(:,t<=tR(r,1) | t>tR(r,2))=0;
    
% For each sample ...
    for s=1:4

% Compute the expected viral titre in each environmental sample
        intV=cumsum(exp(repmat(d(:,s),1,length(t)).*...
                            repmat(t,size(d,1),1)).*...
                        VRoom,2);
        muE=repmat(a(:,s),1,length(t)).*...
            exp(repmat(-d(:,s),1,length(t)).*repmat(t,size(d,1),1)).*...
            intV;

% ... and, hence, the contribution to the challenge
        muEChal(:,r)=muEChal(:,r)+...
                     0.5*(muE(:,t==tEC(r))+muE(:,t==tEC(r)+1));


         
    medP=(prctile(muE,50,1))';
    medPmax(r,s)=max(medP);
    
    %OvP=[day, detected PCR/VI,  median expected viral titre, type]
    OvP=[OvP; [E(E(:,1)==r & E(:,3)==s,2) E(E(:,1)==r & E(:,3)==s,4) (medP(E(E(:,1)==r & E(:,3)==s,2)+1)) E(E(:,1)==r & E(:,3)==s,3)]];
    
    end

end
 b=zeros(10000,2,4);
 
for s=1:4
    %For each sample type, PosSamp=[median expected viral titre, positive PCR/VI]
    PosSamp=[OvP(OvP(:,4)==s,3) OvP(OvP(:,4)==s,2)>1];

    %Rearrange to give number of positive samples at each contamination level
for i=1:length(PosSamp)
    ProbPos(i,1)=max(0,PosSamp(i,1)); % Contamination level
    ProbPos(i,2)=sum(PosSamp(PosSamp(:,1)==PosSamp(i,1),2)) / sum(PosSamp(:,1)==PosSamp(i,1)); %Ratio of success
    ProbPos(i,3)=sum(PosSamp(:,1)==PosSamp(i,1)); %Times the contamination level has been sampled
end
    ProbPos = unique(ProbPos,'rows');
    
    
for k=1:2
        disp([k s])
% bsamp=zeros(10000,2,4);
b(1,k,s)=rand*0.1;
PrPS=1-exp(-b(1,k,s).*ProbPos(:,1));
LH2=prod( binopdf(ProbPos(:,2).*ProbPos(:,3),ProbPos(:,3),PrPS) );
Suc(k,s)=0; Fail(k,s)=0;
stdev=0.05;

xi=0:0.001:0.1;
for xx=1:length(xi)
    LHb(xx,s)=prod( binopdf(ProbPos(:,2).*ProbPos(:,3),ProbPos(:,3),1-exp(-xi(xx).*ProbPos(:,1))) );
end
    
while Suc(k,s)<11000
    b1=b(Suc(k,s)+1,k,s)+randn*stdev;
    
    PrPS=1-exp(-b1.*ProbPos(:,1));

    LH1=prod( binopdf(ProbPos(:,2).*ProbPos(:,3),ProbPos(:,3),PrPS) );
    
    if rand<LH1/LH2
        Suc(k,s)=Suc(k,s)+1;
        b(Suc(k,s)+1,k,s)=b1;
        LH2=LH1;
%         if Suc(k,s)>1000
%         bsamp(Suc(k,s)-1000,k,s)=b1;
%         end    
        if mod(Suc(k,s)+Fail(k,s),1000)==0
        if Suc(k,s)./(Suc(k,s)+Fail(k,s))>0.5
            stdev=stdev*2;
        else if Suc(k,s)./(Suc(k,s)+Fail(k,s))<0.3
                stdev=stdev/2;
            end
        end
        end
    else Fail(k,s)=Fail(k,s)+1;
    end

    
end

end
PP{s}=ProbPos;
end

bsamp=reshape([b(1002:end,1,:); b(1002:end,2,:)],20000,4);

%% Plot figures

figure
xp=0:0.1:25;
tlabel=["Floor (PFU/ml)" "Walls (PFU/ml)" "Trough (PFU/ml)" "Faeces (PFU/ml)"];

for j=1:4
    hold off
subplot(1,4,j)
plot(PP{j}(:,1),PP{j}(:,2),'*','color','red','markersize',5)
hold on
 patch([xp, flip(xp,2)],...
              [1-exp(-prctile(bsamp(:,j),95).*xp), flip(1-exp(-prctile(bsamp(:,j),5).*xp),2)],...
              'b','linestyle','none','FaceVertexAlpha',0.5,...
              'FaceAlpha',0.5);
hold on
plot(xp,1-exp(-median(bsamp(:,j)).*xp),'color','blue','linewidth',5);
xlim([0 25])


ax = gca;
ax.FontSize = 16;


xlabel(tlabel(j),'interpreter','latex','FontSize',18);
if j==1
ylabel('Prob. of detection','interpreter','latex','FontSize',18);
end

end


figure(1)
set(1,'paperunits','centimeters');
set(1,'papersize',[29 7]);
set(1,'paperposition',[-2 -0 33 7]);
% fig1='Env_detection'
% print('-dtiff','-r300',[fig1,'.tif'])
% print(1,'-dpdf',[fig1,'.pdf']);
% savefig([fig1,'.fig']);


