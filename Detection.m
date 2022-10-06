function [Dtime, AI, TIp, Detprob]=detection(E,C,Ni,TotInf,Insp,Ninsp,Nsamp,bsamp,N,dt)
%Function to calculate the time from first infection to detection (Dtime),
% the number of animals infected at that time (AI), the proportion of
% infectiousness occurring before that time (TIp), and the probability of
% detection at any time during the outbreak (Detprob)


%% Initialise arrays
Ndays=size(Ni,1)*dt; Nk=size(Ni,2); Etype=size(E,2);
Dtime=NaN(Etype+2,Nk); AI=zeros(Etype+2,Nk);
TIp=NaN(6,Nk,2); Detprob=zeros(5,Ndays/dt,Nk); 

Esum=permute(sum(E,2),[1,3,2]); %Calculate total environmental contamination

FI=randi(Insp*2,Nk,1)./2; %Day of first inspection

%% Calculate for each simulation
for k=1:Nk
    DetectedC=zeros(Ndays,1); DetectedE=zeros(Ndays,Etype,Nsamp); DetectedEall=zeros(Ndays,Nsamp);

    for j=1:Ndays/dt
        if mod(j-FI(k)/dt,Insp/dt)==0 %Day of inspection
            DetectedC(j)=sum(randsample(N,Ninsp)<=C(j,k));
            
            for s=1:Nsamp
                DetectedE(j,:,s)=rand(1,Etype)<1-exp(-bsamp(k,:).*E(j,:,k));
                EtypeR=randi(4);
                DetectedEall(j,s)=rand<1-exp(-bsamp(k,EtypeR).*E(j,EtypeR,k));
            end
            
        end
        if min(sum([DetectedC,sum(DetectedE,3)]))>0
            break
        end
        
    end
    
    for jj=1:Ndays/dt
        DetprobE(jj,:)=1-exp(-bsamp(k,:).*E(jj,:,k));
        DetprobC(jj)=C(jj,k)/N;
    end

    if max(DetectedC>0)
        Dtime(1,k)=find(DetectedC>=1,1)*dt;
        AI(1,k)=Ni(Dtime(1,k)/dt,k); %number of animals infected at detection time
        TIp(1,k,1)=trapz(TotInf(1:Dtime(1,k)/dt,k))/trapz(TotInf(1:end,k));
        TIp(1,k,2)=trapz(Esum(1:Dtime(1,k)/dt,k))/trapz(Esum(1:end,k)); 
    else
        Dtime(1,k)=nan;
        AI(1,k)=nan;
    end
%     
        DetectedEall=sum(DetectedEall,2);
    if max(DetectedEall>0)
        Dtime(6,k)=find(DetectedEall>=1,1)*dt;
        AI(6,k)=Ni(Dtime(6,k)/dt,k); %number of animals infected at detection time
        TIp(6,k,1)=trapz(TotInf(1:Dtime(6,k)/dt,k))/trapz(TotInf(1:end,k));
        TIp(6,k,2)=trapz(Esum(1:Dtime(6,k)/dt,k))/trapz(Esum(1:end,k)); 
    else
        Dtime(6,k)=nan;
        AI(6,k)=nan;
    end
    
    for e=2:Etype+1
        if max(max(DetectedE(:,e-1,:)))==1
            Dtime(e,k)=find(sum(DetectedE(:,e-1,:),3)>0,1)*dt;
            AI(e,k)=Ni(Dtime(e,k)/dt,k); %number of animals infected at detection time
            TIp(e,k,1)=trapz(TotInf(1:Dtime(e,k)/dt,k))/trapz(TotInf(1:end,k)); 
            TIp(e,k,2)=trapz(Esum(1:Dtime(e,k)/dt,k))/trapz(Esum(1:end,k)); 
        else
            Dtime(e,k)=nan;
            AI(e,k)=nan;
            TIp(e,k,1)=nan;
            TIp(e,k,2)=nan;
        end
    end
    
    Detprob(:,:,k)=[DetprobC; DetprobE'];
end
