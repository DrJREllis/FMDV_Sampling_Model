tic
clear variables

%% Changeable parameters
Ntrials=10000; %No. of simulations
Ndays=30; %No. of days
dt=1/24; %Time interval for calculating shedding levels

N=50; %No. of animals

Esize=N; %Size of the environment (=N if linearly scaled with herd size)

Insp=1; %Days between inspections
Ninsp=3; %Number of animals inspected
Nsamp=5; %Number of samples on each part of the environment

beta=1; %Probability of infection through interaction
betaE=1; %Probability of infection through the environment

%% Load parameters
ParVp=zeros(20000,11,5); ParTr=zeros(20000,10,5); ParNi0=zeros(20000,5);
for k=1:5
    [ParVp(:,:,k), ParTr(:,:,k), ParNi0(:,k)]=Load_parameters(k); 
end

E0=[0 0 0 0]; %Initial contamination: [floor, wall, trough, faeces]

%% Inititalise arrays
E=zeros(Ndays/dt+1,4,Ntrials); Cs=zeros(Ndays/dt,Ntrials); Cfirst=zeros(1,Ntrials);
Ni=zeros(Ndays/dt,Ntrials); TotInf=zeros(Ndays/dt+1,Ntrials); dailyI=zeros(Ndays,Ntrials);
InfCause=zeros(N,2,Ntrials); Infectiousness=zeros(N,Ndays/dt+1,Ntrials); 
ProbE=zeros(Ndays/dt,Ntrials); ProbD=zeros(Ndays/dt,Ntrials); Ni0=zeros(1,Ntrials);

%% Run model of transmission

k=1;
while k<=Ntrials    
%     ParVp(:,2)=5+ones(20000,1)*k/1000;
 
    %Choose a dataset (data fitted to different IPs) at random
    dataset=randi(5);
    
    %Create a viral profile for all individuals
    [VP, Ti, samp]=Viral_profile(ParVp(:,:,dataset),N,Ndays,dt); %Create the viral profile of all animals
    
    %Find the initial number of infections
    Ni0(k)=ParNi0(samp,dataset);
    
    %Run the transmission model
    [E(:,:,k),Cs(:,k),Cfirst(k),Ni(:,k),TotInf(:,k),dailyI(:,k),InfCause(:,:,k),...
        Infectiousness(:,:,k),ProbE(:,k),ProbD(:,k)]=Transmission_model(N,Ndays,dt,ParTr(samp,:),Ni0(k),E0,VP,Ti,Esize); %Run the transmission model

%     %Run the model with direct transmission only
%     [Cs(:,k),Cfirst(k),Ni(:,k),TotInf(:,k),dailyI(:,k),InfCause(:,:,k),...
%         Infectiousness(:,:,k)]=Direct_model(N,Ndays,dt,ParTr(samp,1),Ni0(k),VP,Ti); %Run the transmission model

%     %Run the model with environmnental transmission only
%     [E(:,:,k),Cs(:,k),Cfirst(k),Ni(:,k),dailyI(:,k),InfCause(:,:,k)]=...
%       Env_model(N,Ndays,dt,ParTr(samp,2:end),Ni0(k),E0,VP,Ti,Esize); %Run the transmission model


    % Reject outbreaks that do not infect at least half of the herd:
    if Ni(end,k)>N/2
        k=k+1;
    end
end

toc 


    