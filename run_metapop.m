%Example of a source-sink metapopulation model that runs population dynamics on
%multiple patches connected by stochastic larval exchanges and affected by
%stochastic disturbance impacts (cycloens, bleaching, COTS, floods; all have 
%stochastic temporal regimes and stochastic damage)

%Basic parameters--------------------------------------
patches=3;
timesteps=100;
patch_populations=zeros(patches,timesteps);%to hold the population numbers

%Basic setup--------------------------------------------------------------
populations=zeros(patches, timesteps);
initial_pop=1000;%initial population
populations(:,1)=initial_pop;
K(1:patches,1)=10^4;%carrying capacity
extinctthresh=10;%when to assume a population collapses due to low density

%Set connectivity type using listed biases and then e.g. adding white
%noise-----------------------------------------------------------------
conncttype='vr';
biases=[0.1 0.01 0;
        0.001 0.1 0.1;
        0 0.01 0.1]; %if 0 then default
connct=connectivity(conncttype, patches, timesteps, biases);


%Set population growth and mortality--------------------------------------
natgrowth =zeros(patches,1);
natgrowth(:,1)=[0.2; 0.2; 0.2];%around 0.2 is balanced wiht 0.05 mortality

natadultmort =zeros(patches,1);
natadultmort(:,1)=0.05;
natsettlmort=zeros(patches,timesteps);

%Set up disturbances etc.------------------------------------------------
dodstb=[1 0 0 1];%put in 1 to include cyclone, bleaching, cots, wet years/floods 
climateC=[0 0 0 0];%include effects of climate change for respective disturbance
anthroWQ=1;%worsening WQ due to human impact
pcycl=0.2;%yearly probability of cyclones
cyccatW=[0.3 0.3 0.2 0.2 0.1];%probability of cyclone categories
cycextnW=[0.5 0.35 0.1 0.05];%probability of cyclone footprint: just O1, both O1 and O2, O1 and I, all three
pwet=[1/3 1/3 1/3];%probability of a year being wet, normal, or dry
pblch=0.05;%probability of bleaching, either given here as fixed prob, or usign climate change - careful so it rounds well, timesteps should be in units of 10
blchsevrW=[0.4 0.3 0.2 0.1];%how many DHW, from Hughes for 2,4,6,8 DHW
blchextnW=[0.5 0.3 0.2];%spatial bleaching coverage: only O1 hit, O1 and I both hit, all three hit - maybe do individual probabilties after all
pcots=0.06;%yearly probability of COTS initiation on O1
cotsyrlmort=0.6;%yearly mortality due to COTS outbreak
indvprobs=1;pcycl=zeros(patches,timesteps);pcycl(1:2,3)=1;cyccatW=zeros(patches,timesteps);cyccatW(1,3)=5;cyccatW(2,3)=5;pcycl(3,5)=1;cyccatW(3,5)=3;
fixedcycl=0;
fixedblch=0;
fixedcots=0;
fixedwety=0;

%Run distrubance regimes--------------------------------------------------
[dstbcmlmort, natsettlmort, dstbdistb] = disturbances(patches, timesteps, dodstb, climateC, natsettlmort, anthroWQ, pcycl, cyccatW, cycextnW, pwet, pblch, blchsevrW, blchextnW, pcots, cotsyrlmort, indvprobs, fixedcycl, fixedblch, fixedcots, fixedwety);


%Run population dynamics---------------------------------------------------
[populations] = pop_dyn( patches, timesteps, populations, connct, dstbcmlmort, natgrowth, natsettlmort, natadultmort, K, compt, extinctthresh );


%Plot stochastic population growth----------------------------------------
lambda=zeros(timesteps,1);
for lm=2:timesteps
    %lambda(lm,1)=(sum(populations(:,(lm)))/(initial_pop*patches))^(1/(lm));
    lambda(lm,1)=(populations(3,(lm))/(initial_pop))^(1/(lm));
end

%Plot individual trajectories for three patches----------------------------
figure;
plot(1:timesteps, populations(3,:),'-o')
hold;
plot(1:timesteps, populations(2,:),'-o')
plot(1:timesteps, populations(1,:),'-o')
figure;
plot(2:timesteps,lambda(2:end),'-o')
