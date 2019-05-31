N=100000; % total simulation time, ms 

sw=100; % duration of one sweep
dsw=100; % change of the sweep

av=round(N/sw); % number of sweeps

dt=0.1; %ms


% model parameters
c=217; i=0; gl=12.8;
el=-55.144; vt=-56.252;
delta=0.77; vreset=-68;

a=35.4; tauw=7.5; b=323;


% external intput
Am=10; % pA, optimal for inhibition
taus1=1.5; % ms, rise constant
taus2=10; % ms, decay constant
ts=50;


% number of spikes
A(1)=0;
tT=1; % delta bin, ms
bin=tT; % initial size of a bin, ms

vspike=0;
Ihold=-90;


sigma=50;
corr=2;
temp=0;
j=0;
k=0;


% random initial conditions
% variance = 5% of the mean
v(1)=-56.92;
w(1)=-63.08;
input(1)=Ihold;

t(1)=0;

% zero initial conditions for external stimuli
in(1)=0;
m(1)=0;


for i=2:1:round(N/dt)
     t(i)=(i-1)*dt;              
      
     % additional stimulus
     
     % delta function approximation
     if t(i)==ts        
         stim(i)=1/dt;
     else
         stim(i)=0;
     end;
    
     
     % generate external stimuli
     m(i)=dt/taus1/taus2*(Am*(1-in(i-1))*stim(i)/K(1/taus1,1/taus2)-in(i-1)-(taus1+taus2)*m(i-1)) + m(i-1);
     in(i)=m(i)*dt + in(i-1);
     
     % generate total stimuli
     temp=temp-dt/corr*temp + sqrt(2*dt/corr)*random('normal',0,1,1,1); 
     input(i)=Ihold + temp*sigma + in(i);
     

    % no dendrite
     v(i)=dt/c*(-gl*(v(i-1)-el)+gl*delta*exp((v(i-1)-vt)/delta)-w(i-1)+input(i)) + v(i-1);
     w(i)=dt/tauw*(a*(v(i-1)-el)-w(i-1)) + w(i-1);

     % binning             
if t(i)>=bin;
    bin=bin + tT;
    j=j+1;
    
    if k==0  % check for the first sweep
        A(j)=0;
    end;
    
end;

% sweep
if  t(i)>=sw;
    ts=ts+dsw;
    sw=sw + dsw;
    j=1; % reset the bin number
    k=1; % marker for the first sweep
    % start from the down-state
    v(i)=-56.92;
    w(i)=-63.08;
end;             


     
             if  v(i)>vspike
                 v(i-1)=0;            % add sticks to the previous step
                 v(i)=vreset;
                 w(i)=w(i) + b;
                 % +1 spike to the bin
                 A(j)=A(j)+1/av/tT;   % average number of spikes
             end
             
                       
             
end


save psthS50.mat
