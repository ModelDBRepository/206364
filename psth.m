N=60001; % total simulation time, ms 

sw=1200;  % duration of one sweep
dsw=1200; % change of the sweep

av=round(N/sw); % number of sweeps

dt=0.1; %ms


% model parameters
c=217; i=0; gl=12.8;
el=-55.144; vt=-56.252;
delta=0.77; vreset=-68;
a=35.4; tauw=7.5; b=323;


% parameters for external biexponential intput
Am=200; % pA, optimal for inhibition
taus1=1.5; % ms, rise constant
taus2=10; % ms, decay constant
ts=600;    % STIMULUS TIME!

% parameters for external gaussian input
% AA=50; % pA
% C=500; % sigma in Gaussian


% number of spikes
tT=20; % delta bin, ms
bin=20; % initial size of a bin, ms
A=zeros(1,round(sw/tT));

vspike=0;
Ihold=-90;


sigma=20;
corr=2;
temp=0;
j=1;
k=0;


% initial conditions
v(1)=-57;
w(1)=0;
input(1)=Ihold;
in(1)=0;
m(1)=0;
t(1)=0;

% zero initial conditions for external stimuli

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
    % in(i)= AA*exp(-(t(i)-ts)^2/2/C);
     %input(i)=Ihold + temp*sigma + in(i);

    % no dendrite
     v(i)=dt/c*(-gl*(v(i-1)-el)+gl*delta*exp((v(i-1)-vt)/delta)-w(i-1)+input(i)) + v(i-1);
     w(i)=dt/tauw*(a*(v(i-1)-el)-w(i-1)) + w(i-1);

     % binning             
if t(i)>=bin;
    bin=bin + tT;
    j=j+1;    
end;

% sweep
if  t(i)>=sw;
    ts=ts+dsw;
    sw=sw + dsw;
    j=1; % reset the bin number
    v(i)=-57;
    w(i)=0;     % reset init. cond. after each sweep
end;             


     
             if  v(i)>vspike
                 v(i-1)=0;            % add sticks to the previous step
                 v(i)=vreset;
                 w(i)=w(i) + b;
                 
                 % +1 spike to the bin
                 A(j)=A(j)+1;%    % average number of spikes in each bin
             end
             
                       
             
end

A=A/av/tT;

%save psth100000.mat