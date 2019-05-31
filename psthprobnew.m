N=1201; % total simulation time, ms (+1 to take into account the last bin)

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
Am=100; % pA, optimal for inhibition
taus1=1.5; % ms, rise constant
taus2=10; % ms, decay constant
ts=600;    % STIMULUS TIME!


% number of spikes
A(1)=0;
tT=20; % delta bin, ms
bin=20; % initial size of a bin, ms


vspike=0;
Ihold=-90;


sigma=3;
corr=2;
temp=0;
j=0;
k=0;


v(1)=-57;
w(1)=0;
input(1)=Ihold;

% t(1)=0;

% zero initial conditions for external stimuli
in(1)=0;
m(1)=0;
tt(1)=1;
time=0;

% bassin of attraction
vb=importdata('vb-90.mat');
wb=importdata('wb-90.mat');




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
     
     % generate the noise and the whole stimulus
     temp=temp-dt/corr*temp + sqrt(2*dt/corr)*random('normal',0,1,1,1); 
     input(i)=Ihold + temp*sigma + in(i);

    % no dendrite
     v(i)=dt/c*(-gl*(v(i-1)-el)+gl*delta*exp((v(i-1)-vt)/delta)-w(i-1)+input(i)) + v(i-1);
     w(i)=dt/tauw*(a*(v(i-1)-el)-w(i-1)) + w(i-1);

    if  v(i)>vspike
    v(i-1)=0;            % add sticks to the previous step
    v(i)=vreset;
    w(i)=w(i) + b;
    end
     
     % binning             
if t(i)>=bin;
    bin=bin + tT;
    j=j+1;  % come to the next bin

if k == 0
    A(j)=time;
else
    A(j)=A(j) + time;
end;
    tt=0;   % reset time after each bin
    time=0;
end;


% sweep
if  t(i)>=sw;
    ts=ts+dsw;
    sw=sw + dsw;
    j=0; % reset the bin number
    k=1; % marker for the first sweep
   % v(i)=-57;      % reset for init. cond. after each sweep
   % w(i)=0;     
end;             
     
% if belongs to the bassin of attraction
 if inpolygon(v(i),w(i),vb,wb) == 1
        tt=tt+1;           % time inside of the attraction bassin for one bin
        time=tt*dt;
 end
    
             
end

% normalization
 A=A/tT/av;

%save psth100000.mat
% plot(ttt,v,ttt,in/10-85,ttt,(input-in)/10-95)