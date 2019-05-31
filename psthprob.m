parfor s=0:1:15  % set number, Am=s*10 %pA

N=6000001;   % total simulation time, ms (+1 to take into account the last bin)

sw=6000;  % duration of one sweep
dsw=6000; % change of the sweep

av=round(N/sw); % number of sweeps

dt=0.1; %ms

% cell 1
c=268; 
gl=8.47;
el=-51.31; 
vt=-53.23;
delta=0.85; 
vreset=-60.35;

a=37.79; tauw=20.76; b=441.12;

% parameters for external biexponential intput
Am=s*10; % pA, optimal for inhibition
taus1=1.5; % ms, rise constant
taus2=10; % ms, decay constant
ts=3000;    % STIMULUS TIME!


% number of spikes
A(1)=0;
tT=20; % delta bin, ms
bin=20; % initial size of a bin, ms


vspike=0;
Ihold=-150;


sigma=60;
corr=2;
temp=0;
j=0;
k=0;


v(1)=-55;
w(1)=0;
input(1)=Ihold;

% t(1)=0;

% zero initial conditions for external stimuli
in(1)=0;
m(1)=0;
tt(1)=1;
time=0;

% bassin of attraction
vb=importdata('vb-150.mat');
wb=importdata('wb-150.mat');


for i=2:1:round(N/dt)
     t=(i-1)*dt;              
     % additional stimulus
     
     % delta function approximation
     if t==ts        
         stim=1/dt;
     else
         stim=0;
     end;
    
     % generate external stimuli
     m_old=m;
     in_old=in;
     m=dt/taus1/taus2*(Am*(1-in_old)*stim/K(1/taus1,1/taus2)-in_old-(taus1+taus2)*m_old) + m_old;
     in=m*dt + in_old;
     
     % generate the noise and the whole stimulus
     temp=temp-dt/corr*temp + sqrt(2*dt/corr)*random('normal',0,1,1,1); 
     input=Ihold + temp*sigma + in;

    % no dendrite
     v_old=v;
     w_old=w;
     v=dt/c*(-gl*(v_old-el)+gl*delta*exp((v_old-vt)/delta)-w_old+input) + v_old;
     w=dt/tauw*(a*(v_old-el)-w_old) + w_old;

    if  v>vspike
    v=vreset;
    w=w + b;
    end
     
     % binning             
if t>=bin;
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
if  t>=sw;
    ts=ts+dsw;
    sw=sw + dsw;
    j=0; % reset the bin number
    k=1; % marker for the first sweep
  %  v=-55;      % reset for init. cond. after each sweep
  %  w=0;     
end;             
     
% if belongs to the bassin of attraction
 if inpolygon(v,w,vb,wb) == 1
        tt=tt+1;           % time inside of the attraction bassin for one bin
        time=tt*dt;
 end
    
             
end

% normalization
 A=A/tT/av;
    
 parsave(sprintf('psthprob%d.mat', s*10), A, Am); % save the variables in file
 
end