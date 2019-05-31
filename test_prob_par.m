N=300000; %ms % to have the same number of averages as in Fig. 6
trial=20;

dt=0.1; %ms

n_spike=0;
sum=0;

dsigma=5;

sigma=-5;  %initial variation

sigmax=200;
corr=2;

% cell 1
c=268; 
gl=8.47;
el=-51.31; 
vt=-53.23;
delta=0.85; 
vreset=-60.35;1

a=37.79; tauw=20.76; b=441.12;

% dendritic filtering
tauc=0.7039;  % ms
taus=32.91;   % ms
c1=1;         % microF/cm^2
p=0.0759;     % somatic S / total S

gc=p*(p-1)*(c1/taus-c1/tauc)*1000; % 1000, cm^2
G=gl +gc/p/(1-p);       % total conductance: coupling conductance + leak conductance

% bassin of attraction
vb=importdata('vb-200.mat');
wb=importdata('wb-200.mat');

% -200, -150, -100

t(1)=0;
Ihold=-150;
input(1)=Ihold;

v(1)=vreset;
w(1)=0;

temp=0;
tt=0;
vspike=0;
time=0; %initially

for z=1:1:round(sigmax/dsigma)
    sigma=sigma+dsigma;
    
          parfor j=1:1:trial
                     
n_spike=0;
temp=0;

% reset after trial
v=-45;           % start from the spiking state
vd=v;
V=0;
w=0;
t=0;
input=Ihold;     
              
 for i=2:1:round(N/dt)
    
    % Ornstein-Uhlenbeck
     temp=temp-dt/corr*temp + sqrt(2*dt/corr)*random('normal',0,1,1,1); 
     input=Ihold + temp*sigma;
 
 % dendrite
% v_old=v; % v(i-1)
% w_old=w; % w(i-1)
% V_old=V;
% v=dt/c*( -gl*(v_old-el) +gl*delta*exp((v_old-vt)/delta) -w_old +input -gc/p*V_old  ) + v_old;
% V=dt/c*( -G*V_old +input ) + Vold;
% w=dt/tauw*(a*(v_old-el)-w_old) + w_old;
     
    %no dendrite
     v_old=v; % v(i-1)
     w_old=w; % w(i-1)
     v=dt/c*(-gl*(v_old-el)+gl*delta*exp((v_old-vt)/delta)-w_old+input) + v_old;
     w=dt/tauw*(a*(v_old-el)-w_old) + w_old;
    
             if  v>vspike
                 v=vreset;
                 w=w + b;
             end
             
     if inpolygon(v,w,vb,wb) == 1
        tt=tt+1;           % total time (in bins) inside the attraction bassin
     end
             
     
   
 end
 
          end
        
    prob(z)=tt*dt/trial/N     % probability of the down state
    tt=0;                     % reset of the total time

end

parsave(sprintf('prob%d.mat', Ihold), A, Am); % save the variables in file
