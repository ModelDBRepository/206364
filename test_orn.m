N=5000; %ms
trial=1;

dt=0.1; %ms

n_spike=0;
sum=0;

dsigma=1;

sigma=-1;  % initial variation

sigmax=100;
cormax=2;

% cell 1
c=268; 
gl=8.47;
el=-51.31; 
vt=-53.23;
delta=0.85; 
vreset=-60.35;
a=37.79; tauw=20.76; b=441.12;


% dendritic filtering
tauc=0.317;       % ms
taus=30.91;       % ms2
c1=66.97;         % microF/cm^2
cm=1;             % microF
S1=c1/cm;         % S_dend, calculated acccording "pulse" fit
p=0.065;          % S_soma / (S_soma + S_dend)
gc=cm*p*(1-p)*(1/tauc-1/taus); % coupling conductance
G=gl +gc/p/(1-p); % total conductance: leak conductance + coupling conductance;


% preallocation
v=zeros(1,round(N/dt));
vd=zeros(1,round(N/dt));
w=zeros(1,round(N/dt));
input=zeros(1,round(N/dt));
freq=zeros(1,round(sigmax/sigma));
sigma=0:dsigma:sigmax;

v(1)=vreset;
vd(1)=vreset;
w(1)=0;
%t(1)=0;
%input(1)=0;% maximum inhibition

vspike=0;
Ihold=-150;
corr=cormax;

for p=1:1:round(sigmax/dsigma)
    
          for j=1:1:trial
            n_spike=0;
            temp=0;

% reset after trial
vd(1)=-45;           % start from the spiking state
V(1)=0;
wd(1)=0;

t(1)=0;
input(1)=0;

sq=sqrt(2*dt/corr); % pre-calculation
            
 for i=2:1:round(N/dt)
     t(i)=(i-1)*dt;              
    
     % Ornstein-Uhlenbeck, generate the noise
     temp=temp-dt/corr*temp + sq*rand(1); 
              
% "pure" input I, with no dendrite
 input(i)=(Ihold + temp*sigma(p))*S1;
 gc=0;

% input with dendrite present
% input(i)=Ihold + temp*sigma(p);
 
% input with compensation for dendrite
% input=S1/(1-gc/p/G)*(Ihold + temp*sigma(p));

    
 vd(i)=dt/c*( -gl*(vd(i-1)-el) +gl*delta*exp((vd(i-1)-vt)/delta) -wd(i-1) + input(i)/S1 -gc/p*V(i-1) ) + vd(i-1);
 V(i)=dt/c*( -G*V(i-1) +input(i)/S1 ) + V(i-1);
 wd(i)=dt/tauw*(a*(vd(i-1)-el)-wd(i-1)) + wd(i-1);
    
 
             if  vd(i)>vspike
                 n_spike=n_spike+1;
                 vd(i)=vreset;
                 wd(i)=w(i) + b;
             end
             
 end
        sum=sum+n_spike;          %total number of spikes in N iterations
  end
        
    freq(p)=sum/trial/N*1000     %averaged frequency after N iteration
    sum=0;

end

% save inh_dendrite.mat
plot(1:dsigma:sigmax,freq);