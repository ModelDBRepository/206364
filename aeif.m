N=100; %ms
trial=5;

dt=0.1; %ms
dcorr=1;
dsigma=1;

n_spike=0;
sum=0;

sigma=-dsigma; %initial variation

sigmax=10;
cormax=2;

c=217; i=0; gl=12.8;
el=-55.144; vt=-56.252;
delta=0.77; vreset=-68;
vs=0; gs=0.01;

a=35.4; tauw=7.5; b=323;

% dendritic filtering
gc=0.1;      % mS
tauc=0.15;   % time constant for coupling


v(1)=-68;
w(1)=0;
t(1)=0;
input(1)=0;

vspike=0;
Ihold=-50;

sigma=50;

temp=0;
corr=cormax;

for i=2:1:round(N/dt)
             t(i)=(i-1)*dt;
             
     temp=temp-dt/corr*temp + sqrt(2*dt/corr)*random('normal',0,1,1,1); 
     input(i)=Ihold + temp*sigma;
     
     % dendrite trapez calculation
   % v(i)=dt/c*(-gl*(v(i-1)-el)+gl*delta*exp((v(i-1)-vt)/delta)-w(i-1)+input(i)-gc/c*exp(-t(i)/tauc)*trapz(t,exp(t/tauc).*input)) + v(i-1);
   
    % no dendrite rectangular calculation
   int=exp(t/tauc).*input;
    v(i)=dt/c*(-gl*(v(i-1)-el)+gl*delta*exp((v(i-1)-vt)/delta)-w(i-1)+input(i)-gc/c*exp(-t(i)/tauc)*sum(t.*int) ) + v(i-1);
    
    % dendrite
    % v(i)=dt/c*(-gl*(v(i-1)-el)+gl*delta*exp((v(i-1)-vt)/delta)-w(i-1)+input(i)) + v(i-1);
    
    
     w(i)=dt/tauw*(a*(v(i-1)-el)-w(i-1)) + w(i-1);
    
             if  v(i)>vspike
                 n_spike=n_spike+1;
                 v(i)=vreset;
                 w(i)=w(i) + b;
             end
                    
end
 
freq=n_spike/N*1000; %averaged frequency after N iteration
