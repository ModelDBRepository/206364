N=5000; %ms
trial=5;

dt=0.1; %ms

n_spike=0;
sum=0;

dsigma=5;

sigma=-5;  % initial variation

sigmax=200;
cormax=2;

cell=1;     % parameter set number
Ihold=-150; 

% cell 0
if cell==0
taum=38.59;
c=234; 
el=-50.26;
delta=0.92;
vt=-52.55;
gl=6.06;
a=34.08;
b=344.33;
tauw=19.02;
vreset=-60.35;
end;

% cell 1
if cell==1
taum=31.66;
c=268; 
el=-51.31; 
delta=0.84; 
vt=-53.23;
gl=8.47;
a=37.79; 
b=441;
tauw=20.76;
vreset=-60.35;
end;

% cell 2
if cell==2
taum=27.76;
c=171;
el=-49.50;
delta=0.90;
vt=-51.40;
gl=6.16;
a=36.88;
b=646.87;
tauw=5.27;
vreset=-60.35;
end;

% cell 4
if cell==4
taum=56.43;
c=216;
el=-51.95;
delta=0.77;
vt=-56.14;
gl=3.83;
a=44.71;
b=476.08;
tauw=10.48;
vreset=-60.35;
end;

% cell 6
if cell==6
taum=22.47;
c=161;
el=-52.52;
delta=1.24;
vt=-52.42;
gl=7.17;
a=27.28;
b=261.48;
tauw=9.32;
vreset=-60.35;
end;

% cell 8
if cell==8
taum=72.74;
c=208;
el=-55.23;
delta=1.08;
vt=-55.11;
gl=2.86;
a=48.95;
b=314.51;
tauw=20.74;
vreset=-60.35;
end;

% cell 9
if cell==9
taum=28.71;
c=110;
el=-52.47;
delta=0.91;
vt=-57.70;
gl=3.83;
a=42.40;
b=371.78;
tauw=18.35;
vreset=-60.35;
end;


% dendritic filtering
% old parameters
% tauc=0.7039;  % ms
% taus=32.91;   % ms
% c1=1;         % microF/cm^2
% p=0.0759;     % somatic S / total S
% p=0.01;

% cell 1 parameters

% dendritic filtering
tauc=0.317;       % ms
taus=30.91;       % ms2
c1=53.87;         % pF
p=0.065;          % S_soma / (S_soma + S_dend)
cm=1e7;             % pF/cm^2
S1=c1/cm;         % S_dend, calculated acccording "pulse" fit

%
%gc=cm*p*(1-p)*(1/tauc-1/taus); % coupling conductance

% from Sarah

gL=c1/taus;
gc=p*(1-p)*gL*(taus/tauc-1);

% gc=1000;

% To make compensation work - reduce the dt = 0.05 otherwise Euler wouldn't work ???

% equation for a dendrite
G=gL +gc/p/(1-p); % total conductance: leak conductance + coupling conductance;

% with gl calculated by cm
% G=cm/taus +gc/p/(1-p);

freq=zeros(1,round(sigmax/sigma));

vspike=0;
corr=cormax;
t(1)=0;


% compenstion  with gc change
% gc=7.5;      % dendrite conductance fitted manually
% G=gl +gc/p/(1-p);

% Ihold=Ihold/(1-gc/p/G)*S1; % compensate only Ihold

gc=0;

for z=1:1:round(sigmax/dsigma)
    sigma=sigma+dsigma;

          parfor j=1:1:trial
            n_spike=0;
            temp=0;

% reset after trial
v=-40;           % start from the spiking state
V=0;
w=0;

 for i=2:1:round(N/dt)
                
     % Ornstein-Uhlenbeck1, randn - is much faster than random
     temp=temp-dt/corr*temp + sqrt(2*dt/corr)*randn(1,1);    
                   
% "pure" input I, with no dendrite
% input=(Ihold + temp*sigma)*S1; 
% gc=0;

% input to the rescaled model with S1
% input=(Ihold + temp*sigma);

% input with dendrite present
%  input=Ihold + temp*sigma;
  
% input with FULL compensation for dendrite
% input=1/(1-gc/p/G)*(Ihold + temp*sigma);
 
 % mean compensation
  input=(1/(1-gc/p/G))*Ihold + temp*sigma;

  v_old=v; % v(i-1)
  w_old=w; % w(i-1)
  V_old=V; % V(i-1)
     
 v=dt/c*( -gl*(v_old-el) +gl*delta*exp((v_old-vt)/delta) -w_old + input -gc/p*V_old ) + v_old;
 V=dt/c*( -G*V_old +input ) + V_old;
 w=dt/tauw*(a*(v_old-el)-w_old) + w_old;
 
  
 % v(i)=dt/c*(-gl*(v(i-1)-el)+gl*delta*exp((v(i-1)-vt)/delta)-w(i-1)+input(i) ) + v(i-1);
 % w(i)=dt/tauw*(a*(v(i-1)-el)-w(i-1)) + w(i-1);
    
             if  v>vspike
                 n_spike=n_spike+1; 
                 v=vreset;
                 w=w + b; 
             end
             
 end
        sum=sum+n_spike;         %total number of spikes in N iterations
          end
        
    freq(z)=sum/trial/N*1000     % averaged frequency after N iteration
    sum=0;

end

% save inh_dendrite.mat
plot(1:dsigma:sigmax,freq);
set(gca,'FontSize',30);             % set the axis with big font
ylabel('<Hz>');
xlabel('\sigma, pA');

% plot(1:dsigma:sigmax,freq,1:dsigma:sigmax,freq_d,1:dsigma:sigmax,freq_d_comp);
% set(gca,'FontSize',30);             % set the axis with big font
% ylabel('<Hz>');
% xlabel('\sigma, pA');