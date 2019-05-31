T=100;  %ms
dt=0.005; %ms

% stimuli
sigma=20; %mV
corr=2;    %ms

C=2;         % microF/cm^2
taus=15;     % ms
tauc=0.15;   % ms                               /10 TO CHANGE GC
p=0.021;     % dimensionlessprob 50 max.fig
El=-60;      % mV
Ihold=30;     % pA

gl=C/taus;  % mS/cm^2
gc=p*(p-1)*(gl-C/tauc);

gc=0;        % 

vspike=0;
vreset=-68;
delta=0.77;
vt=-55;


v1(1)=-60;
v2(1)=-60;

temp=0;
time(1)=0;
input(1)=0;

for i=2:1:round(T/dt)
    
    time(i)=(i-1)*dt;
    
    % passive membrane
    % equations                
    v1(i)=dt/C*(-gl*(v1(i-1)-El)-gc/p*(v1(i-1)-v2(i-1))+input(i)) + v1(i-1);
    v2(i)=dt/C*(-gl*(v2(i-1)-El)-gc/(1-p)*(v2(i-1)-v1(i-1))) + v2(i-1);
    
   % v1(i)=dt/C*(-gl*(v1(i-1)-El)-gc/p*(v1(i-1)-v2(i-1))+input(i)+gl*delta*exp((v1(i-1)-vt)/delta)) + v1(i-1);
   % v2(i)=dt/C*(-gl*(v2(i-1)-El)-gc/(1-p)*(v2(i-1)-v1(i-1))) + v2(i-1);
    
    % threshold
    
   % if  v1(i)>=vspike
   %              v1(i-1)=0;            % add sticks to the previous step
   %              v1(i)=vreset;
                % w(i)=w(i) + b;
   %          end
    
end
   
plot(time,v1)