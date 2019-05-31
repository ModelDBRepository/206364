T=100;  %ms
dt=0.005; %ms

% stimuli
sigma=20; %mV
corr=2;    %ms

C=2;        % microF/cm^2
taus=15; % ms
tauc=0.15;   % ms
p=0.021;     % dimensionless
El=-60;     % mV
Ihold=30;    % pA

gl=C/taus;  % mS/cm^2
gc=p*(p-1)*(gl-C/tauc);

vspike=0;
vreset=-68;
delta=0.77;
vt=-55;

v1n(1)=-60;
v3(1)=0;

temp=0;
time(1)=0;
input(1)=0;

% delta function input

ts=50;

for i=2:1:round(T/dt)
    
    time(i)=(i-1)*dt;
    
    % noisy stimuli
     temp=temp-dt/corr*temp+sqrt(2*dt/corr)*random('normal',0,1,1,1); 
     input(i)=Ihold + temp*sigma*gl;
    
    % delta function input
    
   % if time(i)==ts
   %     input(i)=1/dt;
   % else
   %     input(i)=0;
   % end
    
    % input(i)=1;

% passive membrane    
% with integral
v1n(i)=dt/C*(-gl*(v1n(i-1)-El)-gc/p/C*exp(-time(i)/tauc)*trapz(time,exp(time/tauc).*input)+input(i))+ v1n(i-1);
f(i)=gc/p/C*exp(-time(i)/tauc)*trapz(time,exp(time/tauc).*input);

% active membrane
%v1n(i)=dt/C*(-gl*(v1n(i-1)-El)-gc/p/C*exp(-time(i)/tauc)*trapz(time,exp(time/tauc).*input)+input(i)+gl*delta*exp((v1n(i-1)-vt)/delta))+ v1n(i-1);
%f(i)=gc/p/C*exp(-time(i)/tauc)*trapz(time,exp(time/tauc).*input);

% threshould
%if  v1n(i)>=vspike
%                 v1n(i-1)=0;            % add sticks to the previous step
%                 v1n(i)=vreset;
                % w(i)=w(i) + b;
%             end

                        
end

plot(time,v1n)
   