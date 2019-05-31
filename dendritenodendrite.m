N=100; %ms

dt=0.01; %ms

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
tauc=0.7039;  % ms
taus=32.91;   % ms
c1=1;         % microF/cm^2
p=0.0759;     % somatic S / total S
gc=p*(p-1)*(c1/taus-c1/tauc)*1000; % 1000, cm^2
G=gl +gc/p/(1-p);       % total conductance: coupling conductance + leak conductance

vspike=0;
Ihold=-500;

% compensation for the dendrite
%Ihold=Ihold/( 1-(1-p)*c1*1000/G*(1/tauc-1/taus) );

sigma=0;
corr=cormax;
n_spike=0;
temp=0;


% start from spiking state
v(1)=-50;
w(1)=0;

vd(1)=-50;
V(1)=0;           % V is vs=vd, V=vs-vd
wd(1)=0;

vn(1)=-50;
wn(1)=0;

t(1)=0;
input(1)=Ihold;

INT(1)=0;        % dendrite filter
 
for i=2:1:round(N/dt)
     t(i)=(i-1)*dt;    
     
temp=temp-dt/corr*temp + sqrt(2*dt/corr)*random('normal',0,1,1,1);
input(i)=Ihold;
%+ temp*sigma;

 % dendrite diff equation
 vd(i)=dt/c*( -gl*(vd(i-1)-el) +gl*delta*exp((vd(i-1)-vt)/delta) -wd(i-1) + input(i) -gc/p*V(i-1) ) + vd(i-1);
 V(i)=dt/c*( -G*V(i-1) +input(i) ) + V(i-1);
 wd(i)=dt/tauw*(a*(vd(i-1)-el)-wd(i-1)) + wd(i-1);

 % dendrite integral
 %INT(i)=-gc/p*exp(-t(i)*G/c)*trapz(t,input.*exp(t*G/c))/c; % dendrite integral
 %v(i)=dt/c*(-gl*(v(i-1)-el)+gl*delta*exp((v(i-1)-vt)/delta)-w(i-1)+input(i)+INT(i)) + v(i-1);
 %w(i)=dt/tauw*(a*(v(i-1)-el)-w(i-1)) + w(i-1);
 
 % no dendrite *( 1-(1-p)*c1*1000/G*(1/tauc-1/taus) )
 vn(i)=dt/c*(-gl*(vn(i-1)-el)+gl*delta*exp((vn(i-1)-vt)/delta)-wn(i-1)+input(i)*( 1-(1-p)*c1*1000/G*(1/tauc-1/taus) )) + vn(i-1);
 wn(i)=dt/tauw*(a*(vn(i-1)-el)-wn(i-1)) + wn(i-1);
 
          %   if  v(i)>=vspike
          %       v(i-1)=0;            % add sticks to the previous step     
          %       v(i)=vreset;
          %       w(i)=w(i) + b;
          %   end

             if  vd(i)>=vspike
                 vd(i-1)=0;            % add sticks to the previous step     
                 vd(i)=vreset;
                 wd(i)=wd(i) + b;
             end

              if vn(i)>=vspike
                 vn(i-1)=0;            % add sticks to the previous step     
                 vn(i)=vreset;
                 wn(i)=wn(i) + b;
             end
             
end

plot(t,vd,t,vn);
legend('v','vd','vn');
xlabel('time, ms');
ylabel('V, mV');