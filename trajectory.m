N=21000; %ms

dt=0.05; %ms

n_spike=0; % current number of spikes
sum=0;

cormax=2; % correlation time

% cell 1
c=268; 
gl=8.47;
el=-51.31; 
vt=-52.23;
delta=0.84; 
vreset=-68;

a=37.79; tauw=20.76; b=441;


% dendritic filtering
tauc=0.317;       % ms
taus=30.91;       % ms2
c1=66.97;         % microF/cm^2
cm=1e6;             % microF
S1=c1/cm;         % S_dend, calculated acccording "pulse" fit
p=0.065;          % S_soma / (S_soma + S_dend)

gc=cm*p*(1-p)*(1/tauc-1/taus); % coupling conductance
% gc=0;

G=gl +gc/p/(1-p); % total conductance: leak conductance + coupling conductance;

Ihold=-150;

sigma=0;
corr=cormax;
n_spike=0;
vspike=0;
temp=0;

% initial conditions
v(1)=-60;
w(1)=-300;
input(1)=Ihold;

V(1)=0; % V=vs-vd; vs=vd at t=0

t(1)=0;

m(1)=0;
in(1)=0;


% parameters for external biexponential intput
Am=0; % pA, optimal for inhibition
taus1=1.5; % ms, rise constant
taus2=10; % ms, decay constant
ts=500;    % STIMULUS TIME!


% bassin of attraction
 vb=importdata('vb-150.mat');
 wb=importdata('wb-150.mat');
 %down=0;                    % auxilary variable
 tt=0;
 
 td=1000;                    % time of one noise/silent period

 n=0;
 q=0;
 
for i=2:1:round(N/dt)
     t(i)=(i-1)*dt;  
     
 
     % increasing noise steps
     
     if t(i)>=td
         n=n+1;            
         sigma=n*5;         % noise amplitude on the each step
         
        % FR(n)=n_spike/1;  % record the rate during the noise stimulation, Hz
        % SG(n)=sigma;
         
         if mod(td/1000,2)==0 % zero variance, "pauses"
                     
             q=q+1;            % counter for increasing noises
             
             FR(q)=n_spike/1;  % record the rate during the noise stimulation, Hz
             SG(q)=sigma;
             
             sigma=0;             
             n_spike=0;       % set number of spikes to zero during the pause
         end
         td=td+1000;         % next period for noise
         
     end          
  
    
    % Sarah stimulus
%{
%if t(i)>600
%    Ihold=-1400;
%end
%if t(i)>601
%    Ihold=-400;
%end
%if t(i)>800
%    Ihold=-700;
%end
%if t(i)>1300
%    Ihold=-400;
%end
%}


% GENERATE EXTERNAL STIMULI

% delta function approximation
%{
     if t(i)==ts        
         stim=1/dt;
         % new stimulus time for randomness
         ts=ts + ts;
     else
         stim=0;
     end;

     m(i)=dt/taus1/taus2*(Am*(1-in(i-1))*stim/K(1/taus1,1/taus2)-in(i-1)-(taus1+taus2)*m(i-1)) + m(i-1);
     in(i)=m(i)*dt + in(i-1);
%}


temp=temp-dt/corr*temp + sqrt(2*dt/corr)*randn(1,1);

% no compensation
 input(i)=Ihold + temp*sigma;
 %+ in(i);

 
%{
% compensation
% input(i)=1/(1-gc/p/G)*(Ihold +temp*sigma);
 
 % dendrite aproximation
 % V=v-vd, applying the subthreshould approximation + mean(gc/p*V)
 % v(i)=dt/c*( -gl*(v(i-1)-el) +gl*delta*exp((v(i-1)-vt)/delta) -w(i-1) + input(i) -gc/p*V(i-1) ) + v(i-1);
 % V(i)=dt/c*( -G*V(i-1) +input(i) ) + V(i-1);
 % w(i)=dt/tauw*(a*(v(i-1)-el)-w(i-1)) + w(i-1);

%} 
    % no dendrite
  v(i)=dt/c*(-gl*(v(i-1)-el)+gl*delta*exp((v(i-1)-vt)/delta)-w(i-1)+input(i) ) + v(i-1);
  w(i)=dt/tauw*(a*(v(i-1)-el)-w(i-1)) + w(i-1);

 
             if  v(i)>=vspike
                 n_spike=n_spike+1;
                 v(i-1)=0;            % add sticks to the previous step     
                 v(i)=vreset;
                 w(i)=w(i) + b;
             end
             
    %  if inpolygon(v(i),w(i),vb,wb) == 1   % time inside of the attraction bassin
    %      tt=tt+1;          
    %      time=tt*dt;
    %  end      
             
end


%% Plot
subplot(3,1,1)
plot(t,v);
xlabel('time, ms');
ylabel('V, mV');

subplot(3,1,2)
plot(t,input);
xlabel('time, ms');
ylabel('In, pA');

subplot(3,1,3)
plot(SG,FR);
xlabel('Sigma, pA');
ylabel('Rate, Hz');
%%

%plot(t,in);
% plot(v,w,vb,wb,'*')