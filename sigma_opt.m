N=1000; %ms
trial=1;

dt=0.1;

% cell 1
c=268; 
gl=8.47;
el=-51.31; 
vt=-53.23;
delta=0.85; 
vreset=-60.35;

a=37.79; tauw=20.76; b=441.12;

vspike=0;

taum=c/gl;

dsigma=5;
sigma=-5;  % initial variation

sum=0;
n_spike=0;

sigmax=20;

% rescaling
TAUW=tauw/taum;
A=a/gl;
B=b/gl/delta;
VRESET=(vreset-vt)/delta;
VSPIKE=(vspike-vt)/delta;

Ihold=-150;    % pA initial current in non-rescaled model
temp=0;
cor=2;         % ms
sigma=5;       % pA


% starting points for variations
 A=1/TAUW;     % start near Bogdanov-takens bifurcation A*TAUW=1 + 0.05 (next step)
% TAUW=1/A;     % start near Bogdanov-takens bifurcation A*TAUW=1 + 0.01 (next step)

for j=1:1:10
    
    % vary A
    A=A + 0.05; % accuracy
    par(j)=A;
    
    % vary TAUW
   % TAUW=TAUW + 0.01; % accuracy
   % par(j)=TAUW
    
    a=gl*A;
    tauw=TAUW*taum;
       
    % rescale initial conditions and input
    V(1)=(-55-vt)/delta;     % Vrest
    W(1)=a*(el-vt)/gl/delta; % Wrest 
        
    % check Andronov-Hopf
    if A<=1/TAUW
        break
    end;
    
    for z=1:1:round(sigmax/dsigma)
    sigma=sigma+dsigma;
    
  for s=1:1:trial
            n_spike=0;
            temp=0;
    
  for i=2:1:round(N/dt)
 t(i)=(i-1)*dt;  

temp=temp-dt/cor*temp + sqrt(2*dt/cor)*random('normal',0,1,1,1); 
input(i)=Ihold + temp*sigma;                    
                    
     % input rescaling           
             input(i)=input(i)/gl/delta + (1+a/gl)*(el-vt)/delta;
     
     % integration        
             V(i)=dt/taum*(-V(i-1)+exp(V(i-1))-W(i-1)+input(i)) + V(i-1);
             W(i)=dt/TAUW/taum*(A*V(i-1)-W(i-1)) + W(i-1);
             
             if  V(i)>=VSPIKE
                 n_spike=n_spike+1; 
                 V(i)=VRESET;
                 W(i)=W(i) + B;
             end             
             
 end
        
        sum=sum+n_spike;         %total number of spikes in N iterations
  end
         
    freq(z,j)=sum/trial/N*1000     % averaged frequency after N iteration
    sum=0;

    end

end;