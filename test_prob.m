N=10000; %ms
trial=1;

dt=0.1; %ms

n_spike=0;
sum=0;

dsigma=5;

sigma=-dsigma;  %initial variation

sigmax=200;
corr=2;

% cell 1
c=268; 
gl=8.47;
el=-51.31; 
vt=-53.23;
delta=0.85; 
vreset=-60.35;

a=37.79; tauw=20.76; b=441.12;

p=0.00759;

% bassin of attraction
vb=importdata('vb-250.mat');
wb=importdata('wb-250.mat');
down=0;                    % auxilary variable


v(1)=vreset;
w(1)=0;

t(1)=0;
Ihold=-250;
input(1)=Ihold;

temp=0;
tt=0;
vspike=0;
time=0; %initially

for p=1:1:round(sigmax/dsigma)
    sigma=sigma+dsigma;
    
          for j=1:1:trial

              v(1)=vreset;
              w(1)=0;
              
 for i=2:1:round(N/dt)
     t(i)=(i-1)*dt;              
    
     % Ornstein-Uhlenbeck input
     temp=temp-dt/corr*temp + sqrt(2*dt/corr)*random('normal',0,1,1,1); 
     input(i)=Ihold + temp*sigma;
    
   
   %no dendrite
     v(i)=dt/c*(-gl*(v(i-1)-el)+gl*delta*exp((v(i-1)-vt)/delta)-w(i-1)+input(i)) + v(i-1);
     w(i)=dt/tauw*(a*(v(i-1)-el)-w(i-1)) + w(i-1);
    
             if  v(i)>vspike
                 v(i)=vreset;
                 w(i)=w(i) + b;
             end
             
     if inpolygon(v(i),w(i),vb,wb) == 1
        tt=tt+1;           % total time inside of the attraction bassin
        time=tt*dt;
     end
             
   
 end
 
          end
        
    prob(p)=time/trial/N     % probability of the down state
    tt=0;                     % reset of the total time

end

%%
plot(1:dsigma:sigmax,1-prob);
set(gca,'Fontsize',20);
xlabel('\sigma, pA');
ylabel('Spiking probability');

%save downProbH-120.mat
