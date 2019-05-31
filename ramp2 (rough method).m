N=10000; %ms

dt=0.05;

c=217;
gl=12.8;
el=-55.144;
vt=-56.252;
delta=0.77;
vreset=-68;
vspike=0;

a=35.4;
tauw=7.5;
b=323;

taum=c/gl;

% rescaling
TAUW=tauw/taum;
A=a/gl;
B=b/gl/delta;



Ihold=-200; %pA
k=0.4;        %nA/ms

% starting points

%A=2.27;

%TAUW=0.37;

s=0;

for j=1:1:100
    
    % vary A
    A=A+0.01;
    par(j)=A;
    
    % vary TAUW
    %TAUW=TAUW + 0.01;
    %par(j)=TAUW;
    
    a=gl*A;
    tauw=TAUW*taum;
    
    VRESET=(vreset-vt)/delta;
    VSPIKE=(vspike-vt)/delta;
    
   % Iup is not rescaled! Rheobased current for II exc
   % Iup(j)=(gl+a)*(vt-el-delta+delta*log(1+taum/tauw))+delta*gl*(a/gl-taum/tauw);
   
    
    % rescale initial conditions and input
    V(1)=(-68-vt)/delta;     % Vrest
    W(1)=a*(el-vt)/gl/delta; % Wrest
    input(1)=Ihold/gl/delta  + (1+a/gl)*(el-vt)/delta;
        
    % check Andronov-Hopf
    if A<=1/TAUW
        break
    end;

    % ramp generation
                for i=2:1:round(N/dt)
                    t(i)=(i-1)*dt;
                    
                 if (t(i)<=round(N/2));
                    input(i)=k*t(i)+Ihold;
                    mem=input(i)-Ihold;
                 else                    
                    input(i)=2*mem - k*t(i) + Ihold;
                 end;
                 
     % ramp rescaling            
             input(i)=input(i)/gl/delta + (1+a/gl)*(el-vt)/delta;
                    
             V(i)=dt/taum*(-V(i-1)+exp(V(i-1))-W(i-1)+input(i)) + V(i-1);
             W(i)=dt/TAUW/taum*(A*V(i-1)-W(i-1)) + W(i-1);
             
             if  V(i)>=VSPIKE
                 V(i)=VRESET;
                 W(i)=W(i) + B;
                 
                 % first spike
                 
                 if s==0
                     Iup(j)=input(i);
                     s=1;
                 end;
                 
                 % inverse rescaling
                 Idown(j)=(input(i)-(1+a/gl)*(el-vt)/delta)*gl*delta;
                 
             end
                    
                end
 s=0;
 
 ISR(j)=Iup(j)-Idown(j) 

 
%input=zeros;


end;
 % save ramp.mat