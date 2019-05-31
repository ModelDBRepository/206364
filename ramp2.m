N=100000; %ms

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

% rescaling
TAUW=tauw/taum;
A=a/gl;
B=b/gl/delta;
VRESET=(vreset-vt)/delta;
VSPIKE=(vspike-vt)/delta;

Ihold=0;    % pA initial current in non-rescaled model

% starting points for variations
% A=1/TAUW;     % start near Bogdanov-takens bifurcation A*TAUW=1 + 0.05 (next step)
TAUW=1/A;     % start near Bogdanov-takens bifurcation A*TAUW=1 + 0.01 (next step)

for j=1:1:350
    
    % vary A
    %A=A + 0.05; % accuracy
    %par(j)=A;
    
    % vary TAUW
    TAUW=TAUW + 0.01; % accuracy
    par(j)=TAUW
    
    a=gl*A;
    tauw=TAUW*taum;
       
    % rescale initial conditions and input
    V(1)=(-68-vt)/delta;     % Vrest
    W(1)=a*(el-vt)/gl/delta; % Wrest 
        
    % check Andronov-Hopf
    if A<=1/TAUW
        break
    end;
    
    dI=0; % input variation

                for i=2:1:round(N/dt)
                    t(i)=(i-1)*dt;  
                    input(i)=Ihold - dI;
                    
     % input rescaling and integration           
             input(i)=input(i)/gl/delta + (1+a/gl)*(el-vt)/delta;          
             V(i)=dt/taum*(-V(i-1)+exp(V(i-1))-W(i-1)+input(i)) + V(i-1);
             W(i)=dt/TAUW/taum*(A*V(i-1)-W(i-1)) + W(i-1);
             
             if  V(i)>=VSPIKE
                 V(i)=VRESET;
                 W(i)=W(i) + B;
                 dI=dI + 0.2;                 % reduce the input, pA
             end
             
             
             % check the rest state
             if (V(i)-V(i-1))/dt == 0 
                 if (W(i)-W(i-1))/dt == 0     
   Idown(j)=(input(i) - (1+a/gl)*(el-vt)/delta)*gl*delta;
   Iup(j)=(gl+a)*(vt-el-delta+delta*log(1+taum/tauw))+delta*gl*(a/gl-taum/tauw);
                     break
                 end
             end
                    
                end

 
 ISR(j)=Iup(j)-Idown(j) 

end;

plot(par,ISR);