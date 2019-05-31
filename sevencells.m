% this script finds the values of delta I for 7 cells (parameters below)
N=100000; %ms

dt=0.1;

cell=9;     % parameter set number

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
delta=0.85; 
vt=-53.23;
gl=8.47;
a=37.79; 
b=441.12;
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

vspike=0;
% rescaling
TAUW=tauw/taum;
A=a/gl;
B=b/gl/delta;
VRESET=(vreset-vt)/delta;
VSPIKE=(vspike-vt)/delta;

Ihold=50;    % pA initial current in non-rescaled model

            % rescale initial conditions and input
            
    V(1)=(-60-vt)/delta;     % Vrest
    W(1)=a*(el-vt)/gl/delta; % Wrest 
        
    % check Andronov-Hopf
    if A<=1/TAUW
        break
    end;

    dI=0; % initial value for reducing factor
    
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
   Idown=(input(i) - (1+a/gl)*(el-vt)/delta)*gl*delta;                             % pA
   Iup=(gl+a)*(vt-el-delta+delta*log(1+taum/tauw))+delta*gl*(a/gl-taum/tauw);      % pA
                     break
                 end
             end
                    
end

 
 ISR=Iup-Idown                                                                   % in pA
 