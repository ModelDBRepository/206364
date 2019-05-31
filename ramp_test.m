N=1000; %ms

dt=0.05;
i=0;

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


v(1)=-68; %Vrest
w(1)=0;

Iup=(gl+a)*(vt-el-delta+delta*log(1+taum/tauw))+delta*gl*(a/gl-taum/tauw);

Ihold=-150; %pA
k=1;        %nA/ms


                for i=2:1:round(N/dt)
                    t(i)=(i-1)*dt;
                 if (t(i)<=round(N/2));
                    input(i)=k*t(i)+Ihold;
                    mem=input(i)-Ihold;
                 else                    
                    input(i)=2*mem - k*t(i) + Ihold;
                 end;
                    
             v(i)=dt/c*(-gl*(v(i-1)-el)+gl*delta*exp((v(i-1)-vt)/delta)-w(i-1)+input(i)) + v(i-1);
             w(i)=dt/tauw*(a*(v(i-1)-el)-w(i-1)) + w(i-1);
             
             if  v(i)>=vspike
                 v(i-1)=0;
                 v(i)=vreset;
                 w(i)=w(i) + b;
                 Idown=input(i);
                 
             end
                    
                end

ISR=Iup-Idown;
 % save ramp.mat
