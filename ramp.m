N=3000; %ms

dt=0.05;

% cell 1
c=268; 
gl=8.47;
el=-51.31; 
vt=-53.23;
delta=0.85; 
vreset=-60.35;

a=37.79; tauw=20.76; b=441.12;

vspike=0;

nup=0; ndown=0;

j=0; s=0; nsp=0; sp=0;

time(1)=0;

Ihold=-640; %pA
k=1.5;

tflat=300; %ms

v(1)=-65; %Vrest
w(1)=-520;
input(1)=Ihold;


                for i=2:1:round(N/dt)
                    t(i)=(i-1)*dt;
                
                    if t(i)<=tflat
                     input(i)=Ihold;
                    else
                        
                 if (t(i)<=round(N/2));
                    input(i)=k*t(i)+Ihold-k*tflat;
                    mem=input(i)-Ihold;
                 end;
                 
                 if (t(i)>=round(N/2));
                    input(i)=2*mem - k*t(i) +Ihold +k*tflat;
                 end;
                 
                 
                 if t(i)>=(N-tflat)
                     input(i)=Ihold;
                 end;
                 
                    end;
                    
             v(i)=dt/c*(-gl*(v(i-1)-el)+gl*delta*exp((v(i-1)-vt)/delta)-w(i-1)+input(i)) + v(i-1);
             w(i)=dt/tauw*(a*(v(i-1)-el)-w(i-1)) + w(i-1);
             
               if t(i) == N/2  % irradiation of all markers
                     nsp=0;  
                     sp=0;
                     j=0;
                 end;
             
             
             if  v(i)>=vspike
                 v(i-1)=0;
                 v(i)=vreset;
                 w(i)=w(i) + b;
                 nsp=nsp + 1;
      
                             
                 
                 if t(i) < N/2; % UP ramp
                     
                     if nsp==1;        % first spike
                         sp=sp+1;
                         tsp(sp)=t(i);
                         insp(sp)=input(i);
                         INUP(1)=input(i);
                     end;
                     
                     if nsp==2;        % second spike, first ISI
                         sp=sp+1;
                         tsp(sp)=t(i);
                         insp(sp)=input(i);
                         j=j+1;
                         FRUP(j)=1000/(tsp(2)-tsp(1));
                     end;
                     
                     if nsp>2
                         sp=sp+1;
                         tsp(sp)=t(i);
                         insp(sp)=input(i);
                        
                         j=j+1;
                         FRUP(j)=2*1000/(tsp(sp)-tsp(sp-2)); % Hz
                         INUP(j)=insp(sp-1);
                     end
                     
                     
                     
                 end;      
      
                 
                 
                 
               if t(i) > N/2 % DOWN ramp
                   
               if nsp==1;               % first spike
                         sp=sp+1;
                         tsp(sp)=t(i);
                         insp(sp)=input(i);
                     end;
                     
                     if nsp==2;        % second spike, first ISI
                         sp=sp+1;
                         tsp(sp)=t(i);
                         insp(sp)=input(i);
                     end;
                     
                     if nsp>2
                         sp=sp+1;
                         tsp(sp)=t(i);
                         insp(sp)=input(i);
                        
                         j=j+1;
                         FRDOWN(j)=2*1000/(tsp(sp)-tsp(sp-2)); % Hz
                         INDOWN(j)=insp(sp-1);          
                     end    
                   
               end
                 
             
                    
                end; 
                
                end; 
 % save ramp.mat

j=j+1;                         % Addition of the last spike manually
INDOWN(j)=insp(sp);
FRDOWN(j)=1000/(tsp(sp)-tsp(sp-1));

% flip input and rate values
INDOWN=fliplr(INDOWN);
FRDOWN=fliplr(FRDOWN);


% making the step
FRDOWN=[0,FRDOWN];
FRUP=[0,FRUP];

INDOWN=[INDOWN(1)-0.1,INDOWN];
INUP=[INUP(1)-0.1,INUP];

% make linear fit
pUP=polyfit(INUP(2:end),FRUP(2:end),1);
pDOWN=polyfit(INDOWN(2:end),FRDOWN(2:end),1);

% plot the F-I curve with steps
plot(INDOWN,[0,polyval(pDOWN,INDOWN(2:end))],'blue',INUP,[0,polyval(pUP,INUP(2:end))],'red')
%plot(t,v);