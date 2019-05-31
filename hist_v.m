


% no spikes in the trace (visual interspection)

%% PC cell Data (CHECK!)
dt=0.02; % sampling rate, 50 kHz
v_rest=[v(6000/dt:1:7000/dt);v(8000/dt:1:9000/dt);v(10000/dt:1:11000/dt);v(12000/dt:1:13000/dt)];
v_spike=[v(1:1:5500/dt);v(7000/dt:1:8000/dt);v(9000/dt:1:10000/dt);v(11000/dt:1:12000/dt);v(13000/dt:1:21000/dt);];
%%

%% aEIF model
dt=0.05; % time step, "sampling rate" 20 kHz
v_rest=[v(5000/dt:1:6000/dt),v(7000/dt:1:8000/dt),v(9000/dt:1:10000/dt),v(20000/dt:1:21000/dt)];
v_spike=[v(1:1:5000/dt),v(11000/dt:1:12000/dt),v(13000/dt:1:14000/dt),v(15000/dt:1:20000/dt),];
%%

%% Eriks model
dt=0.05; % time step, "sampling rate" 20 kHz
v_rest=[v(7000/dt:1:8000/dt);v(9000/dt:1:10000/dt)];
v_spike=[v(1:1:6000/dt);v(10000/dt:1:11000/dt);v(11000/dt:1:12000/dt);v(12000/dt:1:14000/dt);v(15000/dt:1:16000/dt);v(17000/dt:1:20000/dt)];
%%

%%
%histNorm({v_spike,v_rest},'Legend',{'Spiking','Silent'});
figure
histNorm({v_spike,v_rest},'Legend',{'Spiking','Silent'},'FracBinSize', 20);
set(gca,'Fontsize',30);
xlabel('Membrane potential [mV]');
ylabel('Fraction of time');
title('This is ...');
axis([-75 -30 0 0.25]);
ylim auto;
box off;
%%