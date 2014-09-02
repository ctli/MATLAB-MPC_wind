%% compare fft of c4=0 & c4 = 55

clear all
close all
clc
format compact


%% ================================================================== %%
load 'Q75_N5_eta95_C4_0';

figure(10); clf;
plot(u_sim, '-');
xlim([0 168]);
ylim([0 800]);
xlabel('Time (hr)');
ylabel('Wind Scheduling (MW)');
title('c4=0');
defaultratio;

%% ==============================
Y = fft(u_sim);
Y(1)=[]; % Y(1), is simply the sum of the data, and can be removed

n=length(Y);
power = abs(Y(1:floor(n/2))).^2;

nyquist = 1/2;
freq = (1:n/2)/(n/2)*nyquist; % [1/hr]

period=1./freq; %[hr]

figure(11); clf;
plot(freq,power);
xlabel('Frequency (1/hr)');
ylabel('Power of wind scheduling');
title('Periodogram');
ylim([-0.025 2.5]*1e11);
xlim([-0.005 0.505]);
set(gca, 'tickdir', 'out');
defaultratio;

figure(12); clf;
plot(period,power);
xlabel('Period (Hr)');
ylabel('Power of wind scheduling');
ylim([-0.025 2.5]*1e11);
xlim([-100 9000]);
set(gca, 'tickdir', 'out');
defaultratio;

%% ==============================
du = diff(u_sim);
U = fft(du);
U(1)=[];

n=length(U);
power = abs(U(1:floor(n/2))).^2;

nyquist = 1/2;
freq = (1:n/2)/(n/2)*nyquist; % [1/hr]

period=1./freq; %[hr]

figure(13); clf;
plot(freq,power);
xlabel('Frequency (1/hr)');
ylabel('Power of du');
title('Periodogram');
ylim([-0.1 10]*1e8);
xlim([-0.005 0.505]);
set(gca, 'tickdir', 'out');
defaultratio;

figure(14); clf;
plot(period,power);
ylabel('Power of du');
xlabel('Period (Hr)');
ylim([-0.1 10]*1e8);
xlim([-100 9000]);
set(gca, 'tickdir', 'out');
defaultratio;


%% ================================================================== %%
load 'Q75_N5_eta95_C4_55';

figure(21); hold on;
plot(u_sim, 'r-');
xlim([0 168]);
ylim([0 800]);
xlabel('Time (hr)');
ylabel('Wind Scheduling (MW)');
title('c4=30/55=0.5455');
defaultratio;

%% ==============================
Y = fft(u_sim);
Y(1)=[];

n=length(Y);
power = abs(Y(1:floor(n/2))).^2;

nyquist = 1/2;
freq = (1:n/2)/(n/2)*nyquist; % [1/hr]

period=1./freq; % [hr]

figure(22); hold on;
plot(freq,power, 'r');
xlabel('Frequency (1/hr)');
ylabel('Power of wind scheduling');
title('Periodogram');
ylim([-0.025 2.5]*1e11);
xlim([-0.005 0.505]);
set(gca, 'tickdir', 'out');
defaultratio;

figure(23); clf;
plot(period,power, 'r');
xlabel('Period (Hr)');
ylabel('Power of wind scheduling');
ylim([-0.025 2.5]*1e11);
xlim([-100 9000]);
set(gca, 'tickdir', 'out');
defaultratio;

%% ==============================
du = diff(u_sim);
U = fft(du);
U(1)=[];

n=length(U);
power = abs(U(1:floor(n/2))).^2;

nyquist = 1/2;
freq = (1:n/2)/(n/2)*nyquist; % [1/hr]

period=1./freq; %[hr]

figure(24); clf;
plot(freq,power, 'r');
xlabel('Frequency (1/hr)');
ylabel('Power of du');
title('Periodogram');
ylim([-0.1 10]*1e8);
xlim([-0.005 0.505]);
set(gca, 'tickdir', 'out');
defaultratio;

figure(25); clf;
plot(period,power, 'r');
ylabel('Power of du');
xlabel('Period (Hr)');
ylim([-0.1 10]*1e8);
xlim([-100 9000]);
set(gca, 'tickdir', 'out');
defaultratio;

