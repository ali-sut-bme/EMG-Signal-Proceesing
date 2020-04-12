%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main.m: Main code for EMG Signal Processing 
% author:Ali Abedi
% Master of Biomedical Engineering From Tehran University of Medical Science
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start
clear;clc;
close all;
%%load emg data
emg=load('Raw.mat');
emg=emg.y;
emg=emg';
emg1=emg(:,2); % Get Data Channel
%% Plot Raw Emg Signal
T=load('time.mat');
T=T.m;
T=T(1:72513,:);
plot(T,emg1);
title('Raw EMG Signal');
xlabel('Time(s)');
ylabel('Amplitude(microvolt)')
det = detrend(emg1); % Remove Any DC Offset of The Signal
L0  = length(emg1); % Signal Length
%% Removing Transition values
%time0 = T(round(.005*L0):L0,1);
%emg0 = det(:,round(.005*L0):L0,1);
L = L0-round(.005*L0)+1;
EMG=emg1;
%% Normalizing data:
M = max(EMG(:));  %Find the maximum amount of data
for i=1:L
    EMG(i)=EMG(i)/M;
end
%% Find sampling frequency
fs = 1/(T(2)-T(1));
%%  Finding envelop:
Rec_EMG = abs(EMG);     %Rectification of the EMG signal
A=0.25;
fc=2.19;
wc=fc/(fs/2);
[b,a]=butter(5,wc,'low');
EMG_envelop0=filter(b,a,Rec_EMG+A);  
EMG_envelop= smooth(EMG_envelop0,'sgolay',4)-A;
%% Normalizing envelop:
[M,I] = max(EMG_envelop);  %Find the maximum amount of data
for i=1:L
    EMG_envelop(i)=EMG_envelop(i)/M;
end
%%
figure;
subplot(3,2,1)
plot(T,EMG)
xlabel('sample')
ylabel('Value')
legend('Original')
subplot(3,2,2)
plot(EMG)
axis([0 L0 -1.1 1.1])
xlabel('sample')
ylabel('Value')
legend('Removed Transition values')
subplot(3,2,3)
plot(Rec_EMG)
axis([0 L0 -1.1 1.1])
xlabel('sample')
ylabel('Value')
legend('Rectificated Values')
subplot(3,2,4)
plot(EMG_envelop)
axis([0 L0 -1.1 1.1])
xlabel('sample')
ylabel('Value')
legend('Envelop')
subplot(3,2,[5 6])
plot(EMG)
axis([0 L0 -1.1 1.1])
xlabel('sample')
ylabel('Value')
hold on
plot(EMG_envelop,'r','LineWidth',2)
hold off
legend('Removed Transition values','Envelop')
%% RMS Changes Algorithm for Onset Offset Detection:
EMG=emg1;
L_window = round(L/100);
for i=1:98
     RMS(i)=sqrt(mean((EMG(i*L_window:(i+1)*L_window).^2)));
end
onset = zeros(L,1);
offset = zeros(L,1);
m=1;
n=1;
while m<96
    if 2.5*RMS(m) < RMS(m+2)
        onset((m+2)*L_window)=200;
        m=m+10;
    else
    m=m+1;
    end
end
[~,On] = findpeaks(onset);
while n<97
    if RMS(n) > 1.8*RMS(n+1)
        offset((n+1)*L_window)=200;
        n=n+5;
    else
    n=n+1;
    end
end
[~,Off] = findpeaks(offset);
%% 
figure
plot(EMG)
xlabel('sample')
ylabel('Value')
hold on
plot(On,0,'o','markerfacecolor','g')
hold on
plot(Off,0,'o','markerfacecolor','r')
legend('EMG','Onset','Offset')
%% Plot RMS Signal
t = 1:length(RMS);
figure;
plot(t,RMS)
xlabel('time');
ylabel('Amplitude');
title('EMG Signal RMS')
%% Frequency Spectrum of EMG Signal
Y = fft(EMG);
L = length(EMG);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;
figure;
plot(f,P1)
xlabel('Frequency (Hz)')
ylabel('Amplitude (microvolts)')
title('Frequency Spectrum of Raw EMG')