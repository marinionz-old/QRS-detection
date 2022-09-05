%%
% IDIM-2019-20 Applications of Digital Signal Processing
%

%% Pan-Tompkins algorithm
% QRS Detection
% 1. Bandpass filter
% 2. Derivative filter
% 3. Square filter
% 4. MA filter
% 5. Threshold
% 6. Detection
% 7. Outliers detection

clearvars;
close all;
clc

%** load your signal
load('Pan-Tompkins.mat');

% Plot the signal
figure(1);
subplot(311)
plot(t,ekg)
title('ECG'); xlabel('Time (s)')
xlim([min(t) max(t)])

%** type the sampling frequency
fs = 500;              % Sampling rate
N = length (ekg);       % Signal length

%** compute DFT
% Watch the DC component
y = fft(ekg,floor(N/200));
subplot(312)

% Make an frequency grid
f= linspace(0,fs,floor(N/200)); % Frequecy grid.
plot(f,abs(y)) % We use the absolute value, as we get complex numbers from
% the Fast Fourier Transform
title('FFT including DC component'); xlabel('Frequency (Hz)')

% Remove the DC and represent the DFT again
ekg_AC = ekg-mean(ekg);
y = fft(ekg_AC,floor(N/200));
subplot(313)
plot(f,abs(y))
title('FFT without DC component'); xlabel('Frequency (Hz)')

% As a result, in the first figure, we can observe that the second subplot 
% and the third one differ in the maximum magnitude of the frequency
% amplitudes, which can be observed in the second subplot as the one
% corresponding to the DC component.

% Once the DC component is removed (mean substraction), all frequency
% magnitudes hold a similar order for low frequencies.

% Both 2nd and 3rd subplot hold the principle of symmetry of the FFT and
% they both span from 0 to the sampling frequency (Fs = 500 Hz in this
% case)


%% Signal treatment

figure(2);clf

% Chunck of the signal. Just for representation purposes
seg_plot = 5*fs;

% 1. Bandpass filter. Select the parameters yourself

minlowp=0.5; % Minimum frequency used in the lab (0.5  Hertz)
maxlowp=150; % Maximum frequency used in the lab (150 Hertz)

[b , a] = butter(1,[minlowp maxlowp]/(fs/2));
ekg_band = filtfilt(b,a,ekg);
subplot(221)
plot(t(1:seg_plot),ekg_band(1:seg_plot))
title('Bandpass filter');xlabel('Time(s)');ylabel('Amplitude')
xlim([t(1) t(seg_plot)])

% 2. Derivative Filter. Recall from seminar 03 in DSP the impulse response

b = [1 -1]; % It detects sudden differences of amplitude in the signal
% (like a derivative does)

ekg_der = filtfilt(b,1,ekg_band); % a coefficient = [1]
subplot(222)
plot(t(1:seg_plot),ekg_der(1:seg_plot))
title('Derivative filter');xlabel('Time(s)');ylabel('Amplitude')
xlim([t(1) t(seg_plot)])

% 3.Squaring
ekg_sq = ekg_der.^2;
subplot(223)
plot(t(1:seg_plot),ekg_sq(1:seg_plot))
title('Square');xlabel('Time(s)');ylabel('Amplitude')
xlim([t(1) t(seg_plot)])

% 4. Moving Average. Recall from seminar 03 in DSP the impulse response
t_ma = 100; % 2*t_ma = ms Window as Fs = 500 Hz. I have tried different 
% windows and this one was proven to be a good one in terms of RR intervals
% time
b = ones(t_ma,1)/t_ma;
ekg_ma = filtfilt(b,1,ekg_sq);
subplot(224)
plot(t(1:seg_plot),ekg_ma(1:seg_plot))
title('MA filter');xlabel('Time(s)');ylabel('Amplitude')
xlim([t(1) t(seg_plot)])

% 5. Threshold. Try different thresholds
ekg_nor = ekg_ma/max(ekg_ma);
thr = mean(ekg_nor)/2; % This one was proven to be a good one after trying
% different thresholds and Moving average widndows

% 6. QRS complex. Detect the QRS value and its times
% Detect peak
signalThr = ekg_nor>thr;
QRSStart = find(diff(signalThr)>0);
QRSOver = find(diff(signalThr)<0);

% RR Time Intervals plot
figure(3)
hold on

% RR times intervals calculation and plot
timesRR=t(QRSStart(2:length(QRSStart))-QRSStart(1:length(QRSStart)-1));
plot(t(QRSStart(2:length(QRSStart))),timesRR)

% 7. Outliers search
% I wanted to show which intervals differed a lot from the actual mean of
% the timesRR and show it in the graph

outlier_max=find(timesRR>1.5*mean(timesRR));
outlier_min=find(timesRR<0.7*mean(timesRR));
scatter(t(QRSStart(outlier_max+1)),timesRR(outlier_max))
scatter(t(QRSStart(outlier_min+1)),timesRR(outlier_min))
title('RR time intervals');xlabel('Time(s)');ylabel('RR interval length(s)')    