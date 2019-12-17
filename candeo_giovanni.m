clearvars
close all
clc

% DSP PROJECT by GIOVANNI CANDEO

%PART 0: read the file containing the input signal
%read a WAVE file 
%y array containing sound samples y(n), n=1:length(y)
%Fp is the sampling frequency and should be 96 kHz
[x, Fp] = audioread('candeo_giovanni.wav');

%PART 1: compute the spectrum of the input signal
%sampling period
T = 1/Fp;
time_duration = T*size(x);
%show the original signal in time and frequency
Nx = length(x);         %length
time_x = T*(0:Nx-1);        %time samples
%TODO check why not T*fft(x)
X = fft(x);           %fft
F = 1/(Nx*T);
frequency_x = (0:Nx-1)*F;   %frequency samples

%plot signal in time domain
figure(1);
subplot(2,1,1) 
plot(time_x,x); grid;
xlabel('time[s]'); %xlim([0.1 0.2]);
title('Original audio signal in time');

%plot signal in frequency domain
subplot(2,1,2) % show frequency content in dB scale
plot(frequency_x/1e3,20*log10(abs(X))); grid; 
xlim([0 Fp/1e3]); ylim([-100 100]) 
xlabel('frequency [kHz]'); title('original audio signal in frequency')
hold on;
%shows Fp/2 
plot([1,1]*Fp/2e3,ylim,'r--');
hold off;

%fprintf('--> Press any key to listen to the frist 3 seconds of the track...  \n\n');
%pause;
%sound(y(1:3*Fp), Fp);
%pause(3);

%PART 3: find the frequencies f1 and f2 of the sinusoidal carriers
%by inspection on the spectrum
%f1 = 13000 f2 = 33600

%{
f1 and f2 are the frequencies of the two sinusoidal carriers, 
where 10000 Hz <= f1 < f2 <= 38000 Hz     
(frequencies f1 and f2 are chosen such that there is no frequency overlap 
btwn the 2 modulated components, and specifically so that f2 - f1 >= 17kHz),
and A1 and A2 are the amplitudes of the carriers. 
f2 >=
%}

%index corresponding to the frequency
i_10k = 10000/F;
i_27k = 27000/F;

%looking f1 btw 10kHz and 27kHz 
[A1,pos1] = max(abs(X(i_10k:i_27k)));
f1=frequency_x(pos1+i_10k);
disp(['f1: ' num2str(f1)]);

%looking f2 between 27kHz and 48kHZs
[A2,pos2] = max(abs(X(i_27k:Nx/2)));
f2=frequency_x(pos2+i_27k);
disp(['f2: ' num2str(f2)]);

hold on;
%shows f1 on the figure
plot([1,1]*f1/1e3,ylim,'g--');
%shows f2 on the figure
plot([1,1]*f2/1e3,ylim,'b--');
hold off;
legend('signal','Fp/2','f1','f2');

%A1,A2 are the amplitudes of the carriers
disp(['A1: ',num2str(A1)]);
disp(['A2: ',num2str(A2)]);

%lets find the coefficients to implement second order iir bandpass filters

%normalized angular frequency theta
theta = frequency_x.*(2*pi*T);
theta0(1) = 2*pi*T*f1;
theta0(2) = 2*pi*T*f2;

r = 1-pi/40;
a1 = 2*r*cos(theta0); 
a2(1) = -r^2;
a2(2) = -r^2;
b0 = (1-r)*2*sin(theta0);
%manual
%H = b0./(1-a1*exp(-1i*theta)-a2*exp(-1i*2*theta));
[H1, w1] = freqz(b0(1),[1 -a1(1) -a2(1)], 0:Fp,'whole', Fp);
[H2, w2] = freqz(b0(2),[1 -a1(2) -a2(2)], 0:Fp,'whole', Fp);
%lets extract the carriers at frequency f1 and f2
%filter frequency takes in input coefficients of a rational trans func
carrier1 = filter(b0(1),[1 -a1(1) -a2(1)],x);
carrier2 = filter(b0(2),[1 -a1(2) -a2(2)],x);
freq_carrier1 = fft(carrier1);
freq_carrier2 = fft(carrier2);

%plot filters and carriers in freq domain
figure(2); 
subplot(2,2,1);
plot(w1/1e3,20*log10(abs(H1)));
grid on; xlim([0 Fp/1e3]);xlabel('frequency [kHz]'); ylabel('|H|'); 
subplot(2,2,2);
plot(w2/1e3,20*log10(abs(H2)),'r');
grid on; xlim([0 Fp/1e3]);xlabel('frequency [kHz]'); ylabel('|H|'); 
subplot(2,2,3)
plot(frequency_x/1e3,20*log10(abs(freq_carrier1)));
grid on; xlim([0 Fp/1e3]);xlabel('frequency [kHz]');
subplot(2,2,4)
plot(frequency_x/1e3,20*log10(abs(freq_carrier2)),'r');
grid on; xlim([0 Fp/1e3]);xlabel('frequency [kHz]');



