clearvars
close all
clc

% DSP PROJECT by GIOVANNI CANDEO

%PART 0: read the file containing the input signal
%read a WAVE file 
%y array containing sound samples y(n), n=1:length(y)
%Fp is the sampling frequency and should be 96 kHz
[y, Fp] = audioread('candeo_giovanni.wav');

%PART 1: compute the spectrum of the input signal
%sampling period
T = 1/Fp;
time_duration = T*size(y);
%show the original signal in time and frequency
Ny = length(y);         %length
ty = T*(0:Ny-1);        %time samples
Y = T*fft(y);           %fft
F = 1/(Ny*T)
fy = (0:Ny-1)*F;   %frequency samples
%plot
figure
subplot(2,1,1) 
plot(ty,y); grid;
xlabel('time[s]'); %xlim([0.1 0.2]);
title('Original audio signal in time');
subplot(2,1,2) % show frequency content in dB scale
plot(fy/1e3,20*log10(abs(Y))); grid; 
xlim([0 Fp/1e3]); ylim([-200 50]) 
xlabel('frequency [kHz]'); title('original audio signal in frequency')
hold on;
%shows Fp/2 
plot([1,1]*Fp/2e3,ylim,'r--');
hold off;


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
[A1,pos1] = max(abs(Y(i_10k:i_27k)));
f1=fy(pos1+i_10k);
disp(['f1: ' num2str(f1)]);

%looking f2 between 27kHz and 48kHZs
[A2,pos2] = max(abs(Y(i_27k:Ny/2)));
f2=fy(pos2+i_27k);
disp(['f2: ' num2str(f2)]);

hold on;
%shows f1
plot([1,1]*f1/1e3,ylim,'g--');
%shows f2
plot([1,1]*f2/1e3,ylim,'b--');
hold off;
%test123=fy(27000/F)
legend('signal','Fp/2','f1','f2');

%A1,A2 are the amplitudes of the carriers
disp(['A1: ',num2str(A1)]);
disp(['A2: ',num2str(A2)]);



