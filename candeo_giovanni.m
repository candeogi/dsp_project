clearvars
close all
clc

% DSP PROJECT by Candeo Giovanni
% Student ID 1206150, Audio Sample n.033.

%read the file containing the input signal
%y array containing sound samples y(n), n=1:length(y)
%Fp is the sampling frequency and should be 96 kHz
[x,Fs] = audioread('signal_033.wav');

%compute the spectrum of the input signal
T = 1/Fs; %sampling period
time_duration = T*size(x);

%show the original signal in time and frequency
Nx = length(x);%length
time_x = T*(0:Nx-1);%time samples
X = fft(x);       
F = 1/(Nx*T);
frequency_x = (0:Nx-1)*F;   %frequency samples

%plot signal in time and frequency domain
figure(1);
subplot(3,1,1) 
plot(time_x,x); grid;
xlabel('time[s]');ylabel('x(nT)');
title('Original audio signal in time');

subplot(3,1,2)
plot(frequency_x/1e3,abs(X)); grid; 
xlim([0 (Fs/2)/1e3]); 
xlabel('frequency [kHz]');ylabel('|X(f)|');
title('original audio signal in frequency')

%plot signal in frequency domain
subplot(3,1,3) % show frequency content in dB scale
plot(frequency_x/1e3,20*log10(abs(X))); grid; 
xlim([0 (Fs/2)/1e3]); ylim([-100 100]) 
xlabel('frequency [kHz]'); ylabel('|X(f)|');
ylabel('db');
title('original audio signal in frequency (db)')


%find the frequencies f1 and f2 of the sinusoidal carriers

%{
f1 and f2 are the frequencies of the two sinusoidal carriers, 
where 10000 Hz <= f1 < f2 <= 38000 Hz     
(frequencies f1 and f2 are chosen such that there is no frequency overlap 
between the 2 modulated components, and specifically so that f2 - f1 >= 17kHz),
A1 and A2 are the amplitudes of the carriers. 
%}

%index corresponding to the frequency
i_10k = 10000/F;
i_27k = 27000/F;

%looking f1 btw 10kHz and 27kHz 
[A1,pos1] = max(abs(X(i_10k:i_27k)));
A1 = A1/(Nx/2);
f1=frequency_x(pos1+i_10k);
%disp(['f1: ' num2str(f1)]);

%looking f2 between 27kHz and 48kHZs
[A2,pos2] = max(abs(X(i_27k:Nx/2)));
A2 = A2/(Nx/2);
f2=frequency_x(pos2+i_27k);
%disp(['f2: ' num2str(f2)]);

hold on;
%shows f1 on the figure
plot([1,1]*f1/1e3,ylim,'g--');
%shows f2 on the figure
plot([1,1]*f2/1e3,ylim,'r--');
hold off;
legend('signal','f1','f2');


%lets find the coefficients to implement second order iir bandpass filters

%normalized angular frequency theta
theta = frequency_x.*(2*pi*T);
theta0(1) = 2*pi*T*f1;
theta0(2) = 2*pi*T*f2;

%suppose DELTA f3db = 10Hz
%very narrow bandwitdh
DELTA_teta3db = 2*pi*10/Fs;
delta = DELTA_teta3db/2;
r = 1-delta;
a1 = 2*r*cos(theta0); 
a2(1) = -r^2;
a2(2) = -r^2;
b0 = (1-r)*2*sin(theta0);
%H = b0./(1-a1*exp(-1i*theta)-a2*exp(-1i*2*theta)); %manual
[H1, w1] = freqz(b0(1),[1 -a1(1) -a2(1)], 0:Fs,'whole', Fs);
[H2, w2] = freqz(b0(2),[1 -a1(2) -a2(2)], 0:Fs,'whole', Fs);
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
grid on; xlim([0 (Fs/2)/1e3]);xlabel('frequency [kHz]');
title('magnitude (db) bp filter at f1');
subplot(2,2,2);
plot(w2/1e3,20*log10(abs(H2)),'r');
grid on; xlim([0 (Fs/2)/1e3]);xlabel('frequency [kHz]');
title('magnitude (db) bp filter at f2');
subplot(2,2,3)
plot(frequency_x/1e3,20*log10(abs(freq_carrier1)));
grid on; xlim([0 (Fs/2)/1e3]);xlabel('frequency [kHz]');
title('magnitude (db) of the first carrier');
subplot(2,2,4)
plot(frequency_x/1e3,20*log10(abs(freq_carrier2)),'r');
grid on; xlim([0 (Fs/2)/1e3]);xlabel('frequency [kHz]'); 
title('magnitude (db) of the second carrier');

%demodulation
signal1 = x.*carrier1;
SIGNAL1 = fft(signal1);
signal2 = x.*carrier2;
SIGNAL2 = fft(signal2);

%fprintf('--> Press any key to listen signal 1\n\n');
%pause;
%sound(signal1(1:5*Fs), Fs);
%pause(3);

%fprintf('--> Press any key to listen signal 2\n\n');
%pause;
%sound(signal2(1:5*Fs), Fs);
%pause(3);

%
%Filters to remove distorsion
%[20-8000 Hz]

f0_lp = 8000; %[Hz] 

%{
Minimax
% number of samples is N+1
N = 100; % must be an even number
% limit frequencies
al = 0.05; % transition bandwidth in percentage
fp = f0_lp*(1-al); % pass band upper limit
fs = f0_lp*(1+al); % stop band lower limit
err_lim = 0.0001; % -80 dB attenuation
[NN,Fo,Ao,W] = firpmord([fp fs],[1 0],[1 1]*err_lim,Fs);
disp(['firpmord suggests a filter of order ' num2str(NN) ...
      ' for guaranteeing a ' num2str(20*log10(err_lim)) ' dB error'])

% define filter
h0 = firpm(N,Fo,Ao,W)/T;
t = T*(-N/2:N/2);
% behaviour in frequency obtained by use of function freqz
[H0,ff] = freqz(h0,1,8*(N+1),Fs);
H0 = T*H0; % normalization factor
% show results
figure(3)
subplot(2,1,1)
stem(t,h0); grid; title('FIR - remez approach - time domain')
subplot(2,1,2)
plot(ff,20*log10(abs(H0))); grid; xlim([0 Fs/2]); ylim([-80, 5]); 
title('frequency domain')
hold on; plot([1,1]*fp,ylim,'r--'); plot([1,1]*fs,ylim,'r--'); hold off;
%}

%IIR low pass using elliptic
%we want stop band attenuation of 80db
N_ellip = 10;
Rp = 0.1; %passband ripple
Rs = 80; %stopband attenuation
Wp = f0_lp/ (Fs/2); %passband edge frequency

[b_lp,a_lp] = ellip(N_ellip,Rp,Rs,Wp,'low');
[H_lp, w_lp] = freqz(b_lp,a_lp,f0_lp,Fs);

figure(3)
subplot(2,1,1);
plot(w_lp/1e3,20.*log10(abs(H_lp)));
grid on; 
xlim([0 (Fs/2)/1e3]);xlabel('frequency [kHz]'); 
ylabel('|H|'); 
title('magnitude (db) of elliptic low pass filter');
subplot(2,1,2);
plot(w_lp/1e3, angle(H_lp)*360);
grid on;
xlim([0 (Fs/2)/1e3]);xlabel('frequency [kHz]');
ylabel('angle(H)');
title('angle of elliptic low pass filter');


%IIR notch filter f0 = 20Hz
%high pass
f3db_hp = 20; %[Hz] i want a hp filter at 20 Hz
DELTA3db_hp = 2*pi*f3db_hp/Fs; %DELTA 3db 
r_hp = 1 - DELTA3db_hp/2;

%lets find the coefficients
b_hp = [1 -2 1]; %f0 is 0 so...
a1 = 2*r_hp; a2= r_hp*r_hp;
a_hp = [1 -a1 a2]; %from formula 
[H_hp, w_hp] = freqz(b_hp,a_hp,2048,'whole',Fs); 

figure(4)
subplot(2,1,1);
plot(w_hp/1e3, 20.*log10(abs(H_hp))); 
grid on; 
xlim([0 (Fs)/1e3]);xlabel('frequency [kHz]');
ylabel('|H|');
title('magnitude (db) of second order IIR notch filter');
subplot(2,1,2);
plot(w_hp/1e3, angle(H_hp)*360);
grid on;
xlim([0 (Fs)/1e3]);xlabel('frequency [kHz]');
ylabel('angle(H)');
title('angle of second order IIR notch filter');

%filter the signals to remove distortions
%signal1
signal1_hp = filter(a_hp,b_hp,signal1);
signal1_clean = filter(a_lp,b_lp,signal1_hp);
signal1_clean = signal1./A1;
SIGNAL1_clean = fft(signal1_clean(1:Nx));

%signal2
signal2_hp = filter(a_hp,b_hp,signal2);
signal2_clean = filter(a_lp,b_lp,signal2_hp);
signal2_clean = signal2./A2;
SIGNAL2_clean = fft(signal2_clean(1:Nx));

%sound(signal1_clean,Fs);
%sound(signal2_clean,Fs);
audiowrite('candeo_giovanni.wav',[signal1_clean,signal2_clean],Fs);

%plots
figure(5)
subplot(3,2,1);
plot(time_x,signal1);grid on;
xlabel('time[s]'); 
title('signal 1 in time');
subplot(3,2,3);
plot(time_x,signal1_clean);grid on;
xlabel('time[s]'); 
title('clean signal 1 in time');
subplot(3,2,5);
plot(frequency_x/1e3,20*log10(abs(SIGNAL1_clean))); grid on;
xlim([0 (Fs/2)/1e3]); xlabel('frequency [KHz]');
title('magnitude of signal1 clean');

subplot(3,2,2);
plot(time_x,signal2);grid on;
xlabel('time[s]'); 
title('signal 2 in time');
subplot(3,2,4);
plot(time_x,signal2_clean);grid on;
xlabel('time[s]'); 
title('clean signal 2 in time');
subplot(3,2,6);
plot(frequency_x/1e3,20*log10(abs(SIGNAL2_clean)));grid on;
xlim([0 (Fs/2)/1e3]); xlabel('frequency [KHz]');
title('magnitude of signal2 clean');

%prints
%print(figure(1),'-bestfit','figure1','-dpdf');
%print(figure(2),'-bestfit','figure2','-dpdf');
%print(figure(3),'-bestfit','figure3','-dpdf');
%print(figure(4),'-bestfit','figure4','-dpdf');
%print(figure(5),'-bestfit','figure5','-dpdf');


