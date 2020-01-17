clearvars
close all
clc

% DSP PROJECT by GIOVANNI CANDEO

%PART 0: read the file containing the input signal
%read a WAVE file 
%y array containing sound samples y(n), n=1:length(y)
%Fp is the sampling frequency and should be 96 kHz
[x,Fs] = audioread('candeo_giovanni.wav');

%PART 1: compute the spectrum of the input signal
%sampling period
T = 1/Fs;
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
xlim([0 (Fs/2)/1e3]); ylim([-100 100]) 
xlabel('frequency [kHz]'); title('original audio signal in frequency')
hold on;
%shows Fp/2 
plot([1,1]*Fs/2e3,ylim,'r--');
hold off;

set(figure(1),'Units','Inches');
pos = get(figure(1),'Position');
set(figure(1),'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(figure(1),'figure1','-dpdf','-r0')

%fprintf('--> Press any key to listen to the frist 3 seconds of the track...  \n\n');
%pause;
%sound(x(1:3*Fs), Fs);
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
grid on; xlim([0 (Fs/2)/1e3]);xlabel('frequency [kHz]'); ylabel('|H|'); 
subplot(2,2,2);
plot(w2/1e3,20*log10(abs(H2)),'r');
grid on; xlim([0 (Fs/2)/1e3]);xlabel('frequency [kHz]'); ylabel('|H|'); 
subplot(2,2,3)
plot(frequency_x/1e3,20*log10(abs(freq_carrier1)));
grid on; xlim([0 (Fs/2)/1e3]);xlabel('frequency [kHz]');
subplot(2,2,4)
plot(frequency_x/1e3,20*log10(abs(freq_carrier2)),'r');
grid on; xlim([0 (Fs/2)/1e3]);xlabel('frequency [kHz]');


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

%test
%Ws = 8800/ (Fs/2); 
%[N_test,Wp] = ellipord(Wp,Ws,Rp,Rs);
%disp(['Order N = ' num2str(N_test)])

[b_lp,a_lp] = ellip(N_ellip,Rp,Rs,Wp,'low');
[H_lp, w_lp] = freqz(b_lp,a_lp,f0_lp,Fs);

figure(4)
subplot(2,1,1);
plot(w_lp/1e3,20.*log10(abs(H_lp)));
grid on; 
xlim([0 (Fs/2)/1e3]);xlabel('frequency [kHz]'); 
ylabel('|H|'); 
title('magnitude of elliptic low pass filter');

%IIR notch filter f0 = 20Hz
%high pass

f0_n = 0;
teta0_n = f0_n*2*pi*T; %will be 0 as well obv.
f3db_n = 20; %[Hz] i want a hp filter at 20 Hz
teta3db_n = f3db_n*2*pi*T; %teta 3db 
r_n = 1 - teta3db_n;

%lets find the coefficients
b_hp = [1 -2 1]; %f0 is 0 so...
a1 = 2*r_n; a2= r_n*r_n;
a_hp = [1 -a1 -a2]; %from formula 
[H_notch, w_notch] = freqz(b_hp,a_hp,2048,'whole',Fs); %2048 arbitrary i guess

figure(5)
subplot(2,1,1);
plot(w_notch/1e3, 20.*log10(abs(H_notch))); 
grid on; 
xlim([0 (Fs/2)/1e3]);xlabel('frequency [kHz]');
ylabel('|H|');
title('magnitude of second order IIR notch filter');

%filter the signals to remove distortions
%signal1
signal1_hp = filter(a_hp,b_hp,signal1);
signal1_clean = filter(a_lp,b_lp,signal1_hp);
signal1_clean = signal1./A1;
SIGNAL1_clean = fft(signal1_clean(1:Nx));

%signal2
signal2_hp = filter(a_hp,b_hp,signal2);
signal2_clean = filter(a_lp,b_lp,signal2_hp);
signal2_clean = signal2./A1;
SIGNAL2_clean = fft(signal2_clean(1:Nx));

%audiowrite('candeo_giovanni.wav',[signal1_clean,signal2_clean],Fs);

%plots
figure(6)
subplot(3,2,1);
plot(time_x,signal1);
xlabel('time[s]'); 
title('signal 1 in time');
subplot(3,2,3);
plot(time_x,signal1_clean);
xlabel('time[s]'); 
title('clean signal 1 in time');
subplot(3,2,5);
plot(frequency_x/1e3,20*log10(abs(SIGNAL1_clean)));
xlim([0 (Fs/2)/1e3]); xlabel('frequency [KHz]');
title('magnitude of signal1 clean');

subplot(3,2,2);
plot(time_x,signal2);
xlabel('time[s]'); 
title('signal 2 in time');
subplot(3,2,4);
plot(time_x,signal2_clean);
xlabel('time[s]'); 
title('clean signal 2 in time');
subplot(3,2,6);
plot(frequency_x/1e3,20*log10(abs(SIGNAL2_clean)));
xlim([0 (Fs/2)/1e3]); xlabel('frequency [KHz]');
title('magnitude of signal2 clean');

