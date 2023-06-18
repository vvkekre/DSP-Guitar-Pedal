%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENGR 362 Project Template File
% student name:
% student number:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all, clear, clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import audio data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = importdata('ENGR_362_guitar_Fs_is_48000_Hz.txt');
Fs = 48000;                             % sampling freq

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Play sound of raw audio data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% soundsc(y,Fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graph in time domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ts = 1/Fs;                              % sampling period       
Length_y = length(y(:,1));              % length of signal
time = (0:Length_y-1)*Ts;               % time vector
figure,plot(time,y), axis tight
title('unfiltered: y(t) vs t')                      % labels
xlabel('time')                         % labels
ylabel('y(t)')                          % labels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graph with DFT/FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = fft(y);                             % Discrete Fourier transform
F1 = abs(Y/Length_y);                   % frequency
F2 = F1(1:Length_y/2+1);                % half of frequency
F2(2:end-1) = 2*F2(2:end-1);            % Discrete Fourier transform
f = Fs*(0:(Length_y/2))/Length_y;       % freq vector [Hz]
f_kHz = f/1000;                         % freq vector [kHz]

figure(1) ;                             
subplot(1,2,1)
plot(time,y)
axis tight
title('Unfiltered: y(t) vs t')           % labels
xlabel('time')                         % labels
ylabel('y(t)')                          % labels
subplot(1,2,2)
plot(f_kHz,F2) 
axis([0 5  0 0.06])
title('Unfiltered: Y(F) vs F ')                           % labels
xlabel('F (kHz)')                             % labels
ylabel('Y(F)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D major chord frequencies for notes D3, A3, D4, F#4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D3 = 146.83;                            % freq of note D3 [Hz]
D3_int = round(D3/max(f)*length(f))     % associated integer to above freq
A3 = 220.00;                            % freq of note A3 [Hz]
A3_int = round(A3/max(f)*length(f))     % associated integer to above freq
D4 = 293.66;                            % freq of note D4 [Hz]
D4_int = round(D4/max(f)*length(f))     % associated integer to above freq
F_sharp_4 = 369.99;                     % freq of note F#4 [Hz]
F_sharp_4_int = ...
    round(F_sharp_4/max(f)*length(f))   % associated integer to above freq

note_freq = [D3 A3 D4 F_sharp_4];       % vector of all note freqs
note_freq_int = ...
    [D3_int A3_int D4_int F_sharp_4_int]; % vector of all int note freqs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop to apply filter bank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please continue writing code from here...

delta = 2.5;

for i = 1:length(note_freq)
    center_frequency = round(note_freq(i));
    [m1,n1] = cheby1(7, 1.25 , (center_frequency-delta)/(Fs/2), 'High');
    [m2,n2] = cheby1(7, 1.25 , (center_frequency+delta)/(Fs/2), 'Low');

    if(i == 1)                                      
        D3_loop = filter(m2, n2, filter(m1, n1, y));
    end

    if(i == 2)                                      
        A3_loop = filter(m2, n2, filter(m1, n1, y));
    end

    if(i == 3)                                      
        D4_loop = filter(m2, n2, filter(m1, n1, y));
    end

    if(i == 4)                                      
        F_sharp_4_loop = filter(m2, n2, filter(m1, n1, y));
    end
end

array = [D3_loop,A3_loop, D4_loop,F_sharp_4_loop];
sum = D3_loop + A3_loop + D4_loop + F_sharp_4_loop;

%ft of signal: D3
D3_ft = fft(D3_loop);
D3ft1 = abs(D3_ft/Length_y);
D3ft2 = D3ft1(1: Length_y/2+1);
D3ft2(2:end-1) = 2*D3ft2(2:end-1);

%ft of signal: A3
A3_ft = fft(A3_loop);
A3ft1 = abs(A3_ft/Length_y);
A3ft2 = A3ft1(1: Length_y/2+1);
A3ft2(2:end-1) = 2*A3ft2(2:end-1);

%ft of signal: D4
D4_ft = fft(D4_loop);
D4ft1 = abs(D4_ft/Length_y);
D4ft2 = D4ft1(1: Length_y/2+1);
D4ft2(2:end-1) = 2*D4ft2(2:end-1);

%ft of signal: F_sharp_4 
F_sharp_4_ft = fft(F_sharp_4_loop);
F_sharp_4ft1 = abs(F_sharp_4_ft/Length_y);
F_sharp_4ft2 = F_sharp_4ft1(1: Length_y/2+1);
F_sharp_4ft2(2:end-1) = 2*F_sharp_4ft2(2:end-1);

% scaling the signals in f domain 

mamp = [max(D3ft2), max(A3ft2),max(D4ft2),max(F_sharp_4ft2)] %maximum amplitude
nf = max(mamp)./mamp  %norm factor

%normalizing

D3ft2 = D3ft2.*nf(1)
A3ft2 = A3ft2.*nf(2)
D4ft2 = D4ft2.*nf(3)
F_sharp_4ft2 = F_sharp_4ft2.*nf(4)

total_fft = D3ft2 + A3ft2 + D4ft2 + F_sharp_4ft2
t_f_n = D3_loop*nf(1) + A3_loop*nf(2) + D4_loop*nf(3) + F_sharp_4_loop(4)  %total filter norm


%Plot for D3

ylabel('dB')
subplot(4,1,1);    
[xD3,yD3] = freqz(D3_loop);
plot(yD3/pi,20*log10(abs(xD3)))


%Plot for A3
ylabel('dB')
subplot(4,1,2);
[xA3,yA3] = freqz(A3_loop); 
plot(yA3/pi,20*log10(abs(xA3)))

%Plot for D4
ylabel('dB')
subplot(4,1,3);
[xD4,yD4] = freqz(D4_loop);
plot(yD4/pi,20*log10(abs(xD4)))

%Plot for F_sharp_4
ylabel('dB')
subplot(4,1,4);
[xF_sharp_4,yF_sharp_4] = freqz(F_sharp_4_loop);
plot(yF_sharp_4/pi,20*log10(abs(xF_sharp_4)))


%FFT plot of all 4 frequencies in one graph.

figure(2);
hold on
plot(f_kHz,D3ft2) 
plot(f_kHz,A3ft2)
plot(f_kHz,D4ft2)
plot(f_kHz,F_sharp_4ft2)
axis([0 0.5 0 0.05])
title('Filtered FFT vs F')                           
xlabel('Frequency(kHz)')                       
ylabel('Y(F)')
hold off


% Filtered signal 
figure(3);
subplot(2,1,1);
plot(time, sum);
axis([0 max(time) -0.7 0.7]);
title('Filtered: y(t) vs time')                           
xlabel('time')                       
ylabel('y(t)')
subplot(2,1,2);
plot(time, t_f_n);
axis([0 max(time) -1.5 1.5]);
title('Normal: y(t) vs time')                           
xlabel('time')                       
ylabel('y(t)')

soundsc(t_f_n,Fs)










