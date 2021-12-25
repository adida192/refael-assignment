clc;
close all;
clear all;

PRI = 5e-6 ; %5*10^-6 
PW = 2e-6 ; %2*10^-6 ;
CPI = 0.5e-3 ; %0.5*10^-3 ;
B = 100e6 ; %100*10^6 ; 
fc = 40e9 ; %40*10^9 ; 
R0 = 400 ;
v_target = 100;
c = 3e8;
delay = 2*R0/c;
Ws = 10*B ; % Stand in Nyquist
Ts = 2*pi / Ws;
Fs = 1/Ts ;
Num_Pulses = CPI/PRI;
t = 0: Ts: Num_Pulses*PRI ;


W = zeros(1, length(t));
W(t<PW) = 1;
h_t = exp(1i*pi*(B/PW)*(t-PW/2).^2) .* W;
size_Tx = sum(W==1);
size_PRI = sum (t<PRI);
for n  = 1:Num_Pulses-1 
    strat_index = n*size_PRI + 1;
    last_index = n*size_PRI + size_Tx;
    h_t(strat_index : last_index )= h_t(1:size_Tx);
end

% Single Pulse in Time
figure(1);
plot(10^6*t(1:size_PRI),h_t(1:size_PRI),'LineWidth',1)
grid on;
xlabel('t(usec)','FontSize',12)
ylabel('h(t)','FontSize',12)
title('Single Pulse','FontSize',12)

% Pulse in Frequency domain
figure(2)
H = fftshift(fft(h_t(1:size_Tx)));
fshift = (-size_Tx/2:size_Tx/2-1)*(Fs/size_Tx); % zero-centered frequency range
% powershift = abs(H).^2/size_TR;     % zero-centered power
%plot(fshift,powershift)
plot(10^-6 * fshift,abs(H),'LineWidth',1)
grid on;
xlabel('f(MHz)','FontSize',12)
ylabel('|h(f)|','FontSize',12)
title('Single Pulse in Frequency Domain','FontSize',12)

% Autocorrelation for single pulse
figure(3)
autocorr(h_t(1:size_Tx))

% Pulse Train
figure(4);
plot(10^3*t , h_t, 'LineWidth',0.1)
grid on;
xlabel('t(msec)','FontSize',12)
ylabel('h(t)','FontSize',12)
title('Pulse Train','FontSize',12)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Stx = h_t.*exp(1i*2*pi*fc*t);



