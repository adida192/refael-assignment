clc;
close all;
clear all;
%%%% Parameters  %%%%

PRI = 5e-6 ; %5*10^-6 
PW = 2e-6 ; %2*10^-6 ;
CPI = 0.5e-3 ; %0.5*10^-3 ;
B = 100e6 ; %100*10^6 ; 
fc = 40e9 ; %40*10^9 ; 
R0 = 400 ;
v_target = 100;
c = 3e8;

Ws = 10*B ; % Stand in Nyquist
Ts = (8/3)*10^-9; %2*pi / Ws;

Fs = 1/Ts ;
Num_PRI = CPI/PRI;
t = 0: Ts: CPI +(Num_PRI-1)*Ts; %- Ts ;

%%%% Pulse train %%%%

W = zeros(1, length(t));
W(t<PW) = 1;
h_t_pulse = exp(1i*pi*(B/PW)*(t-PW/2).^2) .* W;
size_Tx = sum(W==1);
size_PRI = sum (t<PRI);
h_t = zeros(1, length(t));
for n  = 0: Num_PRI-1 
    start_index = n*size_PRI + 1;
    last_index = n*size_PRI + size_Tx;
    h_t(start_index : last_index )= h_t_pulse(1:size_Tx);
end
 
% Single Pulse in Time
figure(1);
plot(10^6*t(1:size_PRI),abs(h_t(1:size_PRI)),'LineWidth',1)
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
plot(10^3*t , abs(h_t), 'LineWidth',0.1)
grid on;
xlabel('t(msec)','FontSize',12)
ylabel('h(t)','FontSize',12)
title('Pulse Train','FontSize',12)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_TX = h_t.*exp(1i*2*pi*fc*t);
delay = 2*R0/c;
%scaled = 1-2*v_target/c;
shift = find(t==(8/3)*10^-6);
h_shifted = circshift(h_t,shift);
N = numel(h_shifted);
ix = (1:N) - shift;
tf = ix < 1 | ix > N;
h_shifted(tf) = 0 ;
g = exp(-1i*2*pi*fc*2*R0/c);
S_RX = h_shifted .*...
    exp(-1i*2*pi*fc*2*v_target*t/c);

PRI_matrix = reshape(S_RX,[Num_PRI,size_PRI]);
r_PRI_matrix = PRI_matrix(:,size_Tx+1:end);
corr_function = zeros(size(r_PRI_matrix ));
h = h_t(1:size_Tx);
b = conj(h(end:-1:1));
for i = 1:Num_PRI
    corr_function(i,:) = filter(b,1,r_PRI_matrix(i,:));
end

dopler_matrix = zeros(size(r_PRI_matrix ));
for j = 1:length(r_PRI_matrix(1,:))
    dopler_matrix(:,j) = fftshift(fft(corr_function(:,j)));
end

imagesc(100,50,mag2db(abs(dopler_matrix)));
colorbar;


% figure(33);
% % plot(10^3*t(1:size_PRI) , h_shifted(1:size_PRI), '*')
% % hold on;
% plot(10^3*t , S_RX, '-')
