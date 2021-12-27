clear all;
shift = 3;
h_t = [1,2,3,4,5];
x = circshift(h_t,shift);
N = numel(x);
ix = (1:N) - shift;
tf = ix < 1 | ix > N;
x(tf) = 0 ;

% syms f(x)
% x = t;
% y = h_t;
% [~, index] = sort(x);
% F = griddedInterpolant(x(index), y(index));
% test = scaled*t - delay;
% test(1:425) = zeros(1,425);
% check = F(test);
% check(1:425) = zeros(1,425);
% ff = check.*exp(1i*2*pi*fc*v_target*t/c);
% figure(22)
% plot(t(1,size_PRI),check(1,size_PRI));