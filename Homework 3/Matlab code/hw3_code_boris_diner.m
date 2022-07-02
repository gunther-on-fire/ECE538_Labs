set(0,'defaultaxesfontsize',12);
clear
clf

% Problem 7.29
% Parameters
a=0.8;
Nfft=1024;

% The DTFT of x[n]
wd=linspace(0,2*pi-2*pi/Nfft,Nfft);
Xw = (1-a^2)./(1-2*a*cos(wd)+a^2);

% Reconstructed versions of the DTFT above
% for different number of its samples N
Xr_21 = rec_dtft(a,21,Nfft,wd);
Xr_101 = rec_dtft(a,101,Nfft,wd);

fig1=figure(1);
% The reconstructed spectrum from the formula for N=21
plot(wd,abs(Xr_21),'color','#0072BD','Linewidth',3);
hold on
% The DTFT
plot(wd,abs(Xw),'r--', 'Linewidth',3);
axis([0 2*pi 0 max(abs(Xw))]);
legend('|X_{r}(\omega)|','|X(\omega)|');
xlabel('\omega');
saveas(fig1, sprintf('lab3fig1.png'));

fig2=figure(2);
% The reconstructed spectrum from the formula for N=101
plot(wd,abs(Xr_101),'color','#0072BD','Linewidth',3);
hold on
% The DTFT
plot(wd,abs(Xw), 'r--', 'Linewidth',3);
axis([0 2*pi 0 max(abs(Xw))]);
legend('|X_{r}(\omega)|','|X(\omega)|');
xlabel('\omega');
saveas(fig2, sprintf('lab3fig2.png'));

% Time domain aliasing
[x_t_21, x_hat_21, x_a_21] = time_domain_aliasing(a,21);
[x_t_101, x_hat_101, x_a_101] = time_domain_aliasing(a,101);

n_21=0:20;

fig3=figure(3);
plot(n_21,x_t_21,'color','#4DBEEE','Linewidth',4);
hold on
plot(n_21,x_hat_21,'color','#0072BD', 'Linewidth',4);
hold on
plot(n_21,x_a_21,'^','MarkerFaceColor','#7E2F8E','MarkerEdgeColor','#7E2F8E');
hold off
legend('x_t[n]','\hat{x}[n]','x_a[n]');
xlabel('n');
saveas(fig3, sprintf('lab3fig3.png'));

n_101=0:100;

fig4=figure(4);
plot(n_101,x_t_101,'color','#4DBEEE','Linewidth',4);
hold on
plot(n_101,x_hat_101,'-.','color','#0072BD', 'Linewidth',4);
hold on
plot(n_101,x_t_101,'^','MarkerFaceColor','#7E2F8E','MarkerEdgeColor','#7E2F8E');
hold off
legend('x_t[n]','x_{hat}[n]','x_a[n]');
xlabel('n');
saveas(fig4, sprintf('lab3fig4.png'));

% Problem 7.30
% Frequencies
f1 = 1/128; 
f2 = 5/128; 
fc = 50/128;

% a)
n = 0:255;
x = cos(2*pi*f1*n) + cos(2*pi*f2*n);
x_c = cos(2*pi*fc*n);
x_am = x.*cos(2*pi*fc*n);

fig5=figure(5);
subplot(3,1,1);
plot(n, x);
xlabel('n');
xlim([0 n(end)]);
title('x[n]=cos(2\pi f_{1}n)+cos(2\pi f_{2} n)');
subplot(3,1,2);
plot(n,x_c);
xlabel('n');
xlim([0 n(end)]);
title('x_c[n]=cos(2\pif_cn)');
subplot(3,1,3);
plot(n,x_am);
xlabel('n');
xlim([0 n(end)]);
title('x[n]=x[n]\cdot x_c[n]');
saveas(fig5, sprintf('lab3fig5.png'));


% b), c), and d)
X_am_128 = fft(x_am(1:128),128);
X_am_128_t = fft(x_am(1:100),128);
X_am_256 = fft(x_am(1:180),256);

fig6=figure(6);
plot((0:127),abs(X_am_128));
title('|X_{am}[k]|, 0 \leq n \leq 127');
xlabel('k');
xlim([0 128]);
saveas(fig6, sprintf('lab3fig6.png'));

fig7=figure(7);
plot((0:127),abs(X_am_128_t));
title('|X_{am}[k]|, 0 \leq n \leq 99');
xlabel('k');
xlim([0 128]);
saveas(fig7, sprintf('lab3fig7.png'));

fig8=figure(8);
plot((0:255),abs(X_am_256));
title('|X_{am}[k]|, 0 \leq n \leq 179');
xlabel('k');
xlim([0 256]);
saveas(fig8, sprintf('lab3fig8.png'));

function Xr = rec_dtft(a,N,Nfft,wd)

D=(N-1)/2;

% initialize the reconstructed DTFT
Xr=zeros(1,Nfft);

for k=0:N-1
    wk=2*pi*k/N;
    wvec(1,k+1)=wk;
    % Sampling the "true" DTFT at N equispaced frequencies
    Xk=(1-a^2)./((1-2*a*cos(wk)+a^2)).*exp(-1i*wk*D);
    Xvec(1,k+1)=Xk;
    % Exercising the reconstruction formula
    Xr=Xr+Xk*(sin((wd-wk)*N/2).*exp(-1i*(wd-wk)*(N-1)/2))./(N*sin((wd-wk)/2));
end

end

function [x_t, x_hat, x_a] = time_domain_aliasing(a,N)

D = (N-1)/2;

% Original x[n]
n = 0:(N-1);
x_t = a.^(abs(n-D));
% Construct x_hat[n]
k = 0:N-1;
wk = 2*pi*k/N;
X_wk = (1-a^2)./(1-2*a*cos(wk)+a^2).*exp((-1i)*D*wk);
x_hat = real(ifft(X_wk,N));
% Construct x_a[n]
x_a = a.^(abs(n+N-D)) + a.^(abs(n-D)) + a.^(abs(n-N-D));

end