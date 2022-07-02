clf
clear all

set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')

M=8; % the number of channels

% write a function that takes on a filter and returns the 8-channels result
% a function to return synthesis and analysis filters in the frequency domain
% a function for orth_check
% a function to display filters


%Set up M-channel DFT filter bank:
%N=2; beta=0.35;

h_mx_a = calc_filters(1);
[H_a, G_a] = uniform_filter_bank(h_mx_a);
%disp(H_a);
    
h_mx_b = calc_filters(16, 0.35);
[H_b, G_b] = uniform_filter_bank(h_mx_b);
%disp(H_b);

h_mx_c = calc_filters(24, 0.1);
[H_c, G_c] = uniform_filter_bank(h_mx_c);

inn_prod_a = orth_check(H_a);
disp(inn_prod_a);

inn_prod_b = orth_check(H_b);
disp(inn_prod_b);

inn_prod_c = orth_check(H_c);
disp(inn_prod_c);

% Generate a random Gaussian process input signal
x=randn(1,128);

% Number of DFT points 
Nfft=1024;

H_filters = {H_a, H_b, H_c};
G_filters = {G_a, G_b, G_c};

for k=1:numel(H_filters)
    H = H_filters{k};
    G = G_filters{k};
    
    y = calc_output(H,G,M,x);
    
    domega=2*pi/Nfft;
    omega=-pi:domega:pi-domega;
    
    figure(1+3*(k-1));
   
    for ind=1:M
        hf(ind,:)=abs(fftshift(fft(H(ind,:), Nfft)));
        hold on
        plot(omega,hf(ind,:));
    end
    
    box on;
        
    axis([-pi pi 0 1.1*max(hf(1,:))])
    legend('H_{0}(\omega)','H_{1}(\omega)', ...
        'H_{2}(\omega)','H_{3}(\omega)', ...
        'H_{4}(\omega)','H_{5}(\omega)', ...
        'H_{6}(\omega)','H_{7}(\omega)');
    title("Figure "+k+"("+char(96+k)+"): The Frequency Response of h_{0}[n] through h_{7}[n]");
    xlabel('\omega (radians/s)');
    ylabel('|H_{m}(\omega)|');

    yf1=abs(fftshift(fft(x,Nfft)));
    yf2=abs(fftshift(fft(y,Nfft)));

    figure(2+3*(k-1));
    plot(omega,yf1,'Linewidth',1)
    axis([-pi pi 0 max(yf1)])
    title('The magnitude of the DTFT of Gaussian random process input signal x[n] (\mu = 0, \sigma^2 = 1)');
    xlabel('\omega (radians/s)');
    ylabel('|X(\omega)|');

    figure(3+3*(k-1));
    plot(omega,yf2,'Linewidth',1)
    axis([-pi pi 0 max(yf2)])
    title('The magnitude of the DTFT of the output y[n] of the M=8 channel uniform PR filter bank');
    xlabel('\omega (radians/s)');
    ylabel('|Y(\omega)|');

end

% Takes on N and beta and returns the matrix with 8 rows
function h_mx = calc_filters(N, beta)
    n=-N:(N-1);
    
    if exist('beta','var')
        n=n+0.5;
        h=2*beta*cos((1+beta)*pi*n/2)./(pi*(1-4*beta^2*n.^2));
        h=h+sin((1-beta)*pi*n/2)./(pi*(n-4*beta^2*n.^3));
        h=h*sqrt(2);
    else
        h=[1 1]/sqrt(2);
    end

    % Obtain all necessary filters for building the system
    h0=h;
    h1=(-1).^(0:(length(n)-1)).*h;
        
    h00=zeros(1,2*length(h)); 
    h10=h00;
        
    h00(1,1:2:length(h00))=h0;
    h10(1,1:2:length(h10))=h1;
        
    h000=zeros(1,4*length(h)); 
    h100=h000;
        
    h000(1,1:4:length(h000))=h0;
    h100(1,1:4:length(h100))=h1;
    
    % A cell array    
    h_mx{1}= h0;
    h_mx{2}= h1;
    h_mx{3}= h00;
    h_mx{4}= h10;
    h_mx{5}= h000;
    h_mx{6}= h100;
    
end

function [H, G] = uniform_filter_bank(h_mx)

    h0=h_mx{1};
    h1=h_mx{2};
    h00=h_mx{3};
    h10=h_mx{4};
    h000=h_mx{5};
    h100=h_mx{6};
     
    H(1,:)=conv(conv(h0,h00), h000);
    H(2,:)=conv(conv(h0,h00), h100);
    H(3,:)=conv(conv(h0,h10), h000);
    H(4,:)=conv(conv(h0,h10), h100);
    H(5,:)=conv(conv(h1,h00), h000);
    H(6,:)=conv(conv(h1,h00), h100);
    H(7,:)=conv(conv(h1,h10), h000);
    H(8,:)=conv(conv(h1,h10), h100);

    G(1,:)=H(1,:);
    G(2,:)=-H(2,:);
    G(3,:)=-H(3,:);
    G(4,:)=H(4,:);
    G(5,:)=-H(5,:);
    G(6,:)=H(6,:);
    G(7,:)=H(7,:);
    G(8,:)=-H(8,:);

end

function inn_prod = orth_check(H)
    inn_prod = H*H';
end

function y = calc_output(H, G, M, x)

    for m=1:M
        % Pass x[n] through M filter
        W(m,:)=conv(x,H(m,:));
        % We decimate the result by the factor of M
        X(m,:)=W(m,1:M:length(W(m,:)));
    end

    for m=1:M
        % Zero inserts
        Z(m,:)=zeros(1,M*length(X(m,:)));
        Z(m,1:M:length(Z(m,:)))=X(m,:);
        % Convolution with synthesis filters
        Y(m,:)=conv(Z(m,:),G(m,:));
    end
    
    % Get the output y[n]
    y=zeros(1,length(Y(1,:)));

    for m=1:M
        y=y+Y(m,:);
    end
end
    