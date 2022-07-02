clear all
clf

set(0,'defaultTextInterpreter','latex');

% Generate Gaussian noise
num_samples = 200;
v = randn(1, num_samples);

% Shift registers
reg_4 = [1 0 0 0];
reg_6 = [1 0 0 0 0 0];
reg_7 = [1 0 0 0 0 0 0];

regs = {reg_4, reg_6, reg_7};

% Set of parameters a_{2} and D_{2} and also time delay D_{1}
params = [1,22; 1,21; -1,21];
D1 = 20;

for k=1:numel(regs)
    reg = regs{k};
    % Generate the input sequence of dimension M
    x = gen_input(reg);
    M = length(x);
    x_pad = [x, zeros(1, num_samples-M)]

    for i=1:length(params)
        
        % Pick the amplitude a_{2} and the delay D_{2}
        % from the set of parameters
        a2 = params(i,1);
        D2 = params(i,2);
        % Calculate the output and the cross-correlation 
        % between the output and the input
        [y, ryx] = calc(x, a2, D1, D2, v);
        
        % Plot the results
        fig=figure(i+3*(k-1));
        
        subplot(3, 1, 1);
        plot(0:num_samples-1,x_pad, 'Linewidth', 2)
        title("The input sequence x[n] for M="+M)
        xlabel("n");
        ylabel("x[n]");
        
        subplot(3, 1, 2);
        plot(0:num_samples-1, y, 'Linewidth', 2);
        title(["The received signal y[n] with additive noise $\mathcal{N}$(0,1)"... 
        "for $a_{2}$="+a2+", $D_{2}$="+D2]);
        xlabel('n');
        ylabel('y[n]');
        
        subplot(3, 1, 3);
        plot(0:59, ryx(M:M + 59), 'Linewidth', 2);
        title(["The cross-correlation $r_{yx}[l]$ between the output y[l]"... 
        "and the input x[l] for parameters $a_{2}$="+a2+", $D_{2}$="+D2]);
        xlabel("$0 \leq l \leq 59$");
        ylabel("$r_{yx}[l]$");
        xline(D1, 'r', 'Linewidth', 1.5);
        xline(D2, 'r', 'Linewidth', 1.5);
        
        saveas(fig, sprintf('fig%d.png', i+3*(k-1)))

    end;
end;    

function x = gen_input(reg)
    N = 2^length(reg)-1;
    for ri=1:N
        x(ri)=reg(1,end);
        reg(2:end)=reg(1:end-1);
        reg(1,1)=rem((reg(1,1)+x(1,ri)),2);
    end
    % Transform 0s and 1s to -1s and 1s 
    x = 2*x-1;
end
    
function [y, ryx] = calc(x, a2, D1, D2, v)
    x_D1 = [zeros(1, D1), x, zeros(1, length(v)-length(x)-D1)];
    x_D2 = [zeros(1, D2), x, zeros(1, length(v)-length(x)-D2)];
    y = x_D1 + a2 .* x_D2 + v;
    ryx = conv(y,x(end:-1:1));
end