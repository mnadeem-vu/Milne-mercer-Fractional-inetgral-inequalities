% FIG_EXAMPLED  Reproduces exampleD.png (s-convex, s=1/2, f(x)=e^x).
%
%   f(x) = e^x, [a,b] = [0,1], lambda = 1/3, mu = 2/3, h(t) = t^(1/2)
%   (s-convex, s = 1/2).
%   Requires LHS_exact.m, Kfun.m, RHS_thm1.m on the MATLAB path (same folder).

clear; close all; clc;

f      = @(x) exp(x);
fprime = @(x) abs(exp(x));
h      = @(t) t.^0.5;  % s-convex, s = 1/2

a = 0; b = 1; lambda = 1/3; mu = 2/3;

beta_grid = linspace(0.05, 0.95, 40);
LHS = zeros(size(beta_grid));
RHS = zeros(size(beta_grid));

for k = 1:numel(beta_grid)
    beta = beta_grid(k);
    LHS(k) = abs(LHS_exact(f, a, b, lambda, mu, beta));
    RHS(k) = RHS_thm1(fprime, a, b, lambda, mu, beta, h);
end
Gap = RHS - LHS;

beta_table = [0.2 0.4 0.6 0.8];
fprintf('  beta      LHS         RHS         Gap\n');
for beta = beta_table
    L = abs(LHS_exact(f, a, b, lambda, mu, beta));
    R = RHS_thm1(fprime, a, b, lambda, mu, beta, h);
    fprintf('  %.1f   %10.6f  %10.6f  %10.6f\n', beta, L, R, R-L);
end

figure;
plot(beta_grid, LHS, 'o-', 'DisplayName', 'LHS'); hold on;
plot(beta_grid, RHS, 's-', 'DisplayName', 'RHS');
plot(beta_grid, Gap, '^-', 'DisplayName', 'Gap');
xlabel('\beta'); ylabel('Value');
title('Sensitivity: f(x)=e^x, h(t)=t^{0.5} (s-convex)');
legend('Location', 'best'); grid on;
