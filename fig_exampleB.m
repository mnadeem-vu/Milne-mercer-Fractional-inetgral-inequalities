% FIG_EXAMPLEB  Reproduces exampleB.png (P-convex case, f(x)=(1+x)^3).
%
%   f(x) = (1+x)^3, [a,b] = [0,3], lambda = 0.5, mu = 2, h(t) = 1 (P-convex).
%   Requires LHS_exact.m, Kfun.m, RHS_thm1.m on the MATLAB path (same folder).

clear; close all; clc;

f      = @(x) (1 + x).^3;
fprime = @(x) abs(3*(1 + x).^2);
h      = @(t) ones(size(t));  % P-convex: h(t) = 1

a = 0; b = 3; lambda = 0.5; mu = 2;

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
title('Sensitivity: f(x)=(1+x)^3, h(t)=1 (P-convex)');
legend('Location', 'best'); grid on;
