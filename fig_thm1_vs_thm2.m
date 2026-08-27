% FIG_THM1_VS_THM2  Reproduces thm1_vs_thm2.png.
%
%   f(x) = e^x on [0,1], lambda = 1/3, mu = 2/3, h(t) = t (convex case).
%   Compares Theorem 3.5 (triangle) against Theorem 3.14 (Hoelder) for
%   p = 2, 3, 5, against the exact LHS.
%   Requires LHS_exact.m, Kfun.m, RHS_thm1.m, RHS_thm2.m (same folder).

clear; close all; clc;

f      = @(x) exp(x);
fprime = @(x) abs(exp(x));
h      = @(t) t;  % Convex case

a = 0; b = 1; lambda = 1/3; mu = 2/3;

beta_grid = linspace(0.05, 0.95, 40);
p_list = [2 3 5];

LHS  = zeros(size(beta_grid));
T1   = zeros(size(beta_grid));
T2   = zeros(numel(p_list), numel(beta_grid));

for k = 1:numel(beta_grid)
    beta = beta_grid(k);
    LHS(k) = abs(LHS_exact(f, a, b, lambda, mu, beta));
    T1(k)  = RHS_thm1(fprime, a, b, lambda, mu, beta, h);
    for j = 1:numel(p_list)
        T2(j,k) = RHS_thm2(fprime, a, b, lambda, mu, beta, h, p_list(j));
    end
end

figure; hold on;
plot(beta_grid, LHS, 'k--', 'LineWidth', 1.6, 'DisplayName', 'LHS (exact)');
plot(beta_grid, T1,  'o-',  'DisplayName', 'Theorem 1 (triangle)');
for j = 1:numel(p_list)
    plot(beta_grid, T2(j,:), '-', 'Marker', 's', ...
         'DisplayName', sprintf('Theorem 2, p=%d', p_list(j)));
end
xlabel('\beta'); ylabel('Bound value');
title('Theorem 1 vs Theorem 2 bounds: f(x)=e^x, h(t)=t, [0,1]');
legend('Location', 'best'); grid on;
