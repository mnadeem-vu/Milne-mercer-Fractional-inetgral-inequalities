% FIG_THM1_VS_THM2_VS_THM3  Reproduces thm1_vs_thm2_vs_thm3.png.
%
%   f(x) = e^x on [0,1], lambda = 1/3, mu = 2/3, h(t) = t (convex case).
%   Compares Theorem 3.5 (triangle), Theorem 3.14 (Hoelder, p=2), and
%   Theorem 3.22 (power-mean, p=2 and p=5), against the exact LHS.
%   Requires LHS_exact.m, Kfun.m, RHS_thm1.m, RHS_thm2.m, RHS_thm3.m
%   (same folder).

clear; close all; clc;

f      = @(x) exp(x);
fprime = @(x) abs(exp(x));
h      = @(t) t;  % Convex case

a = 0; b = 1; lambda = 1/3; mu = 2/3;

beta_grid = linspace(0.05, 0.95, 40);

LHS = zeros(size(beta_grid));
T1  = zeros(size(beta_grid));
T2  = zeros(size(beta_grid));   % Hoelder, p = 2
T3a = zeros(size(beta_grid));   % power-mean, p = 2
T3b = zeros(size(beta_grid));   % power-mean, p = 5

for k = 1:numel(beta_grid)
    beta = beta_grid(k);
    LHS(k) = abs(LHS_exact(f, a, b, lambda, mu, beta));
    T1(k)  = RHS_thm1(fprime, a, b, lambda, mu, beta, h);
    T2(k)  = RHS_thm2(fprime, a, b, lambda, mu, beta, h, 2);
    T3a(k) = RHS_thm3(fprime, a, b, lambda, mu, beta, h, 2);
    T3b(k) = RHS_thm3(fprime, a, b, lambda, mu, beta, h, 5);
end

figure; hold on;
plot(beta_grid, LHS, 'k--', 'LineWidth', 1.6, 'DisplayName', 'LHS (exact)');
plot(beta_grid, T1,  'o-', 'DisplayName', 'Theorem 1 (triangle)');
plot(beta_grid, T2,  's-', 'DisplayName', 'Theorem 2 (Hoelder, p=2)');
plot(beta_grid, T3a, '^-', 'DisplayName', 'Theorem 3 (power-mean, p=2)');
plot(beta_grid, T3b, 'v-', 'DisplayName', 'Theorem 3 (power-mean, p=5)');
xlabel('\beta'); ylabel('Bound value');
title('Theorem 1 vs 2 vs 3: f(x)=e^x, h(t)=t, [0,1]');
legend('Location', 'best'); grid on;
