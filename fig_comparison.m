% FIG_COMPARISON  Reproduces comparison.png.
%
%   RHS bound of Theorem 3.5 across convexity classes, for f(x) = e^x on
%   [0,1], plotted against beta, for P-convex, s-convex (s=0.3,0.5,0.7),
%   and (ordinary) convex h.
%   Requires Kfun.m, RHS_thm1.m on the MATLAB path (same folder).

clear; close all; clc;

fprime = @(x) abs(exp(x));
a = 0; b = 1; lambda = 1/3; mu = 2/3;

beta_grid = linspace(0.05, 0.95, 40);

h_list = {@(t) ones(size(t)), ...        % P-convex
          @(t) t.^0.3, ...               % s-convex, s=0.3
          @(t) t.^0.5, ...               % s-convex, s=0.5
          @(t) t.^0.7, ...               % s-convex, s=0.7
          @(t) t};                       % Convex
labels = {'P-convex (h=1)', 's-convex, s=0.3', 's-convex, s=0.5', ...
          's-convex, s=0.7', 'Convex (h=t)'};

RHS = zeros(numel(h_list), numel(beta_grid));
for j = 1:numel(h_list)
    for k = 1:numel(beta_grid)
        RHS(j,k) = RHS_thm1(fprime, a, b, lambda, mu, beta_grid(k), h_list{j});
    end
end

figure; hold on;
for j = 1:numel(h_list)
    plot(beta_grid, RHS(j,:), 'LineWidth', 1.6, 'DisplayName', labels{j});
end
xlabel('\beta'); ylabel('RHS bound');
title('RHS bound across convexity classes: f(x)=e^x on [0,1]');
legend('Location', 'best'); grid on;
