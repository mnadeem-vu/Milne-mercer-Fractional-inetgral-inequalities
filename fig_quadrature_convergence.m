% FIG_QUADRATURE_CONVERGENCE  Reproduces quadrature_convergence.png
% (Section 5.4, Milne's Quadrature Formula application).
%
%   f(x) = e^x on [a,b] = [0,1]. Composite Milne rule on a uniform
%   partition into n subintervals; compares the actual quadrature error
%   against the bound of Proposition 5.20 (eq. 11), C(1) = 5/6 (s=1,
%   ordinary convexity). Self-contained.

clear; close all; clc;

f      = @(x) exp(x);
fprime = @(x) abs(exp(x));
a = 0; b = 1;

I_exact = integral(f, a, b, 'RelTol', 1e-13);

n_list = [1 2 4 8 16 32 64];
E_actual = zeros(size(n_list));
Bound    = zeros(size(n_list));

for idx = 1:numel(n_list)
    n = n_list(idx);
    w = linspace(a, b, n+1);
    M = 0; bnd = 0;
    for i = 1:n
        wi = w(i); wi1 = w(i+1);
        M   = M + (wi1-wi)/3 * (2*f(wi) - f((wi+wi1)/2) + 2*f(wi1));
        bnd = bnd + (wi1-wi)^2 * (fprime(wi) + fprime(wi1));
    end
    bnd = (5/24) * bnd;   % C(1) = 5/6, bound = C(s)/4 * sum(...) => 5/24 * sum(...)
    E_actual(idx) = abs(M - I_exact);
    Bound(idx)    = bnd;
end

fprintf('   n     |E_M(f,d)|      Bound (11)\n');
for idx = 1:numel(n_list)
    fprintf('  %3d   %.6f      %.6f\n', n_list(idx), E_actual(idx), Bound(idx));
end

figure;
loglog(n_list, E_actual, 'o-', 'DisplayName', 'Actual error |E_M(f,d)|'); hold on;
loglog(n_list, Bound,    's--', 'DisplayName', 'Bound (11)');
xlabel('Number of subintervals n'); ylabel('Value');
title('Composite Milne quadrature: error vs. bound, f(x) = e^x on [0,1]');
legend('Location', 'best'); grid on;
