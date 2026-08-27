% FIG_MEANS_VS_P  Reproduces means_vs_p.png (Special means application).
%
%   Data: [a,b] = [1,6], lambda = 2.2, mu = 4.4, s = 0.5.
%   Uses the closed-form classical (beta -> 0+) Propositions 5.7-5.11
%   directly (arithmetic mean A, generalized log-mean Ls). Self-contained.

clear; close all; clc;

a = 1; b = 6; lambda = 2.2; mu = 4.4; s = 0.5;
E = a + b - lambda;
F = a + b - mu;

A  = @(x,y) (x + y)/2;
Ls = @(x,y,k) (y^(k+1) - x^(k+1)) / ((k+1)*(y - x));   % generalized log-mean L_k^k(x,y)
Cs = @(s) 4*(1+2^(1-s))/(3*(s+1)) - 1/((s+1)*(s+2)) - 2^(1-s)/(s+2);

% Proposition 5.7 (triangle), f(u) = u^s
LHS1 = (4*A(E^s,F^s) - A(E,F)^s)/3 - Ls(F,E,s);
RHS_triangle = (mu-lambda)*s/2 * A(E^(s-1), F^(s-1)) * Cs(s);

% Proposition 5.9 (Hoelder), f(u) = u^(-s)
LHS2 = (4*A(E^(-s),F^(-s)) - A(E,F)^(-s))/3 - Ls(F,E,-s);

p_list = [2 3 5 10];
RHS_holder = zeros(size(p_list));
RHS_pm     = zeros(size(p_list));

Omega0 = 5/6;
R1_0 = 4/(3*(s+1)) - 1/((s+1)*(s+2));
R2_0 = 4/(3*(s+1)) - 1/(s+2);

for i = 1:numel(p_list)
    p = p_list(i);
    q = p/(p-1);

    part1 = ((4/3)^(q+1) - (1/3)^(q+1)) / (q+1);
    Sigma = ((1+2^(-s))*E^(-p*(s+1)) + 2^(-s)*F^(-p*(s+1)))^(1/p) ...
          + ((1+2^(-s))*F^(-p*(s+1)) + 2^(-s)*E^(-p*(s+1)))^(1/p);
    RHS_holder(i) = (mu-lambda)*s/4 * part1^(1/q) * (1/(s+1))^(1/p) * Sigma;

    Sig3 = ((R1_0+2^(-s)*R2_0)*E^(p*(s-1)) + 2^(-s)*R2_0*F^(p*(s-1)))^(1/p) ...
         + ((R1_0+2^(-s)*R2_0)*F^(p*(s-1)) + 2^(-s)*R2_0*E^(p*(s-1)))^(1/p);
    RHS_pm(i) = (mu-lambda)*s/4 * Omega0^(1-1/p) * Sig3;
end

fprintf('LHS1 (Prop 5.7/5.11) = %.6f, LHS2 (Prop 5.9) = %.6f, Triangle RHS = %.6f\n', ...
        LHS1, LHS2, RHS_triangle);
fprintf('  p     Hoelder      power-mean\n');
for i = 1:numel(p_list)
    fprintf('  %2d   %.6f    %.6f\n', p_list(i), RHS_holder(i), RHS_pm(i));
end

figure; hold on;
plot(p_list, RHS_triangle*ones(size(p_list)), 'o-', 'DisplayName', 'Prop. 5.7 (triangle)');
plot(p_list, RHS_holder, 's-', 'DisplayName', 'Prop. 5.9 (Hoelder)');
plot(p_list, RHS_pm,     '^-', 'DisplayName', 'Prop. 5.11 (power-mean)');
yline(LHS1, 'k--', 'DisplayName', 'LHS, Props. 5.7/5.11 (f=u^s)');
yline(LHS2, 'Color', [0.5 0.5 0.5], 'LineStyle', ':', ...
      'DisplayName', 'LHS, Prop. 5.9 (f=u^{-s})');
xlabel('p'); ylabel('Bound value');
title('Special means application: RHS bounds vs p');
legend('Location', 'best'); grid on;
