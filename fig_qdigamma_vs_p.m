% FIG_QDIGAMMA_VS_P  Reproduces qdigamma_vs_p.png (q-digamma application).
%
%   Data: sigma = 0.6, [a,b] = [1,4], lambda = 1.5, mu = 3.0, s = 0.5.
%   psi_sigma and psi_sigma' are evaluated by truncating their defining
%   series at N = 300 terms (well past geometric convergence for
%   sigma = 0.6). Self-contained.

clear; close all; clc;

sigma = 0.6;
a = 1; b = 4; lambda = 1.5; mu = 3.0; s = 0.5;
N = 300;   % series truncation

k = (0:N)';
psi  = @(x) -log(1-sigma) + log(sigma) * sum(sigma.^(k+x) ./ (1 - sigma.^(k+x)));
psip = @(x) (log(sigma))^2 * sum(sigma.^(k+x) ./ (1 - sigma.^(k+x)).^2);

E = a + b - lambda;
F = a + b - mu;
mid = a + b - (lambda + mu)/2;

term1 = (2*psi(E) - psi(mid) + 2*psi(F)) / 3;
integral_term = integral(@(x) arrayfun(psi, x), F, E, 'RelTol', 1e-12, 'AbsTol', 1e-14);
LHS = abs(term1 - integral_term/(mu - lambda));

Cs = @(s) 4*(1+2^(1-s))/(3*(s+1)) - 1/((s+1)*(s+2)) - 2^(1-s)/(s+2);
RHS_A = (mu-lambda)/4 * (psip(E) + psip(F)) * Cs(s);

p_list = [2 3 5 10];
RHS_B = zeros(size(p_list));
RHS_C = zeros(size(p_list));

Omega0 = 5/6;
R1_0 = 4/(3*(s+1)) - 1/((s+1)*(s+2));
R2_0 = 4/(3*(s+1)) - 1/(s+2);

pE = psip(E); pF = psip(F);

for i = 1:numel(p_list)
    p = p_list(i);
    q = p/(p-1);

    part1 = ((4/3)^(q+1) - (1/3)^(q+1)) / (q+1);
    t1 = ((1+2^(-s))*pE^p + 2^(-s)*pF^p)^(1/p);
    t2 = ((1+2^(-s))*pF^p + 2^(-s)*pE^p)^(1/p);
    RHS_B(i) = (mu-lambda)/4 * part1^(1/q) * (1/(s+1))^(1/p) * (t1+t2);

    u1 = ((R1_0+2^(-s)*R2_0)*pE^p + 2^(-s)*R2_0*pF^p)^(1/p);
    u2 = ((R1_0+2^(-s)*R2_0)*pF^p + 2^(-s)*R2_0*pE^p)^(1/p);
    RHS_C(i) = (mu-lambda)/4 * Omega0^(1-1/p) * (u1+u2);
end

fprintf('LHS = %.6f, Triangle RHS_A = %.6f\n', LHS, RHS_A);
fprintf('  p     Hoelder(B)   power-mean(C)\n');
for i = 1:numel(p_list)
    fprintf('  %2d   %.6f    %.6f\n', p_list(i), RHS_B(i), RHS_C(i));
end

figure; hold on;
plot(p_list, RHS_A*ones(size(p_list)), 'o-', 'DisplayName', 'Prop. 5.14 (triangle)');
plot(p_list, RHS_B, 's-', 'DisplayName', 'Prop. 5.15 (Hoelder)');
plot(p_list, RHS_C, '^-', 'DisplayName', 'Prop. 5.16 (power-mean)');
yline(LHS, 'k--', 'DisplayName', 'LHS (exact)');
xlabel('p'); ylabel('Bound value');
title('q-digamma application: RHS bounds vs p');
legend('Location', 'best'); grid on;
