% FIG_MATRIX_VS_LAMBDA  Reproduces matrix_vs_lambda.png
% (Section 5.5, Matrix Inequality application).
%
%   Commuting positive-definite A = diag(a1,a2), B = diag(b1,b2),
%   X = diag(x1,x2), with a1=2,b1=1,x1=1.5 and a2=0.5,b2=3,x2=0.7.
%   psi(theta) = sum_i x_i (a_i^theta b_i^(1-theta) + a_i^(1-theta) b_i^theta).
%   Compares the exact LHS of Proposition 5.23 against its bound as
%   lambda varies, with mu = 0.8 fixed. Self-contained.

clear; close all; clc;

a_vals = [2, 0.5];
b_vals = [1, 3];
x_vals = [1.5, 0.7];

phi   = @(theta,ai,bi) ai.^theta .* bi.^(1-theta) + ai.^(1-theta) .* bi.^theta;
phip  = @(theta,ai,bi) log(ai/bi) .* (ai.^theta.*bi.^(1-theta) - ai.^(1-theta).*bi.^theta);

psi   = @(theta) sum(arrayfun(@(i) x_vals(i)*phi(theta,a_vals(i),b_vals(i)), 1:2));
psip  = @(theta) sum(arrayfun(@(i) x_vals(i)*phip(theta,a_vals(i),b_vals(i)), 1:2));

mu = 0.8;
lambda_list = 0.05:0.05:0.70;

LHS = zeros(size(lambda_list));
RHS = zeros(size(lambda_list));

for k = 1:numel(lambda_list)
    lambda = lambda_list(k);
    E = 1 - lambda;
    F = 1 - mu;
    mid = 1 - (lambda+mu)/2;

    term1 = (2*psi(E) - psi(mid) + 2*psi(F)) / 3;
    integral_term = integral(@(t) arrayfun(psi, t), F, E, ...
                              'RelTol', 1e-12, 'AbsTol', 1e-14);
    LHS(k) = abs(term1 - integral_term/(mu-lambda));

    RHS(k) = 5*(mu-lambda)/24 * (abs(psip(E)) + abs(psip(F)));
end

fprintf(' lambda    Exact LHS     Bound (Prop 5.23)\n');
for k = 1:numel(lambda_list)
    fprintf('  %.2f     %.6f      %.6f\n', lambda_list(k), LHS(k), RHS(k));
end

figure;
plot(lambda_list, LHS, 'o-', 'DisplayName', 'Exact LHS'); hold on;
plot(lambda_list, RHS, 's--', 'DisplayName', 'Bound (Prop. 5.23)');
xlabel('\lambda (with \mu = 0.8 fixed)'); ylabel('Value');
title('Matrix application: exact LHS vs. bound as \lambda varies');
legend('Location', 'best'); grid on;
