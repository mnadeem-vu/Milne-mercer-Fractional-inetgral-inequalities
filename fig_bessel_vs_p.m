% FIG_BESSEL_VS_P  Reproduces bessel_vs_p.png (Modified Bessel application).
%
%   Uses the closed-form classical (beta -> 0+) propositions 5.1-5.3
%   directly: J_nu(u) = 2^nu * Gamma(nu+1) * u^(-nu) * I_nu(u),
%   J_nu'(u) = u/(2(nu+1)) J_{nu+1}(u),
%   J_nu''(u) = (1/(4(nu+1))) [2 J_{nu+1}(u) + u^2/(nu+2) J_{nu+2}(u)].
%
%   Data: nu = 1.2, [a,b] = [1,4], lambda = 1.5, mu = 3.0, s = 0.6.
%   Self-contained: no external files required beyond MATLAB's besseli.

clear; close all; clc;

nu = 1.2;
a = 1; b = 4; lambda = 1.5; mu = 3.0; s = 0.6;

Jnu    = @(u,n) 2^n .* gamma(n+1) .* u.^(-n) .* besseli(n, u);
Jnup   = @(u,n) u ./ (2*(n+1)) .* Jnu(u, n+1);
Jnupp  = @(u,n) (2*Jnu(u,n+1) + u.^2/(n+2).*Jnu(u,n+2)) / (4*(n+1));

E = a + b - lambda;
F = a + b - mu;
mid = a + b - (lambda+mu)/2;

% Exact LHS (via the closed-form antiderivative J_nu, Proposition 5.1's LHS)
% Using f(u) = J_nu'(u) directly: (1/3)[2f(E)-f(mid)+2f(F)] - [J_nu(E)-J_nu(F)]/(mu-lambda)
LHS = abs( (2*Jnup(E,nu) - Jnup(mid,nu) + 2*Jnup(F,nu))/3 ...
           - (Jnu(E,nu) - Jnu(F,nu))/(mu-lambda) );

% Triangle bound (Proposition 5.1)
Cs = @(s) 4*(1+2^(1-s))/(3*(s+1)) - 1/((s+1)*(s+2)) - 2^(1-s)/(s+2);
% Bracket in Prop. 5.1 is |2 J_{nu+1}(u) + u^2/(nu+2) J_{nu+2}(u)|:
absterm = @(u) abs(2*Jnu(u,nu+1) + u.^2/(nu+2).*Jnu(u,nu+2));
RHS_triangle = (mu-lambda)/(16*(nu+1)) * (absterm(E) + absterm(F)) * Cs(s);

p_list = [2 3 5 10];
RHS_holder = zeros(size(p_list));
RHS_pm     = zeros(size(p_list));

Jpp_E = Jnupp(E, nu);
Jpp_F = Jnupp(F, nu);

for i = 1:numel(p_list)
    p = p_list(i);
    q = p/(p-1);

    % Hoelder bound (Proposition 5.2)
    part1 = ((4/3)^(q+1) - (1/3)^(q+1)) / (q+1);
    term_int = part1^(1/q);
    term_s   = (1/(s+1))^(1/p);
    sig1 = ((1+2^(-s))*abs(Jpp_E)^p + 2^(-s)*abs(Jpp_F)^p)^(1/p);
    sig2 = ((1+2^(-s))*abs(Jpp_F)^p + 2^(-s)*abs(Jpp_E)^p)^(1/p);
    RHS_holder(i) = (mu-lambda)/4 * term_int * term_s * (sig1 + sig2);

    % Power-mean bound (Proposition 5.3); classical limits Omega0=5/6,
    % R1_0, R2_0 as beta -> 0+ (closed forms from Corollary "classical-pm-s")
    Omega0 = 5/6;
    R1_0 = 4/(3*(s+1)) - 1/((s+1)*(s+2));
    R2_0 = 4/(3*(s+1)) - 1/(s+2);
    t1 = ((R1_0 + 2^(-s)*R2_0)*abs(Jpp_E)^p + 2^(-s)*R2_0*abs(Jpp_F)^p)^(1/p);
    t2 = ((R1_0 + 2^(-s)*R2_0)*abs(Jpp_F)^p + 2^(-s)*R2_0*abs(Jpp_E)^p)^(1/p);
    RHS_pm(i) = (mu-lambda)/4 * Omega0^(1-1/p) * (t1 + t2);
end

fprintf('LHS = %.6f, Triangle RHS = %.6f\n', LHS, RHS_triangle);
fprintf('  p     Hoelder      power-mean\n');
for i = 1:numel(p_list)
    fprintf('  %2d   %.6f    %.6f\n', p_list(i), RHS_holder(i), RHS_pm(i));
end

figure; hold on;
plot(p_list, RHS_triangle*ones(size(p_list)), 'o-', 'DisplayName', 'Prop. 5.1 (triangle)');
plot(p_list, RHS_holder, 's-', 'DisplayName', 'Prop. 5.3 (Hoelder)');
plot(p_list, RHS_pm,     '^-', 'DisplayName', 'Prop. 5.5 (power-mean)');
yline(LHS, 'k--', 'DisplayName', 'LHS (exact)');
xlabel('p'); ylabel('Bound value');
title('Modified Bessel application: RHS bounds vs p');
legend('Location', 'best'); grid on;
