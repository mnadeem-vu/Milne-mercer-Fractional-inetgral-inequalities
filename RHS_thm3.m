function R = RHS_thm3(fprime, a, b, lambda, mu, beta, h, p)
% RHS_THM3  Bound of Theorem 3.22 (improved power-mean Milne-Mercer bound).
%
%   R = RHS_thm3(fprime, a, b, lambda, mu, beta, h, p)
%
%   fprime - function handle for |f'|
%   h      - function handle for h(t)
%   p      - power-mean exponent, p >= 1

    E = a + b - lambda;
    F = a + b - mu;
    h_half = h(0.5);

    Omega = integral(@(t) Kfun(t,beta), 0, 1, 'RelTol',1e-12,'AbsTol',1e-14);
    R1 = integral(@(t) Kfun(t,beta).*h(t),     0, 1, 'RelTol',1e-12,'AbsTol',1e-14);
    R2 = integral(@(t) Kfun(t,beta).*h(1-t),   0, 1, 'RelTol',1e-12,'AbsTol',1e-14);

    fE = fprime(E); fF = fprime(F);
    term1 = ((R1 + h_half*R2)*fE^p + h_half*R2*fF^p)^(1/p);
    term2 = ((R1 + h_half*R2)*fF^p + h_half*R2*fE^p)^(1/p);

    R = (mu - lambda)/4 * Omega^(1 - 1/p) * (term1 + term2);
end
