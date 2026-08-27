function R = RHS_thm2(fprime, a, b, lambda, mu, beta, h, p)
% RHS_THM2  Bound of Theorem 3.14 (Hoelder-type Milne-Mercer bound).
%
%   R = RHS_thm2(fprime, a, b, lambda, mu, beta, h, p)
%
%   fprime - function handle for |f'|
%   h      - function handle for h(t)
%   p      - Hoelder exponent, p > 1;  q = p/(p-1)

    q = p/(p-1);
    E = a + b - lambda;
    F = a + b - mu;
    h_half = h(0.5);

    Kq  = integral(@(t) Kfun(t,beta).^q, 0, 1, 'RelTol',1e-12,'AbsTol',1e-14)^(1/q);
    Hp  = integral(@(t) h(t), 0, 1, 'RelTol',1e-12,'AbsTol',1e-14)^(1/p);

    fE = fprime(E); fF = fprime(F);
    term1 = ((1+h_half)*fE^p + h_half*fF^p)^(1/p);
    term2 = ((1+h_half)*fF^p + h_half*fE^p)^(1/p);

    R = (mu - lambda)/4 * Kq * Hp * (term1 + term2);
end
