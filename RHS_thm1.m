function R = RHS_thm1(fprime, a, b, lambda, mu, beta, h)
% RHS_THM1  Bound of Theorem 3.5 (plain triangle-inequality Milne-Mercer bound).
%
%   R = RHS_thm1(fprime, a, b, lambda, mu, beta, h)
%
%   fprime - function handle for |f'|  (pass e.g. @(x) abs(f_prime(x)))
%   h      - function handle for h(t), the convexity-generating function
%
%   R = (mu-lambda)/4 * [|f'(E)|+|f'(F)|] * int_0^1 K(t)[h(t)+2 h(1/2) h(1-t)] dt

    E = a + b - lambda;
    F = a + b - mu;

    h_half = h(0.5);
    integrand = @(t) Kfun(t, beta) .* (h(t) + 2*h_half.*h(1 - t));
    I = integral(integrand, 0, 1, 'RelTol', 1e-12, 'AbsTol', 1e-14);

    R = (mu - lambda)/4 * (fprime(E) + fprime(F)) * I;
end
