function L = LHS_exact(f, a, b, lambda, mu, beta)
% LHS_EXACT  Exact left-hand side of the Milne-Mercer identity (Lemma 3.1).
%
%   L = LHS_exact(f, a, b, lambda, mu, beta)
%
%   f      - function handle, f(x)
%   a,b    - interval endpoints, 0 <= a < b
%   lambda,mu - reflection points, a <= lambda < mu <= b
%   beta   - fractional order, 0 < beta < 1
%
%   Computes
%     (1/3)[2f(E) - f(mid) + 2f(F)]
%       - ((1-beta)/2^beta)*(mu-lambda)^(beta-1) * [J1 + J2]
%   where E = a+b-lambda, F = a+b-mu, mid = a+b-(lambda+mu)/2, and
%     J1 = int_{mid}^{E} (E-t)^{-beta} f(t) dt   ( N3 J^beta_{mid+}{f(E)} )
%     J2 = int_{F}^{mid} (t-F)^{-beta} f(t) dt   ( N3 J^beta_{mid-}{f(F)} )

    E   = a + b - lambda;
    F   = a + b - mu;
    mid = a + b - (lambda + mu)/2;

    term1 = (2*f(E) - f(mid) + 2*f(F)) / 3;

    J1 = integral(@(t) (E - t).^(-beta) .* f(t), mid, E, ...
                  'RelTol', 1e-12, 'AbsTol', 1e-14);
    J2 = integral(@(t) (t - F).^(-beta) .* f(t), F, mid, ...
                  'RelTol', 1e-12, 'AbsTol', 1e-14);

    term2 = ((1 - beta) / 2^beta) * (mu - lambda)^(beta - 1) * (J1 + J2);

    L = term1 - term2;
end
