function y = Kfun(t, beta)
% KFUN  Kernel K(t) = 4/3 - (1-t)^(1-beta) from Lemma 3.1.
    y = 4/3 - (1 - t).^(1 - beta);
end
