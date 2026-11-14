function m = exaktheitsgrad(xh, wh, tol)
% EXAKTHEITSGRAD  Gemessener Exaktheitsgrad einer Formel auf [-1,1].
% Integriert x^k, k = 0,1,2,... und vergleicht mit dem exakten Integral
%   int_{-1}^{1} x^k dx = 2/(k+1) (k gerade), 0 (k ungerade).
% Sobald die Abweichung tol überschreitet, ist m = k-1.

kmax = 300;
for k = 0:kmax
    exakt = (mod(k, 2) == 0) * 2 / (k + 1);
    Q     = quadratur(xh(:).^k, wh);
    if abs(Q - exakt) > tol
        m = k - 1;
        return
    end
end
m = kmax;
end
