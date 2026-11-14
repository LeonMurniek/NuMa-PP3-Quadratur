function [x, w] = summierte_formel(xh, wh, a, b, M)
% SUMMIERTE_FORMEL  Formel (xh, wh) auf M gleich lange Teilintervalle von [a,b] anwenden.
% Punkte, die in zwei Teilintervallen liegen (Randpunkte bei Newton-Cotes), werden
% nur einmal ausgewertet; ihre Gewichte addieren sich.

h = (b - a) / M;
X = []; W = [];
for j = 1:M
    [xj, wj] = transformiere(xh, wh, a + (j-1)*h, a + j*h);
    X = [X; xj]; W = [W; wj]; %#ok<AGROW>
end

[X, idx] = sort(X);
W = W(idx);
neu = [true; diff(X) > 1e-12 * (b - a)];   % neuer Punkt, falls nicht gleich dem Vorgänger
grp = cumsum(neu);
x   = X(neu);
w   = accumarray(grp, W);
end
