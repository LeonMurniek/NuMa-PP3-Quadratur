function [x, w] = transformiere(xh, wh, a, b)
% TRANSFORMIERE  Formel (xh, wh) von [-1,1] auf [a,b] übertragen (Lemma 3.3).
x = (b - a) / 2 * xh + (a + b) / 2;
w = (b - a) / 2 * wh;
end
