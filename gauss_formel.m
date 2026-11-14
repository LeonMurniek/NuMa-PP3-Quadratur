function [x, w] = gauss_formel(N)
% GAUSS_FORMEL  Gauß-Legendre-Formel mit N Punkten auf [-1,1] (Golub-Welsch).
%
% Die Knoten sind die Eigenwerte der Jacobi-Matrix J (symmetrisch, tridiagonal,
% Diagonale 0, Nebendiagonale beta_k = 1/sqrt(4 - k^-2), k = 1..N-1), siehe Prop. 3.24.
% Die Gewichte sind w_i = 2 * (erste Komponente des normierten Eigenvektors i)^2.

if N == 1
    x = 0; w = 2; return
end

k    = 1:N-1;
beta = 1 ./ sqrt(4 - k.^-2);
J    = diag(beta, 1) + diag(beta, -1);

[V, D] = eig(J);
[x, idx] = sort(diag(D));
w = 2 * (V(1, idx).^2)';

% Symmetrie erzwingen (Knoten x_i = -x_{N+1-i}, gleiche Gewichte), entfernt Rundungsrauschen
x = (x - flipud(x)) / 2;
w = (w + flipud(w)) / 2;
end
