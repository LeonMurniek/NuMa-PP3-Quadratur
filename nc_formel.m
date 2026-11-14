function [x, w] = nc_formel(N)
% NC_FORMEL  Abgeschlossene Newton-Cotes-Formel mit N >= 2 Punkten auf [-1,1].
%
% Knoten: N äquidistante Punkte inkl. Randpunkte. Gewichte w_i = int l_i(x) dx
% mit den Lagrange-Polynomen l_i (Satz 3.5). l_i wird baryzentrisch ausgewertet,
% das Integral (Grad N-1) macht eine Gauß-Formel mit N Punkten exakt.

x = linspace(-1, 1, N)';

% Baryzentrische Gewichte für äquidistante Knoten: lambda_i = (-1)^i * binom(N-1, i)
% (Faktor h^(N-1)/(N-1)! hebt sich im Quotienten weg). Skalierung gegen Überlauf.
i      = (0:N-1)';
lambda = (-1).^i .* arrayfun(@(j) nchoosek(N-1, j), i);
lambda = lambda / max(abs(lambda));

[xg, wg] = gauss_formel(N);

w = zeros(N, 1);
for j = 1:N
    % l_j an den Gauß-Knoten:  l_j(t) = (lambda_j/(t-x_j)) / sum_m lambda_m/(t-x_m)
    lj = zeros(N, 1);
    for q = 1:N
        d = xg(q) - x;
        treffer = find(d == 0, 1);
        if ~isempty(treffer)          % Gauß-Knoten fällt auf Stützstelle (z. B. t = 0 bei ungeradem N)
            lj(q) = (treffer == j);
        else
            t = lambda ./ d;
            lj(q) = t(j) / sum(t);
        end
    end
    w(j) = wg' * lj;
end
w = (w + flipud(w)) / 2;   % Symmetrie
end
