%% Programmierprojekt 3 – Aufgabe 3: Summierte Quadratur
% Budget: höchstens 48 Funktionsauswertungen auf [0,1].
% N Punkte pro Teilintervall -> M Teilintervalle:
%   Gauß (keine Randpunkte):          N*M          <= 48  ->  M = floor(48/N)
%   Newton-Cotes (Randpunkte geteilt): M*(N-1) + 1 <= 48  ->  M = floor(47/(N-1))

clear; close all; clc;

a = 0; b = 1; budget = 48;

f = {@(x) exp(x), ...
     @(x) 1 ./ (1 + 400*(x - 0.3).^2), ...
     @(x) abs(sin(5*x))};
namen = {'f_1 = e^x', 'f_2 = 1/(1+400(x-0.3)^2)', 'f_3 = |sin(5x)|'};

% Exakte Integrale (analytisch)
I = [exp(1) - 1, ...
     (atan(20*(1 - 0.3)) + atan(20*0.3)) / 20, ...
     (3 + cos(5)) / 5];

Ns = 1:budget;
errG = nan(budget, 3);  errN = nan(budget, 3);
MG = nan(budget, 1);    MN = nan(budget, 1);
nG = nan(budget, 1);    nN = nan(budget, 1);   % tatsächlich verwendete Auswertungen

for N = Ns
    % --- Gauß ---
    M = floor(budget / N);  MG(N) = M;
    [xh, wh] = gauss_formel(N);
    [x, w]   = summierte_formel(xh, wh, a, b, M);
    nG(N) = numel(x);
    for j = 1:3, errG(N, j) = abs(quadratur(f{j}(x), w) - I(j)); end

    % --- Newton-Cotes (ab N = 2) ---
    if N >= 2
        M = floor((budget - 1) / (N - 1));  MN(N) = M;
        [xh, wh] = nc_formel(N);
        [x, w]   = summierte_formel(xh, wh, a, b, M);
        nN(N) = numel(x);
        for j = 1:3, errN(N, j) = abs(quadratur(f{j}(x), w) - I(j)); end
    end
end
assert(all(nG <= budget) && all(nN(2:end) <= budget));

boden = 1e-17;    % Fehler 0 lässt sich nicht logarithmisch plotten
figure('Name', 'Aufgabe 3', 'Position', [100 100 1300 420]);
for j = 1:3
    subplot(1, 3, j);
    semilogy(Ns, max(errG(:, j), boden), 'o-', 'LineWidth', 1.2); hold on;
    eNp = errN(:, j); eNp(~isnan(eNp)) = max(eNp(~isnan(eNp)), boden);
    semilogy(Ns, eNp, 's-', 'LineWidth', 1.2);
    grid on; xlabel('N (Punkte pro Teilintervall)'); ylabel('|Q - I|');
    title(namen{j}); ylim([1e-17 1]);
    if j == 1, legend('Gauß', 'Newton-Cotes', 'Location', 'northeast'); end
end
speichere_plot('a3_summierte_quadratur.png');

% Beste Aufteilung je Funktion
fprintf('Beste Aufteilung (Budget %d):\n', budget);
for j = 1:3
    [eg, ig] = min(errG(:, j));  [en, in] = min(errN(:, j));
    fprintf('  %-26s Gauß: N=%2d, M=%2d, Fehler %.2e | NC: N=%2d, M=%2d, Fehler %.2e\n', ...
        namen{j}, ig, MG(ig), eg, in, MN(in), en);
end
fprintf('\nFehler für ausgewählte N (Gauß / NC):\n  N  M_G M_NC |');
fprintf('      f1 G      f1 NC      f2 G      f2 NC      f3 G      f3 NC\n');
for N = [2 3 4 5 6 8 12 16 24 48]
    fprintf('%3d %4d %4d | %9.1e %9.1e %9.1e %9.1e %9.1e %9.1e\n', N, MG(N), MN(N), ...
        errG(N,1), errN(N,1), errG(N,2), errN(N,2), errG(N,3), errN(N,3));
end
