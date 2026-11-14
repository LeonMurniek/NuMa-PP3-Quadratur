%% Programmierprojekt 3 – Aufgabe 1 und 2: Exaktheitsgrad und Kondition
% Numerische Mathematik für Elektrotechnik, WS 2026/27
%
%  Gauß-Legendre (N = 1..40) und Newton-Cotes (N = 2..40):
%  gemessener Exaktheitsgrad und Konditionszahl sum_i |w_i| über N.

clear; close all; clc;

Nmax = 40;
tols = [1e-10 1e-14];         % Abweichung, ab der ein Monom "nicht mehr exakt" integriert wird

% Formeln einmal vorab berechnen
G = cell(Nmax, 1);  NC = cell(Nmax, 1);
for N = 1:Nmax
    [G{N}(:,1),  G{N}(:,2)]  = gauss_formel(N);
    if N >= 2, [NC{N}(:,1), NC{N}(:,2)] = nc_formel(N); end
end
save('formeln.mat', 'G', 'NC');

mG = nan(Nmax, 2);  mN = nan(Nmax, 2);
kG = nan(Nmax, 1);  kN = nan(Nmax, 1);
for N = 1:Nmax
    kG(N) = sum(abs(G{N}(:,2)));
    if N >= 2, kN(N) = sum(abs(NC{N}(:,2))); end
    for t = 1:2
        mG(N,t) = exaktheitsgrad(G{N}(:,1), G{N}(:,2), tols(t));
        if N >= 2, mN(N,t) = exaktheitsgrad(NC{N}(:,1), NC{N}(:,2), tols(t)); end
    end
end

% Theorie: Gauß 2N-1; Newton-Cotes N (N ungerade) bzw. N-1 (N gerade)
Nv   = (1:Nmax)';
thG  = 2*Nv - 1;
thN  = Nv - mod(Nv + 1, 2);       % N ungerade -> N, N gerade -> N-1

% Kontrolle: Summe der Gewichte = 2 (Exaktheit für f = 1)
fprintf('max |sum(w) - 2|:  Gauß %.1e,  NC %.1e\n', ...
    max(cellfun(@(g) abs(sum(g(:,2)) - 2), G)), max(cellfun(@(g) abs(sum(g(:,2)) - 2), NC(2:end))));

fprintf('\nGemessener Grad bei tol = %.0e / %.0e\n', tols);
fprintf('   N | Gauß (Theorie)    | NC (Theorie)     | sum|w| Gauß | sum|w| NC\n');
for N = [1 2 3 4 5 6 8 10 12 15 20 25 30 35 40]
    fprintf('%4d | %3d /%3d (%3d)   | %3d /%3d (%3d)  | %9.3f   | %9.3e\n', N, mG(N,:), thG(N), ...
        mN(N,:), thN(N), kG(N), kN(N));
end

% Welche Genauigkeit ist durch Rundung bei NC überhaupt noch möglich?
% Datenfehler-Verstärkung ~ eps * sum|w|
fprintf('\neps * sum|w|:  Gauß N=40: %.1e,  NC N=40: %.1e,  NC N=20: %.1e\n', ...
    eps*kG(40), eps*kN(40), eps*kN(20));

figure('Name', 'Aufgabe 1/2: Exaktheitsgrad', 'Position', [100 100 1000 420]);
subplot(1, 2, 1);
plot(Nv, mG(:,1), 'o-', Nv, mG(:,2), 'o--', 'LineWidth', 1.2, 'Color', [0 0.447 0.741]); hold on;
plot(Nv, mN(:,1), 's-', Nv, mN(:,2), 's--', 'LineWidth', 1.2, 'Color', [0.85 0.325 0.098]);
plot(Nv, thG, 'k-', Nv, thN, 'k:', 'LineWidth', 1);
grid on; xlabel('N (Punkte)'); ylabel('Exaktheitsgrad m');
legend('Gauß tol 1e-10', 'Gauß tol 1e-14', 'NC tol 1e-10', 'NC tol 1e-14', ...
    'Theorie Gauß 2N-1', 'Theorie NC N bzw. N-1', 'Location', 'northwest');
title('Gemessener Exaktheitsgrad');

subplot(1, 2, 2);
semilogy(Nv, kG, 'o-', 'LineWidth', 1.3); hold on;
semilogy(Nv, kN, 's-', 'LineWidth', 1.3);
yline(2, 'k--');
grid on; xlabel('N (Punkte)'); ylabel('\Sigma_i |w_i|');
legend('Gauß', 'Newton-Cotes', 'Minimum 2 (alle w_i > 0)', 'Location', 'northwest');
title('Konditionszahl der Quadraturformel');
speichere_plot('a1_exaktheitsgrad_kondition.png');
