%% Programmierprojekt 3 – Aufgabe 4: Energie des Kondensators aus Messwerten
% C(v) = 10 uF / (1 + (v/8V)^2),   W = int_0^12 v*C(v) dv   (in uJ)
% Messgerät: C(v_i) gerundet auf 0,01 uF.  Ziel: |Q - W| < 1 uJ mit möglichst wenigen Messpunkten.
%
% Fehlerzerlegung:  Q(f~) - W = [Q(f) - W]  +  [Q(f~) - Q(f)]
%                              Quadraturfehler   Datenfehler

clear; close all; clc;

a = 0; b = 12;
C   = @(v) 10 ./ (1 + (v/8).^2);        % uF
f   = @(v) v .* C(v);                   % exakte Leistung pro dv   [V*uF = uJ/V]
ft  = @(v) v .* (round(C(v), 2));       % was das Messgerät liefert (auf 0,01 uF gerundet)
W   = 320 * log(1 + 144/64);            % exakt: 10*32*ln(1 + 12^2/8^2) = 377.1703... uJ
dC  = 0.005;                            % maximaler Rundungsfehler von C (uF)
fprintf('W exakt = %.6f uJ  (Kontrolle mit integral: %.6f)\n\n', W, integral(f, a, b, 'AbsTol', 1e-12));

% Verfahren: Name, Funktion für Formel (xh, wh) auf [-1,1], N pro Teilintervall
verfahren = {
    'Gauß einfach',        @(N) gauss_formel(N),  'N',  []
    'Newton-Cotes einfach',@(N) nc_formel(N),     'N',  []
    'Trapez summiert',     @(~) nc_formel(2),     'M',  2
    'Simpson summiert',    @(~) nc_formel(3),     'M',  3
    'Boole summiert (NC5)',@(~) nc_formel(5),     'M',  5
    'Gauß summiert (N=3)', @(~) gauss_formel(3),  'M',  3
};

tab = [];   % Spalten: Verfahren, Punkte, e_Q, e_D, Schranke
for v = 1:size(verfahren, 1)
    if verfahren{v, 3} == 'N', param = 2:40;  else, param = 1:60; end
    for p = param
        if verfahren{v, 3} == 'N'
            [xh, wh] = verfahren{v, 2}(p);  M = 1;
        else
            [xh, wh] = verfahren{v, 2}(p);  M = p;
        end
        [x, w] = summierte_formel(xh, wh, a, b, M);
        Qf  = quadratur(f(x),  w);
        Qft = quadratur(ft(x), w);
        eQ  = Qf - W;
        eD  = Qft - Qf;
        schr = dC * sum(abs(w) .* abs(x));   % Worst-Case-Datenfehler: |f~-f| <= v*0.005 uF
        tab = [tab; v, numel(x), eQ, eD, schr]; %#ok<AGROW>
    end
end

% Minimale Punktzahl je Verfahren: erstes P, ab dem Ziel erreicht ist (und für größere P dieser Familie bleibt)
fprintf('%-22s | Punkte (Ziel real) | Punkte (Ziel garantiert) |  e_Q       e_D     Schranke_D\n', 'Verfahren');
for v = 1:size(verfahren, 1)
    T = tab(tab(:,1) == v, :);
    real_ok = abs(T(:,3) + T(:,4)) < 1;
    gar_ok  = abs(T(:,3)) + T(:,5) < 1;
    ir = find(real_ok, 1);  ig = find(gar_ok, 1);
    s1 = '-'; s2 = '-'; zeile = '';
    if ~isempty(ir), s1 = sprintf('%d', T(ir,2)); end
    if ~isempty(ig), s2 = sprintf('%d', T(ig,2)); zeile = sprintf('%9.2e %9.2e %9.2e', T(ig,3), T(ig,4), T(ig,5)); end
    fprintf('%-22s | %18s | %24s | %s\n', verfahren{v,1}, s1, s2, zeile);
end

% Detailtabelle Gauß einfach
fprintf('\nGauß einfach: N | e_Q (uJ) | e_D (uJ) | Worst-Case-Schranke | gesamt\n');
T = tab(tab(:,1) == 1, :);
for r = 1:12
    fprintf('  %2d | %10.3e | %10.3e | %10.3e | %10.3e\n', T(r,2), T(r,3), T(r,4), T(r,5), T(r,3)+T(r,4));
end
fprintf('\nNewton-Cotes einfach: N | e_Q | e_D | Schranke\n');
T = tab(tab(:,1) == 2, :);
for r = [3 5 7 9 11 13 15 19 25 29 39]
    fprintf('  %2d | %10.3e | %10.3e | %10.3e\n', T(r,2), T(r,3), T(r,4), T(r,5));
end

% ---- Plots ----
farben = lines(6);
figure('Name', 'Aufgabe 4', 'Position', [100 100 1300 450]);
subplot(1, 2, 1);
for v = 1:size(verfahren, 1)
    T = tab(tab(:,1) == v, :);
    loglog(T(:,2), max(abs(T(:,3) + T(:,4)), 1e-9), '.-', 'Color', farben(v,:), 'LineWidth', 1.1, 'MarkerSize', 8); hold on;
end
yline(1, 'k--', '1 \muJ', 'LineWidth', 1.2);
grid on; xlabel('Anzahl Messpunkte'); ylabel('|Q(f~) - W|  in \muJ');
legend(verfahren(:,1), 'Location', 'southwest'); xlim([1 120]);
title('Gesamtfehler mit gerundeten Messwerten');

subplot(1, 2, 2);
for v = [1 2]
    T = tab(tab(:,1) == v, :);
    loglog(T(:,2), max(abs(T(:,3)), 1e-12), '-', 'Color', farben(v,:), 'LineWidth', 1.3); hold on;
    loglog(T(:,2), max(abs(T(:,4)), 1e-12), '--', 'Color', farben(v,:), 'LineWidth', 1.3);
    loglog(T(:,2), T(:,5), ':', 'Color', farben(v,:), 'LineWidth', 1.3);
end
yline(1, 'k--', '1 \muJ');
grid on; xlabel('Anzahl Messpunkte N'); ylabel('Fehler in \muJ'); ylim([1e-12 1e4]);
legend('Gauß: Quadraturfehler', 'Gauß: Datenfehler', 'Gauß: Schranke Datenfehler', ...
       'NC: Quadraturfehler', 'NC: Datenfehler', 'NC: Schranke Datenfehler', 'Location', 'southwest');
title('Zerlegung in Quadratur- und Datenfehler');
speichere_plot('a4_kondensator.png');
