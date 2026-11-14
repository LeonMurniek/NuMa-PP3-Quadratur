function Q = quadratur(y, w)
% QUADRATUR  Q = sum_i w_i * y_i  (y = Funktionswerte an den Knoten).
Q = sum(w(:) .* y(:));
end
