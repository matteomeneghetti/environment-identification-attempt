
data = VarName3;
out = VarName2;

Ak = data(1);
b = out(1);
xk = [];

% Step iniziale (k = 0)
Pk = inv(Ak.' * Ak);
xk = Pk * (Ak.' * b);

% Step iterativo (k = 1,2,3,...,n)

% Concatenare nuovi dati b e misurazioni A


%%

M = size(Ak, 1);
for i = 2:height(data)
    Ak = data(i, :);
    b = out(i, :);
    if isnan(Ak) || isnan(b)
        continue
    end
    M = size(Ak, 1);
    % Usare SMW per calcolare Pk
    Pk = Pk - Pk * Ak.' * inv(eye(M) + Ak*Pk*Ak.') * Ak * Pk;
    if isnan(Pk)
        continue
    end
    K = Pk*Ak.';
    xk = xk + K * (b - Ak*xk);
end