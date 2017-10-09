%% Este script calculara las expresiones de la ecuacion de potencia para los generadores

% dv = sym('dv', [1 n_gen]);

Gp = real(YKronPre);
Bp = imag(YKronPre);
gp = real(YShuntKronPre);
bp = imag(YShuntKronPre);

% PeEqPre = sym(zeros(n_gen, n_gen));
% for i = 1:n_gen
%     for j = 1:n_gen
%         if(i ~= j)
%             PeEqPre(i, j) = (-Gp(i,j) + gp(i,j))*Em(i)^2 + Em(i)*Em(j)*(Gp(i,j)*cos(dv(i)-dv(j)) + Bp(i,j)*sin(dv(i)-dv(j)));
%         end
%     end
% end

%% Aprovechamos y calculamos las potencias pre-falla que deben dar igual a lo que da el flujo de carga
Pe0m = zeros(n_gen, n_gen);
Pe0 = zeros(1, n_gen);
for i = 1:n_gen
    for j = 1:n_gen
        if(i ~= j)
            Pe0m(i, j) = (-Gp(i,j) + gp(i,j))*Em(i)^2 + Em(i)*Em(j)*(Gp(i,j)*cos(d0(i)-d0(j)) + Bp(i,j)*sin(d0(i)-d0(j)));
        end
    end
end

for i = 1:n_gen
    Pe0(i) = sum(Pe0m(i, 1:n_gen));
end
