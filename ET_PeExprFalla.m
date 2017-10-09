%% Este script calculara las expresiones de la ecuacion de potencia para los generadores

% dv = sym('dv', [1 n_gen]);

Gf = real(YKronFalla);
Bf = imag(YKronFalla);
gf = real(YShuntKronFalla);
bf = imag(YShuntKronFalla);

% PeEqFalla = sym(zeros(n_gen, n_gen));
% for i = 1:n_gen
%     for j = 1:n_gen
%         if(i ~= j)
%             PeEqFalla(i, j) = (-Gf(i,j) + gf(i,j))*Em(i)^2 + Em(i)*Em(j)*(Gf(i,j)*cos(dv(i)-dv(j)) + Bf(i,j)*sin(dv(i)-dv(j)));
%         end
%     end
% end

%% Aprovechamos y calculamos las potencias pre-falla que deben dar igual a lo que da el flujo de carga
Pe0Falla = zeros(n_gen, n_gen);
for i = 1:n_gen
    for j = 1:n_gen
        if(i ~= j)
            Pe0Falla(i, j) = (-Gf(i,j) + gf(i,j))*Em(i)^2 + Em(i)*Em(j)*(Gf(i,j)*cos(d0(i)-d0(j)) + Bf(i,j)*sin(d0(i)-d0(j)));
        end
    end
end