%% Este script calculara las expresiones de la ecuacion de potencia para los generadores

function Pe = ET_Pe(E, d, YKron, YShuntKron)

    %% Aprovechamos y calculamos las potencias pre-falla que deben dar igual a lo que da el flujo de carga
    
    G = real(YKron);
    B = imag(YKron);
    g = real(YShuntKron);
    b = imag(YShuntKron);
    
    n_gen = length(E);
    Pe0m = zeros(n_gen, n_gen);
    Pe = zeros(1, n_gen);
    for i = 1:n_gen
        for j = 1:n_gen
            if(i ~= j)
                Pe0m(i, j) = (-G(i,j) + g(i,j))*E(i)^2 + E(i)*E(j)*(G(i,j)*cos(d(i)-d(j)) + B(i,j)*sin(d(i)-d(j)));
            end
        end
    end

    for i = 1:n_gen
        Pe(i) = sum(Pe0m(i, 1:n_gen));
    end
end