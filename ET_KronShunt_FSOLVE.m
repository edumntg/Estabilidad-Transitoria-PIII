%% Se tendran tantos shunts como lineas en el equivalente de Kron

function F = ET_KronShunt_FSOLVE(x, Ykron)

    % Habra tantas ecuaciones como lineas en el sistema equivalente
    % Habra tantas variables como lineas en el sistema equivalente
    
    % Ejemplo: En un sistema equivalente de 3 barras, habra 3 lineas
    % Y por tanto habra 3 variables (3 shunts) que son: B12, B13 y B23
    
    % Ejemplo: En un sistema equivalente de 2 barras, habra una linea
    % Y por tanto habra dos variable B12 y B21
    

    ns = size(Ykron, 1);
    Shuntij = zeros(1, ns);
    kv = 1;
    kf = 1;
    if(ns < 3)
        for i = 1:ns
            for j = 1:ns
                if(i ~= j)
                    Shuntij(i, j) = x(kv);
                    kv = kv + 1;
                end
            end
        end
        
        for i = 1:ns
            F(kf) = Ykron(i, i);
            for j = 1:ns
                if(i ~= j)
                    F(kf) = F(kf) + Ykron(i, j) - Shuntij(i, j);
                end
            end
            kf = kf + 1;
        end
    else
        for i = 1:ns
            for j = 1:ns
                if(i ~= j && i < j)
                    Shuntij(i, j) = x(kv);
                    Shuntij(j, i) = x(kv);
                    kv = kv + 1;
                end
            end
        end
        
        for i = 1:ns
            F(kf) = Ykron(i, i);
            for j = 1:ns
                if(i ~= j)
                    F(kf) = F(kf) + Ykron(i, j) - Shuntij(i, j);
                end
            end
            kf = kf + 1;
        end
    end
end