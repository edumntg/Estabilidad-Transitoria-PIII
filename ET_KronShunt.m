% Eduardo Montilva 12-10089

% Funcion para determinar los shunts del sistema equivalente Kron

function Y = ET_KronShunt(x)
    kv = 1;
    n = length(x);
    YShuntKron = zeros(n, n);
    if(n < 3)
        for i = 1:n
            for j = 1:n
                if(i ~= j)
                    YShuntKron(i, j) = x(kv);
                    kv = kv + 1;
                end
            end
        end
    else
        for i = 1:n
            for j = 1:n
                if(i ~= j && i < j)
                    YShuntKron(i, j) = x(kv);
                    YShuntKron(j, i) = x(kv);
                    kv = kv + 1;
                end
            end
        end
    end
    Y = YShuntKron;
end


