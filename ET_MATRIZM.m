% Eduardo Montilva 12-10089

% Programa para realizar la reduccion de Kron

function T = ET_MATRIZM(Ra, Xq, Xd, ng)
    %% Debemos reducir el sistema para eliminar las barras PQ.
    
    T = zeros(2*ng, 2*ng);
    
    %% Ahora esta matriz la convertimos a RM
    for i = 1:ng
        for j = 1:ng
            if(i == j)
                T(2*i-1, 2*i-1) = Ra(i);
                T(2*i, 2*i) = Ra(i);
                
                T(2*i-1, 2*i) = -Xd(i);
                T(2*i, 2*i-1) = Xq(i);
            end
        end
    end
end