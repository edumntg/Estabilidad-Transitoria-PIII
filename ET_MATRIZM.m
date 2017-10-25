% Eduardo Montilva 12-10089

% Programa para realizar la reduccion de Kron

function M = ET_MATRIZM(Ra, Xq, Xd, ng)
    %% Debemos reducir el sistema para eliminar las barras PQ.
    
    M = zeros(2*ng, 2*ng);
    
    %% Ahora esta matriz la convertimos a RM
    for i = 1:ng
        for j = 1:ng
            if(i == j)
                M(2*i-1, 2*i-1) = Ra(i);
                M(2*i, 2*i) = Ra(i);
                
                M(2*i-1, 2*i) = -Xd(i);
                M(2*i, 2*i-1) = Xq(i);
            end
        end
    end
end