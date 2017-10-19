% Eduardo Montilva 12-10089

% Programa para realizar la reduccion de Kron

function T = ET_TPARK(d, ng)
    %% Debemos reducir el sistema para eliminar las barras PQ.
    
   	T = zeros(2*ng, 2*ng);
    
    %% Ahora esta matriz la convertimos a RM
    for i = 1:ng
        for j = 1:ng
            if(i == j)
                T(2*i-1, 2*i-1) = cos(d(i));
                T(2*i, 2*i) = cos(d(i));
                
                T(2*i-1, 2*i) = sin(d(i));
                T(2*i, 2*i-1) = -sin(d(i));
            end
        end
    end
end