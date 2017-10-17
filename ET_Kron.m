% Eduardo Montilva 12-10089

% Programa para realizar la reduccion de Kron

function [Ykron, Ya, Yb, Yc, Yd] = ET_Kron(Ybus, BUSDATA)
    %% Debemos reducir el sistema para eliminar las barras PQ.
    
    n = size(Ybus, 1);
    nPQ = 0;
    for i = 1:n
        if(BUSDATA(i, 2) == 0) %% Es una barra PQ
            nPQ = nPQ + 1;
        end
    end
    nPVS = n - nPQ;
    
    % La matriz Ya sera una matriz de tamaño = nPVS (tamaño igual al num de barras SLACK + PV)
    Ya = Ybus(1:nPVS, 1:nPVS);
    % La matriz Yb sera una matriz que va desde nPVS+1 hasta n
    Yb = Ybus(1:nPVS, nPVS+1:n);
    % La matriz Yc siempre es igual a la traspuesta de Yb
    Yc = transpose(Yb);
    % La matriz Yd es igual a una matriz que va desde nPVS+1 hasta n
    Yd = Ybus(nPVS+1:n, nPVS+1:n);
    
    Ykron = Ya - Yb*inv(Yd)*Yc;
end