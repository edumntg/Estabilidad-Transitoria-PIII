% Eduardo Montilva 12-10089

% Programa para realizar la reduccion de Kron

function [Ykron, Ya, Yb, Yc, Yd] = ET_KronExtendida(Ybusext, ng)
    %% Debemos reducir el sistema para eliminar las barras PQ.
    
    n = size(Ybusext, 1);
    
    % La matriz Ya sera una matriz de tamaño = ng (numero de gen)
    Ya = Ybusext(1:ng, 1:ng);
    
    % La matriz Yb sera una matriz que va desde ng+1 hasta n
    Yb = Ybusext(1:ng, ng+1:n);
    
    % La matriz Yc siempre es igual a la traspuesta de Yb
    Yc = transpose(Yb);
    
    % La matriz Yd es igual a una matriz que va desde ng+1 hasta n
    
    Yd = Ybusext(ng+1:n, ng+1:n);
    
    Ykron = Ya - Yb*inv(Yd)*Yc;
end