%% Este script arma la matriz extendida del sistema

    % Las barras donde se ubican los generadores pasaran a ser las primeras
    % barras enumeradas y por ende los primeros elementos de la matriz
    % extendida
    
    % Las barras donde estan conectadas los generadores pasaran a ser las
    % primeras barras enumeradas
    
%     n_orig = size(Ybus, 1);                 % tamaño original del sistema
    n_nuevo = n + n_gen;               % nuevo tamaño del sistema: tamaño original + nro de generadores
    YbusExtPre = zeros(n_nuevo, n_nuevo);
    for g = 1:n_gen
        Zgen = GENDATA(g, 2) + 1i*GENDATA(g, 3);
        YbusExtPre(g,g) = YbusExtPre(g,g) + 1/Zgen;

        YbusExtPre(g + n_gen, g + n_gen) = YbusExtPre(g + n_gen, g + n_gen) + 1/Zgen;

        YbusExtPre(g, n_gen + g) = -1/Zgen;
        YbusExtPre(n_gen + g, g) = -1/Zgen;
    end

    % Agregamos las impedancias de cargas conectadas a las barras
    for i = 1:n
        if(Pload(i) ~= 0) %Si hay Pload, entonces hay carga en esa barra
            Zcarga(i) = (V(i)^2)/conj(abs(Pload(i)) + 1i*abs(Qload(i)));
            YbusExtPre(i + n_gen, i + n_gen) = YbusExtPre(i + n_gen, i + n_gen) + 1/Zcarga(i);
        end
    end

    YbusExtPre((n_gen + 1):n_nuevo, (n_gen + 1):n_nuevo) = YbusExtPre((n_gen + 1):n_nuevo, (n_gen + 1):n_nuevo) + Ybus(1:n, 1:n);
