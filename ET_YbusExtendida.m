%% Este script arma la matriz extendida del sistema

    % Las barras donde se ubican los generadores pasaran a ser las primeras
    % barras enumeradas y por ende los primeros elementos de la matriz
    % extendida
    
    % Las barras donde estan conectadas los generadores pasaran a ser las
    % primeras barras enumeradas
    
    % La matriz extendida RM va a ser de tama�o 2*n
function [YbusExt, YbusRM] = ET_YbusExtendida(GENDATA, BUSDATA, V, Ybus, falla)

    n = size(Ybus, 1);                   % tama�o original del sistema
    n_gen = size(GENDATA, 1);
    n_nuevo = n + n_gen;                    % nuevo tama�o del sistema: tama�o original + nro de generadores
    
    G = real(Ybus);
    B = imag(Ybus);
    
    YbusExt = zeros(n_nuevo, n_nuevo);
    for g = 1:n_gen
        Zgen = GENDATA(g, 2) + 1i*GENDATA(g, 3);
        YbusExt(g,g) = YbusExt(g,g) + 1/Zgen;

        YbusExt(g + n_gen, g + n_gen) = YbusExt(g + n_gen, g + n_gen) + 1/Zgen;

        YbusExt(g, n_gen + g) = -1/Zgen;
        YbusExt(n_gen + g, g) = -1/Zgen;
    end

    % Agregamos las impedancias de cargas conectadas a las barras
    for i = 1:n
        if(abs(BUSDATA(i, 5) + 1i*BUSDATA(i, 6)) ~= 0) %Si hay Pload, entonces hay carga en esa barra
            Zcarga(i) = (V(i)^2)/conj(BUSDATA(i, 5) + 1i*BUSDATA(i, 6));
            YbusExt(i + n_gen, i + n_gen) = YbusExt(i + n_gen, i + n_gen) + 1/Zcarga(i);
        end
    end
    
    % Si nos encontramos en situacion de falla, debemos agregar la
    % impedancia de falla a la barra y eliminar las impedancias de carga y
    % shunts
    if(falla)
        for i = 1:n
            if(GENDATA(i, 10) ~= 0)
                bus = GENDATA(i, 1);
                Zfalla = GENDATA(i, 11);
                if(Zfalla ~= 0)
                	YbusExt(bus + n_gen, bus + n_gen) = YbusExt(bus + n_gen, bus + n_gen) + 1/Zfalla;
                end
            end
        end
    end

    YbusExt((n_gen + 1):n_nuevo, (n_gen + 1):n_nuevo) = YbusExt((n_gen + 1):n_nuevo, (n_gen + 1):n_nuevo) + Ybus(1:n, 1:n);


    YbusRM = zeros(2*n, 2*n);
    
    % La matriz YbusRMPre estara conformada por [Ya Yb; Yc Ya] donde Ya
    % lleva los elementos propios
    
    Ya = zeros(n, n);
    Yb = zeros(n, n);
    for i = 1:n
        for j = 1:n
            Ya(i, j) = G(i, j);
            Yb(i, j) = -B(i, j);
        end
    end
    Yc = -Yb;
    
    YbusRM = [Ya Yb; Yc Ya];
end