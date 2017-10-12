%% Este script arma la matriz extendida del sistema

    % Las barras donde se ubican los generadores pasaran a ser las primeras
    % barras enumeradas y por ende los primeros elementos de la matriz
    % extendida
    
    % Las barras donde estan conectadas los generadores pasaran a ser las
    % primeras barras enumeradas
    
    % La matriz extendida RM va a ser de tamaño 2*n
function [YbusExt, YbusRM] = ET_YbusExtendida(GENDATA, BUSDATA, LINEDATA, V, Ybus, falla)

    n = size(Ybus, 1);                   % tamaño original del sistema
    n_gen = size(GENDATA, 1);
    n_nuevo = n + n_gen;                    % nuevo tamaño del sistema: tamaño original + nro de generadores
    
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

    YbusExt((n_gen + 1):n_nuevo, (n_gen + 1):n_nuevo) = YbusExt((n_gen + 1):n_nuevo, (n_gen + 1):n_nuevo) + Ybus(1:n, 1:n);

    % Si nos encontramos en situacion de falla, debemos agregar la
    % impedancia de falla a la barra y eliminar las impedancias de carga y
    % shunts
    if(falla)
        for i = 1:n
            if(BUSDATA(i, 10) ~= 0) % barra con falla especificada
                
                % Asumimos en primera instancia que es una falla trifasica,
                % de esta forma se cortocircuitan las cargas y los shunts
                bus = BUSDATA(i, 1);
                Zfalla = BUSDATA(i, 11);
                
                % Se agrega Zfalla
                if(Zfalla ~= 0)
                	YbusExt(bus + n_gen, bus + n_gen) = YbusExt(bus + n_gen, bus + n_gen) + 1/Zfalla;
                end
                
                %Se elimina Zcarga
                
                Scarga = BUSDATA(i, 5) + 1i*BUSDATA(i, 6);
                if(Scarga ~= 0)
                    Zcarga = V(bus)^2 /conj(Scarga);
                    YbusExt(bus + n_gen, bus + n_gen) = YbusExt(bus + n_gen, bus + n_gen) - 1/Zcarga;
                end
                
                % Se elimina el shunt conectado
                for j = 1:size(LINEDATA, 1)
                    if(LINEDATA(j, 1) == LINEDATA(j, 2) && LINEDATA(j, 1) == bus)
                        % Es shunt y ademas esta conectado a la barra de
                        % estudio
                        Zshunt = 1i*LINEDATA(j, 4);
                        YbusExt(bus + n_gen, bus + n_gen) = YbusExt(bus + n_gen, bus + n_gen) - 1/Zshunt;
                    end
                end
            end
        end
    end

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