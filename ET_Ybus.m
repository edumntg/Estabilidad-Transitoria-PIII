%% Eduardo Montilva 12-10089
% Script el cual tiene como funcion armar la Ybus del sistema, incluyendo
% los shunts conectados a las barras

% global n nl ns Ybus

function [Ybus, G, B, g, b] = ET_Ybus(BUSDATA, LINEDATA, n, V, theta, FALLADATA)

    TIPO_FALLA = -1;
    Ybus = zeros(n,n);          % Se inicializa como una matriz de ceros
    if(isempty(FALLADATA) == 0) % Hay una falla especificada
        TIPO_FALLA = FALLADATA(1);
        FALLA_FASES = FALLADATA(2);
        FALLA_Bi = FALLADATA(3);
        FALLA_Bj = FALLADATA(4);
        if(TIPO_FALLA == 1) % Falla longitudinal
            Ybus = zeros(n+1,n+1); % Si la falla es longitudinal, la matriz Ybus tendra una barra adicional
        end
    end


    g = zeros(n, n);
    b = zeros(n, n);
    
    nl = size(LINEDATA, 1);

    % Se agregan las impedancias de las lineas
    for i = 1:nl 
        from = LINEDATA(i, 1);                      % Barra de inicio de la linea
        to = LINEDATA(i, 2);                        % Barra de fin de la linea
        Zl = LINEDATA(i, 3) + 1i*LINEDATA(i, 4);    % Impedancia de la linea
        a = LINEDATA(i, 7);                         % Tap del transformador
        Bl = 1i*LINEDATA(i, 5);                     % Shunt x2 de la linea

        Ybus(from, from) = Ybus(from, from) + (1/Zl)/(a^2);
        if(from ~= to)
            Ybus(to, to) = Ybus(to, to) + (1/Zl)/(a^2);
            Ybus(from, to) = Ybus(from, to) - (1/Zl)/a;
            Ybus(to, from) = Ybus(to, from) - (1/Zl)/a;
        end

        % Se agregan los shunts de las lineas
        Ybus(from, from) = Ybus(from, from) + Bl/2;

        if(from ~= to)
            Ybus(to, to) = Ybus(to, to) + Bl/2;
            g(from, to) = real(Bl/2);
            g(to, from) = real(Bl/2);

            b(from, to) = imag(Bl/2);
            b(to, from) = imag(Bl/2);
        end
    end
    
    %% Si V y theta no estan vacios, significa que ya se realizo el FDC
     % Por tanto ya se pueden agregar las impedancias de carga a las barras
     % especificas
    if(isempty(V) == 0) % Vector voltajes no esta vacio
        for i = 1:size(BUSDATA, 1)
            bus = BUSDATA(i, 1);
            Scarga = BUSDATA(i, 5) + 1i*BUSDATA(i, 6);
            if(abs(Scarga ~= 0))
                Zcarga = V(i)^2 /conj(Scarga);

                Ybus(bus, bus) = Ybus(bus, bus) + 1/Zcarga;
            end
        end
    end
    
    % Si se definio un tipo de falla entonces hay que aplicar
    if(TIPO_FALLA ~= -1) % Falla definida
        Zfalla = 1e-10;
        if(TIPO_FALLA == 0) % Falla en barra
            bus = FALLA_Bi;
            % Agregamos la impedancia de falla a la barra fallada, y
            % borramos las impedancias de carga y shunts conectadas a esa
            % barra

            Ybus(bus, bus) = Ybus(bus, bus) + 1/Zfalla;
            %% Eliminamos impedancia de carga
            for i = 1:size(BUSDATA, 1)
                if(bus == BUSDATA(i, 1))
                    Scarga = BUSDATA(i, 5) + 1i*BUSDATA(i, 6);
                    if(abs(Scarga ~= 0))
                        Zcarga = V(i)^2 /conj(Scarga);

                        Ybus(bus, bus) = Ybus(bus, bus) - 1/Zcarga;
                    end
                    break;
                end
            end
            
            %% Se elimina el shunt conectado
            for j = 1:size(LINEDATA, 1)
                if(LINEDATA(j, 1) == LINEDATA(j, 2) && LINEDATA(j, 1) == bus)
                    % Es shunt y ademas esta conectado a la barra de
                    % estudio
                    Zshunt = 1i*LINEDATA(j, 4);
                    Ybus(bus, bus) = Ybus(bus, bus) - 1/Zshunt;
                end
            end
            
        %% Falla longitudinal
        elseif(TIPO_FALLA == 1)
            % Ya se ha definido la barra adicional y las impedancias para
            % esa barra en el archivo Excel, solo nos queda agregar la
            % impedancia de falla a la barra
            
            % Esta barra de falla siempre se encontrara en el ultimo
            % elemento de la matriz Ybus
            ny = size(Ybus, 1);
            Ybus(ny, ny) = Ybus(ny, ny) + 1/Zfalla;
            
        end
    end

    G = real(Ybus);
    B = imag(Ybus);
end



