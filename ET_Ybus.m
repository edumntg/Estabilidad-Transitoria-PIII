%% Eduardo Montilva 12-10089
% Script el cual tiene como funcion armar la Ybus del sistema, incluyendo
% los shunts conectados a las barras

% global n nl ns Ybus

function [Ybus, G, B, g, b] = ET_Ybus(BUSDATA, LINEDATA, n, falla, V, theta)

    Ybus = zeros(n,n);          % Se inicializa como una matriz de ceros

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
    
    % Si nos encontramos en situacion de falla, debemos agregar la
    % impedancia de falla a la barra y eliminar las impedancias de shunts
    if(falla)
        for i = 1:n
            if(BUSDATA(i, 10) ~= 0) % barra con falla especificada
                
                % Asumimos en primera instancia que es una falla trifasica,
                % de esta forma se cortocircuitan las cargas y los shunts
                bus = BUSDATA(i, 1);
                Zfalla = BUSDATA(i, 11);
                
                % Se agrega Zfalla
                if(Zfalla ~= 0)
                	Ybus(bus, bus) = Ybus(bus, bus) + 1/Zfalla;
                end
                
                % Se elimina el shunt conectado
                for j = 1:size(LINEDATA, 1)
                    if(LINEDATA(j, 1) == LINEDATA(j, 2) && LINEDATA(j, 1) == bus)
                        % Es shunt y ademas esta conectado a la barra de
                        % estudio
                        Zshunt = 1i*LINEDATA(j, 4);
                        Ybus(bus, bus) = Ybus(bus, bus) - 1/Zshunt;
                    end
                end
                
                %% Si V y theta no estan vacios, significa que ya se realizo el FDC
                % Por tanto ya se pueden remover las impedancias de carga a las barras
                % especificas con falla
                if(isempty(V) == 0) % Vector voltajes no esta vacio
                    for i = 1:size(BUSDATA, 1)
                        bus = BUSDATA(i, 1);
                        Scarga = BUSDATA(i, 5) + 1i*BUSDATA(i, 6);
                        if(abs(Scarga ~= 0))
                            Zcarga = V(i)^2 /conj(Scarga);

                            Ybus(bus, bus) = Ybus(bus, bus) - 1/Zcarga;
                        end
                    end
                end
                
            end
        end
    end
    

    G = real(Ybus);
    B = imag(Ybus);
end



