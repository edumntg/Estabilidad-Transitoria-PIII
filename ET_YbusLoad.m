%% Eduardo Montilva 12-10089
% Script el cual tiene como funcion armar la Ybus del sistema, incluyendo
% los shunts conectados a las barras

% global n nl ns Ybus

function [Ybus, G, B, g, b] = ET_YbusLoad(BUSDATA, LINEDATA, FALLADATA, n, nl, V, event)

    
    Ybus = zeros(n,n);          % Se inicializa como una m


    g = zeros(n, n);
    b = zeros(n, n);

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
    
    %% Se agrega la impedancia de carga
    for i = 1:n
        bus = BUSDATA(i, 1);
        Sload = BUSDATA(i, 5) + 1i*BUSDATA(i, 6);
        if(abs(Sload) ~= 0)
            Zload = V(i)^2/conj(Sload);
            Ybus(bus, bus) = Ybus(bus, bus) + 1/Zload;
        end
    end
 
    if(event == 1) % se esta en condicion falla
        Zf = 1e-10 + 1i*1e-10;
        for i = 1:n
            if(FALLADATA(2) == i)
                Ybus(i, i) = Ybus(i, i ) + 1/Zf;
                
                % Si hay carga acoplada a esa barra, se elimina
                Sload = BUSDATA(i, 5) + 1i*BUSDATA(i, 6);
                if(abs(Sload) ~= 0)
                    Zload = V(i)^2/conj(Sload);
                    Ybus(i, i) = Ybus(i, i) - 1/Zload;
                end
                
                for j = 1:nl
                    if(LINEDATA(j, 1) == LINEDATA(j, 2) && LINEDATA(j, 1) == i) %shunt
                        Zshunt = 1i*LINEDATA(j, 4);
                        Ybus(i, i) = Ybus(i, i) - 1/Zshunt;
                    end
                end
            end
        end
    end
    
    G = real(Ybus);
    B = imag(Ybus);
end



