%% En este script se declarara todo lo necesario para la simulacion de una falla

    %Falla 3F en barra
    
bi = input('Ingrese la barra donde ocurre la falla: ');
barrai = bi + n_gen;
%% La barra de falla se mueve a la ultima fila y columna de la matriz extendida
ny = size(YbusExtPre, 1);
fila_mover = YbusExtPre(barrai, 1:ny);
fila_sust = YbusExtPre(ny, 1:ny);

YbusExtPre(barrai, 1:ny) = fila_sust;
YbusExtPre(ny, 1:ny) = fila_mover;

%% Al ultimo elemento diagonal (correspondiente a la barra fallada) se agrega la impedancia de falla

Zf = 1e-10;

YbusExtFalla = YbusExtPre;
YbusExtFalla(barrai, barrai) = YbusExtFalla(barrai, barrai) + 1/Zf;

% Ya que se cortocircuita la carga en esa barra (si existe) y los shunt (si
% existen), estos se eliminan del elemento de la diagonal correspondiente a la barra

if(abs(Zcarga(bi)) ~= 0) % Existe impedancia de carga conectada a esa barra
    YbusExtFalla(barrai, barrai) = YbusExtFalla(barrai, barrai) - 1/Zcarga(bi);
end

if(isempty(SHUNTDATA) == 0)
    for i = 1:size(SHUNTDATA, 1)
        if(SHUNTDATA(i, 1) == bi) % Existe shunt conectado a esa barra
            Zshunt = 1i*SHUNTDATA(i, 2);
            YbusExtFalla(barrai, barrai) = YbusExtFalla(barrai, barrai) - 1/Zshunt;
        end
    end
end

% Se calcula el Kron falla

[YKronFalla, YaF, YbF, YcF, YdF] = ET_Kron(YbusExtFalla, n_gen);

% Calculamos los shunts del sistema equivalente
x0 = zeros(1, size(YKronFalla, 1));
x = fsolve(@(x)ET_KronShunt_FSOLVE(x, YKronFalla), x0);

YShuntKronFalla = ET_KronShunt(x);

ET_PeExprFalla;

% Se asume que la topologia del sistema post-falla es igual al pre-falla

YbusExtPost = YbusExtPre;
[YKronPost, YaPost, YbFPost, YcPost, YdPost] = ET_Kron(YbusExtPost, n_gen);

% Calculamos los shunts del sistema equivalente
x0 = zeros(1, size(YKronPost, 1));
x = fsolve(@(x)ET_KronShunt_FSOLVE(x, YKronPost), x0);

YShuntKronPost = ET_KronShunt(x);

ET_PeExprPost;

%% Para este tipo de falla, las potencias mecanicas no varian

Pmpre = Pm;
Pmfalla = Pm;
Pmpost = Pm;