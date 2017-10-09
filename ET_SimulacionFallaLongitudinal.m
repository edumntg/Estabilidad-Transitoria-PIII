%% En este script se declarara todo lo necesario para la simulacion de una falla

    % Primero se probara con fallas longitudinales a mitad de linea
    
lni = input('Ingrese la barra de partida de la linea donde ocurre la falla: ');
lnj = input('Ingrese la barra de llegada de la linea donde ocurre la falla: ');
lnperc = input('Ingrese a que porcentaje de la linea (desde la barra de partida) ocurre la falla: ');
lout = input('Ingrese si la linea sale al despejar la filla (1: si, 0: no): ');
lnperc = abs(lnperc);
lnperc = min(0.9999, lnperc);
lnperc = max(0.0001, lnperc);
LineaFalla = 0;

for i = 1:size(LINEDATA, 1)
    if((LINEDATA(i, 1) == lni && LINEDATA(i, 2) == lnj) || (LINEDATA(i, 1) == lnj && LINEDATA(i, 2) == lni))
        LineaFalla = LINEDATA(i, 1:size(LINEDATA, 2));
    end
end
if(LineaFalla == 0)
    error('No existe linea entre las barras especificadas');
end

Zlinea = LineaFalla(3) + 1i*LineaFalla(4);
Zlineai = Zlinea*lnperc;
Zlineaj = Zlinea*(1-lnperc);

Blinea = 1i*LineaFalla(5);
%% Al ocurrir la falla, a la Ybus extendida se le agrega una barra adicional que corresponde al punto de la falla
nextpre = size(YbusExtPre, 1);
nextfalla = nextpre + 1;
YbusExtFalla = zeros(nextfalla, nextfalla);
YbusExtFalla(1:nextpre, 1:nextpre) = YbusExtPre(1:nextpre, 1:nextpre);

% A los elementos diagonales de la Ybus extendida pre-falla hay que
% eliminarles la impedancia de la linea fallada y agregarle la mitad de esa
% misma impedancia

% El numero del elemento en la diagonal de la matriz correspondiente a la
% linea sera el numero de la barra de la linea mas el numero de generadores
% (Recordando que el sistema se extiende agregando una barra adicional por
% generador)

YbusExtFalla(lni+n_gen, lni+n_gen) = YbusExtFalla(lni+n_gen, lni+n_gen) - 1/Zlinea;
YbusExtFalla(lni+n_gen, lni+n_gen) = YbusExtFalla(lni+n_gen, lni+n_gen) + 1/(Zlineai);
YbusExtFalla(lnj+n_gen, lnj+n_gen) = YbusExtFalla(lnj+n_gen, lnj+n_gen) - 1/Zlinea;
YbusExtFalla(lnj+n_gen, lnj+n_gen) = YbusExtFalla(lnj+n_gen, lnj+n_gen) + 1/(Zlineaj);

% Ahora, la impedancia mutua entre las barras i e j se hace cero
YbusExtFalla(lni+n_gen, lnj+n_gen) = 0;
YbusExtFalla(lnj+n_gen, lni+n_gen) = 0;

% Ahora, esta mutua se agrega entre barras i e nueva-barra, j y nueva-barra
% (nota: la mitad)
YbusExtFalla(lni+n_gen, nextfalla) = -1/(Zlineai);
YbusExtFalla(nextfalla, lni+n_gen) = -1/(Zlineai);
YbusExtFalla(lnj+n_gen, nextfalla) = -1/(Zlineaj);
YbusExtFalla(nextfalla, lnj+n_gen) = -1/(Zlineaj);

% A la diagonal de la nueva barra se agrega la impedancia de la linea
% completa y los shunts (en este casi si se agrega completo el shunt)
% Tambien hay que agregar la impedancia de falla

Zfalla = 1e-10;

YbusExtFalla(nextfalla, nextfalla) = (1/Zlineai) + (1/Zlineaj) + 2*Blinea + 1/Zfalla;

% Se calcula la Kron de falla

[YKronFalla, YaF, YbF, YcF, YdF] = ET_Kron(YbusExtFalla, n_gen);

% Calculamos los shunts del sistema equivalente
x0 = zeros(1, size(YKronFalla, 1));
x = fsolve(@(x)ET_KronShunt_FSOLVE(x, YKronFalla), x0);

YShuntKronFalla = ET_KronShunt(x);

ET_PeExprFalla;

% Si la topologia del sistema post falla es igual al pre falla, su Ybus y
% Kron seran iguales a las pre-falla
YbusExtPost = YbusExtPre;
[YKronPost, YaPost, YbFPost, YcPost, YdPost] = ET_Kron(YbusExtPost, n_gen);

%% Si luego del despeje se decide despejar la linea, hay que calcular la topologia post falla
if(lout == 1)
    %% Al despejar la linea, la matriz YbusExtendida-PostFalla sera igual a la YBusExtendida-PreFalla, pero sin la linea fallada
    %  Esto significa que los elementos diagonales de la matriz PreFalla
    %  donde estaba ubicada la linea ya no tendran la impedancia de esa
    %  linea, ni sus shunts, ademas de que no habra mutua entre estas
    %  barras (Se elimina el elemento mutuo de la linea)
    
    YbusExtPost = YbusExtPre;
    YbusExtPost(lni+n_gen, lni+n_gen) = YbusExtPost(lni+n_gen, lni+n_gen) - 1/Zlinea - Blinea;
    YbusExtPost(lnj+n_gen, lnj+n_gen) = YbusExtPost(lnj+n_gen, lnj+n_gen) - 1/Zlinea - Blinea;
    
    % Se elimina la mutua, como en las mutuas se agrega de forma negativa,
    % para eliminarla hay que sumarla
    YbusExtPost(lni+n_gen, lnj+n_gen) = YbusExtPost(lni+n_gen, lnj+n_gen) + 1/Zlinea;
    YbusExtPost(lnj+n_gen, lni+n_gen) = YbusExtPost(lnj+n_gen, lni+n_gen) + 1/Zlinea;
    
    % Se calcula la Kron post-falla

    [YKronPost, YaPost, YbFPost, YcPost, YdPost] = ET_Kron(YbusExtPost, n_gen);
    
end

% Calculamos los shunts del sistema equivalente
x0 = zeros(1, size(YKronPost, 1));
x = fsolve(@(x)ET_KronShunt_FSOLVE(x, YKronPost), x0);

YShuntKronPost = ET_KronShunt(x);

ET_PeExprPost;

%% Para este tipo de falla, las potencias mecanicas no varian

Pmpre = Pm;
Pmfalla = Pm;
Pmpost = Pm;