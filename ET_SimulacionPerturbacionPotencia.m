%% En este script se declara todo lo necesario para la simulacion de una perturbacion en potencia mecanica de una maquina

busi = input('Ingrese el Nro. de barra donde se ubica el generador que sufrira la perturbacion: ');
GEN = [];
NGen = 0;
GFound = 0;
for i = 1:n_gen
    if(GENDATA(i, 1) == busi)
        GEN = GENDATA(i, 1:size(GENDATA, 2));
        GFound = 1;
        NGen = i;
        break;
    end
end

if(GFound == 0)
    error('No existen un generador conectado a la barra especificada');
end

Pmnew = input('Ingrese (en p.u) el nuevo valor de la potencia mecanica para la maquina: ');

%% Para una variacion de la potencia mecanica, la topologia del sistema no cambia. 

YbusExtFalla = YbusExtPre;
YbusExtPost = YbusExtPre;

% Se calcula el Kron falla

[YKronFalla, YaF, YbF, YcF, YdF] = ET_Kron(YbusExtFalla, n_gen);

% Calculamos los shunts del sistema equivalente
x0 = zeros(1, size(YKronFalla, 1));
x = fsolve(@(x)ET_KronShunt_FSOLVE(x, YKronFalla), x0);

YShuntKronFalla = ET_KronShunt(x);

ET_PeExprFalla;

% Se calcula el Kron post-falla

[YKronPost, YaPt, YbPt, YcPt, YdPt] = ET_Kron(YbusExtPost, n_gen);

% Calculamos los shunts del sistema equivalente
x0 = zeros(1, size(YKronPost, 1));
x = fsolve(@(x)ET_KronShunt_FSOLVE(x, YKronPost), x0);

YShuntKronPost = ET_KronShunt(x);

ET_PeExprPost;


%% Para este tipo de falla, las potencias mecanicas SI varian

Pmpre = Pm;

Pmfalla = Pmpre;
Pmfalla(NGen) = Pmnew;

Pmpost = Pmfalla;