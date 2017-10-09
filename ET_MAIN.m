%% Programa generico para estabilidad transitoria
clc, clear all

DATAFILE = 'BUSDATA.xlsx';
% DATAFILE = 'Datos9barras.xlsx';

f = 60;
w0 = 2*pi*f;

%% Tiempos de perturbacion y despeje
dt = 0.5/f;
ti = 0;
tp = 4/f;
td = 18/f;
tf = 250/f;

%%

FDC_LoadData;
FDC_Main_Excel;

Pm = Pgen; 

%% Calculos de tensiones internas de los generadores
n_gen = size(GENDATA, 1);
for g = 1:n_gen
    
    H(g) = GENDATA(g, 6);
    Em(g) = GENDATA(g, 4);
    Eang(g) = GENDATA(g, 5);
    E(g) = Em(g)*(cos(Eang(g)) + 1i*sin(Eang(g)));
    
    if(Em(g) == -1) % No posee tension interna declarada
        Sgen = (Pgen(g) + 1i*Qgen(g));
        if(abs(Sgen) > 0)
            Zgen = GENDATA(g, 2) + 1i*GENDATA(g, 3);
            barra = GENDATA(g, 1);
            Vbarra = V(barra)*(cos(d(barra)) + 1i*sin(d(barra)));
            Ibarra(g) = conj(Sgen/Vbarra);
            Ef = Vbarra + Ibarra(g)*Zgen;
            GENDATA(g, 4) = abs(Ef);
            GENDATA(g, 5) = angle(Ef);
            d0(g) = angle(Ef);
            E(g) = Ef;
            Em(g) = abs(Ef);
            Eang(g) = angle(Ef);
        end
    end
end

ET_YbusExtendida;                           % Se arma la matriz extendida del sistema
[YKronPre, YaPre, YbPre, YcPre, YdPre] = ET_Kron(YbusExtPre, n_gen);         % Se obtiene la matriz equivalente de Kron para el sistema pre-falla

x0 = zeros(1, size(YKronPre, 1));
x = fsolve(@(x)ET_KronShunt_FSOLVE(x, YKronPre), x0);

YShuntKronPre = ET_KronShunt(x);

ET_PeExprPre;

%% Una vez obtenidas las expresiones pre-falla, pasamos a simular la falla.
tipo_perturbacion = input('Ingrese el tipo de perturbacion (0 = Potencia Mecanica, 1 = Falla): ');
if(tipo_perturbacion == 0)
    ET_SimulacionPerturbacionPotencia;
else
    tipo_falla = input('Ingrese el tipo de falla (1 = longitudinal, 2 = barra): ');
    if(tipo_falla == 1)
        ET_SimulacionFallaLongitudinal;
    elseif(tipo_falla == 2)
        ET_SimulacionFallaBarra;
    end
end
% Por ahora, se asume que la falla se despeja y la linea queda en servicio,
% por lo que las ecuaciones de potencia post-falla seran iguales a la
% pre-falla

% PeEqPost = PeEqPre;
ET_Oscilacion;

ET_Plot;