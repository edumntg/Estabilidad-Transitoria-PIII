%% Programa generico para estabilidad transitoria
clc, clear all

%% Los datos del archivo Excel se especifican de la siguiente manera
    % Hoja 1: Datos de barras
    % Hoja 2: Datos de ramas prefalla
    % Hoja 3: Datos de ramas falla
    % Hoja 4: Datos de ramas postfalla
    % Hoja 5: Generadores

DATAFILE = 'BUSDATA.xlsx';
% DATAFILE = 'Datos9barras.xlsx';

f = 60;
w0 = 2*pi*f;

ShowUnits = 0;
Vb = 115;
Sb = 100;

[BUSDATA, LINEDATA_PRE, LINEDATA_FALLA, LINEDATA_POST, GENDATA, SIMULATIONDATA] = ET_LoadData(DATAFILE);

%% No tocar
ti = 0;
f = SIMULATIONDATA(1);
tp = SIMULATIONDATA(2)/f;
td = SIMULATIONDATA(3)/f;
tf = SIMULATIONDATA(4)/f;
dt = SIMULATIONDATA(5)/f;
%%

n = size(BUSDATA, 1);                   % el numero de filas en el archivo excel es igual al numero de barras
nl_pre = size(LINEDATA_PRE, 1);         % el numero de filas en el archivo excel es igual al numero de ramas
nl_falla = size(LINEDATA_FALLA, 1);     % el numero de filas en el archivo excel es igual al numero de ramas
nl_post = size(LINEDATA_POST, 1);       % el numero de filas en el archivo excel es igual al numero de ramas
ng = size(GENDATA, 1);                  % el numero de filas en el archivo excel es igual al numero de generadores

[Ybusp, Gp, Bp, gp, bp] = ET_Ybus(LINEDATA_PRE, n);
[Ybusf, Gf, Bf, gf, bf] = ET_Ybus(LINEDATA_FALLA, n);
[Ybuspt, Gpt, Bpt, gpt, bpt] = ET_Ybus(LINEDATA_POST, n);

[V, theta, Pgen, Qgen, Pneta, Qneta, Sshunt, Pflow, Pflow_bus, ...
Qflow, Qflow_bus, Ploss, Qloss] = ET_FDC(BUSDATA, LINEDATA_PRE, Gp, Bp, gp, bp);

ET_PrintFDC;

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
            Vbarra = V(barra)*(cos(theta(barra)) + 1i*sin(theta(barra)));
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

[YbusExtp, YbusRMp] = ET_YbusExtendida(GENDATA, BUSDATA, V, Ybusp, 0);                     % Se arma la matriz extendida del sistema
[YbusExtf, YbusRMf] = ET_YbusExtendida(GENDATA, BUSDATA, V, Ybusf, 1);                     % Se arma la matriz extendida del sistema
[YbusExtpt, YbusRMpt] = ET_YbusExtendida(GENDATA, BUSDATA, V, Ybuspt, 0);                  % Se arma la matriz extendida del sistema

[YKronp, Yap, Ybp, Ycp, Ydp] = ET_Kron(YbusExtp, n_gen);                                % Se obtiene la matriz equivalente de Kron para el sistema pre-falla
[YKronf, Yaf, Ybf, Ycf, Ydf] = ET_Kron(YbusExtf, n_gen);                                % Se obtiene la matriz equivalente de Kron para el sistema pre-falla
[YKronpt, Yapt, Ybpt, Ycpt, Ydpt] = ET_Kron(YbusExtpt, n_gen);                                % Se obtiene la matriz equivalente de Kron para el sistema pre-falla

x0 = zeros(1, size(YKronp, 1));
x = fsolve(@(x)ET_KronShunt_FSOLVE(x, YKronp), x0);
YShuntKronp = ET_KronShunt(x);

x0 = zeros(1, size(YKronf, 1));
x = fsolve(@(x)ET_KronShunt_FSOLVE(x, YKronf), x0);
YShuntKronf = ET_KronShunt(x);

x0 = zeros(1, size(YKronpt, 1));
x = fsolve(@(x)ET_KronShunt_FSOLVE(x, YKronpt), x0);
YShuntKronpt = ET_KronShunt(x);

Pe0p = ET_Pe(Em, d0, YKronp, YShuntKronp);





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