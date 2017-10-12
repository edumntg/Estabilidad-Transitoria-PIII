%% Programa generico para estabilidad transitoria
clear all

%% Los datos del archivo Excel se especifican de la siguiente manera
    % Hoja 1: Datos de barras
    % Hoja 2: Datos de ramas prefalla
    % Hoja 3: Datos de ramas falla
    % Hoja 4: Datos de ramas postfalla
    % Hoja 5: Generadores

% DATAFILE = 'BUSDATA.xlsx';
DATAFILE = 'Datos9barras_confalla.xlsx';

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

[YbusExtp, YbusRMp] = ET_YbusExtendida(GENDATA, BUSDATA, LINEDATA_PRE, V, Ybusp, 0);                     % Se arma la matriz extendida del sistema
[YbusExtf, YbusRMf] = ET_YbusExtendida(GENDATA, BUSDATA, LINEDATA_FALLA, V, Ybusf, 1);                     % Se arma la matriz extendida del sistema
[YbusExtpt, YbusRMpt] = ET_YbusExtendida(GENDATA, BUSDATA, LINEDATA_POST, V, Ybuspt, 0);                  % Se arma la matriz extendida del sistema

[YKronp, Yap, Ybp, Ycp, Ydp] = ET_Kron(YbusExtp, n_gen);                                % Se obtiene la matriz equivalente de Kron para el sistema pre-falla
[YKronf, Yaf, Ybf, Ycf, Ydf] = ET_Kron(YbusExtf, n_gen);                                % Se obtiene la matriz equivalente de Kron para el sistema pre-falla
[YKronpt, Yapt, Ybpt, Ycpt, Ydpt] = ET_Kron(YbusExtpt, n_gen);                                % Se obtiene la matriz equivalente de Kron para el sistema pre-falla

x0 = zeros(1, size(YKronp, 1));
[x,fvalp,exitflagp,outputp,jacobianp] = fsolve(@(x)ET_KronShunt_FSOLVE(x, YKronp), x0);
YShuntKronp = ET_KronShunt(x);

x0 = zeros(1, size(YKronf, 1));
[x,fvalf,exitflagf,outputf,jacobianf] = fsolve(@(x)ET_KronShunt_FSOLVE(x, YKronf), x0);
YShuntKronf = ET_KronShunt(x);

x0 = zeros(1, size(YKronpt, 1));
[x,fvalpt,exitflagpt,outputpt,jacobianpt] = fsolve(@(x)ET_KronShunt_FSOLVE(x, YKronpt), x0);
YShuntKronpt = ET_KronShunt(x);

Pe0p = ET_Pe(Em, d0, YKronp, YShuntKronp);

%% Caso pre falla
w0 = zeros(1, ng);
[wp, dp, Pep] = ET_Integracion(H, Em, d0, w0, Pe0p, Pm, YKronp, YShuntKronp, f, [ti dt tp]);

%% Caso falla
for gi = 1:ng
    w0(gi) = wp(gi, length(wp));
    d0(gi) = dp(gi, length(dp));
    Pe0(gi) = Pep(gi, length(Pep));
end
[wf, df, Pef] = ET_Integracion(H, Em, d0, w0, Pe0, Pm, YKronf, YShuntKronf, f, [tp dt td]);
pause;
%% Caso post-falla
for gi = 1:ng
    w0(gi) = wf(gi, length(wf));
    d0(gi) = df(gi, length(df));
    Pe0(gi) = Pef(gi, length(Pef));
end
[wpt, dpt, Pept] = ET_Integracion(H, Em, d0, w0, Pe0, Pm, YKronpt, YShuntKronpt, f, [td dt tf]);

ET_Plot(wp, dp, Pep, wf, df, Pef, wpt, dpt, Pept, [ti tp td tf dt]);
% plot(ti:dt:tp, Pep(1, :), tp:dt:td, Pef(1, :), td:dt:tf, Pept(1, :)), line([tp tp], [Pef(1, 1) Pep(1, length(Pep))]), grid minor
