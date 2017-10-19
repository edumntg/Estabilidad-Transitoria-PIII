%% Programa generico para estabilidad transitoria
clc, clear all, close all

%% Los datos del archivo Excel se especifican de la siguiente manera
    % Hoja 1: Datos de barras
    % Hoja 2: Datos de ramas prefalla
    % Hoja 3: Datos de ramas falla
    % Hoja 4: Datos de ramas postfalla
    % Hoja 5: Generadores

% DATAFILE = 'BUSDATA.xlsx';
% DATAFILE = 'Datos9barras_confalla.xlsx';
DATAFILE = 'BUSDATA_fallabarra.xlsx';
f = 60;
w0 = 2*pi*f;

ShowUnits = 0;
Vb = 115;
Sb = 100;

[BUSDATA, LINEDATA_PRE, LINEDATA_FALLA, LINEDATA_POST, GENDATA, SIMULATIONDATA, FALLADATA] = ET_LoadData(DATAFILE);

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

%% Aqui se cargaran en vectores los parametros de las maquinas como lo son resistencias, reactancias, constantes de
% tiempo y constantes de inercia
for i = 1:ng
    Ra(i, 1) = GENDATA(i, 6);
    Xq(i, 1) = GENDATA(i, 7);
    Xqp(i, 1) = GENDATA(i, 8);
    Xqpp(i, 1) = GENDATA(i, 9);
    Xd(i, 1) = GENDATA(i, 10);
    Xdp(i, 1) = GENDATA(i, 11);
    Xdpp(i, 1) = GENDATA(i, 12);
    Tq0p(i, 1) = GENDATA(i, 13);
    Td0p(i, 1) = GENDATA(i, 14);
    H(i, 1) = GENDATA(i, 15);
end


[Ybus_pre, G_pre, B_pre, g_pre, b_pre] = ET_Ybus(BUSDATA, LINEDATA_PRE, n, [], [], []);
% [Ybusf, Gf, Bf, gf, bf] = ET_Ybus(BUSDATA, LINEDATA_FALLA, n, [], [], []);
% [Ybuspt, Gpt, Bpt, gpt, bpt] = ET_Ybus(BUSDATA, LINEDATA_POST, n, [], [], []);

[V, theta, Pgen, Qgen, Pneta, Qneta, Sshunt, Pflow, Pflow_bus, ...
Qflow, Qflow_bus, Ploss, Qloss] = ET_FDC(BUSDATA, LINEDATA_PRE, G_pre, B_pre, g_pre, b_pre);

ET_PrintFDC;

%% Una vez que se realizo el FDC, se pueden modelar las cargas como impedancia y se agregan a la Ybus
[Ybusc_pre, Gc_pre, Bc_pre, gc_pre, bc_pre] = ET_Ybus(BUSDATA, LINEDATA_PRE, n, V, theta, []);
[Ybuscf, Gcf, Bcf, gcf, bcf] = ET_Ybus(BUSDATA, LINEDATA_FALLA, n, V, theta, FALLADATA);
[Ybuscpt, Gcpt, Bcpt, gcpt, bcpt] = ET_Ybus(BUSDATA, LINEDATA_POST, n, V, theta, []);

Pm = Pgen; 

%% Estos son los equivalentes de Kron sin impedancias de generadores. En este equivalente se eliminan barras PQ
[YKron_pre, YaR_pre, YbR_pre, YcR_pre, YdR_pre] = ET_Kron(Ybusc_pre, ng);                                % Se obtiene la matriz equivalente de Kron para el sistema pre-falla
[YKronRf, YaRf, YbRf, YcRf, YdRf] = ET_Kron(Ybuscf, ng);                                % Se obtiene la matriz equivalente de Kron para el sistema pre-falla
[YKronRpt, YaRpt, YbRpt, YcRpt, YdRpt] = ET_Kron(Ybuscpt, ng);                                % Se obtiene la matriz equivalente de Kron para el sistema pre-falla

%% Matrices RM
YKrm_pre = ET_YRM(YKron_pre);
YRMf = ET_YRM(YKronRf);
YRMpt = ET_YRM(YKronRpt);

%% Calculos de los VRM e IRM
[Vrm_pre, Irm_pre] = ET_VRMIRM(V, theta, YKrm_pre, BUSDATA);
[VRMf, IRMf] = ET_VRMIRM(V, theta, YRMf, BUSDATA);
[VRMp, IRMpt] = ET_VRMIRM(V, theta, YRMpt, BUSDATA);

%% Calculos de las tensiones internas ficticas, a fin de obtener los angulos delta para cada generador
[Efp, d0p] = ET_EFICT(Vrm_pre, Irm_pre, GENDATA);

%% Calculo de la potencia mecanica de los generadores
Pm_pre = ET_PMEC(GENDATA, Vrm_pre, Irm_pre);

%% Calculo de las potencias electricas pre-falla en bornes de los generadores
Pe0p = ET_Pe(GENDATA, Pm_pre, Irm_pre);

%% Calculo de la matriz T de Park para calcular la Eexc en el momento pre falla
Mp = ET_MATRIZM(Ra, Xqp, Xdp, ng);
T = ET_TPARK(d0p, ng);
A = inv(T)*Mp*T;

%% Calculo de tensiones y corrientes qd
[Vqdp, Iqd, Vqp, Vdp, Iq, Id] = ET_VIQD(Vrm_pre, Irm_pre, T, ng);
[Eqp, Edp] = ET_EQD(Vqp, Vdp, Iq, Id, Ra, Xdp, Xqp, ng);
Eexc = ET_EEXC(Eqp, Edp);
%% Se empieza con las integraciones
% Se integrara por separado: prefalla, falla y postfalla
   
% Se definen los valores iniciales para caso pre-falla
wo = zeros(ng, 1);
do = d0p;
Peo = Pe0p;

[w_pre, d_pre, Pe_pre, Eqp_pre, Edp_pre, Vt_pre] = ET_Integracion(wo, do, Peo, Iq, Id, Ra, Xq, Xqp, Xd, Xdp, Eqp, Edp, Tq0p, Td0p, Pm_pre, YKrm_pre, Mp, Eexc, w0, H, ng, [ti dt tp]);

%% Caso falla
for gi = 1:ng
    w0(gi) = wp(gi, length(wp));
    d0(gi) = dp(gi, length(dp));
    Pe0(gi) = Pep(gi, length(Pep));
end
[wf, df, Pef] = ET_Integracion(H, Em, d0, w0, Pe0, Pm, YKronf, YShuntKronf, f, [tp dt td]);

%% Caso post-falla
for gi = 1:ng
    w0(gi) = wf(gi, length(wf));
    d0(gi) = df(gi, length(df));
    Pe0(gi) = Pef(gi, length(Pef));
end
[wpt, dpt, Pept] = ET_Integracion(H, Em, d0, w0, Pe0, Pm, YKronpt, YShuntKronpt, f, [td dt tf]);

ET_Plot(wp, dp, Pep, wf, df, Pef, wpt, dpt, Pept, [ti tp td tf dt]);
% plot(ti:dt:tp, Pep(1, :), tp:dt:td, Pef(1, :), td:dt:tf, Pept(1, :)), line([tp tp], [Pef(1, 1) Pep(1, length(Pep))]), grid minor
