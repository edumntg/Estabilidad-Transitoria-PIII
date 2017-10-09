%% Eduardo Montilva 12-10089
%% Flujo de Carga 
%% Metodo de solucion: Newton Raphson
% clc, clear all;
% global REF SLACK PV PQ LINEDATA SHUNTDATA Pdesbalance

Vb = 115;         % Voltaje base (kV)
Sb = 100;         % Potencia base (MVA)

%% Tipos de barras
PQ = 0;
SLACK = 1;
PV = 2;
REF = 3;

ShowUnits = 0;    % Mostrar resultados en unidades? (MW, kV, MVAr) (1 = SI, 0 = NO) Si elige NO, se muestran en p.u

PrintResults = 1; % Imprimir resultados?
%% Declaracion data de shunt en barras
%           busid X MVAr
% SHUNTDATA = [];

%% Prueba con desbalance
Pdesbalance = -0.00;

%% Ejecucion del FDC
FDC_LoadData;
FDC_Ybus;
FDC_Solver;
if(PrintResults)
    FDC_Print;
end