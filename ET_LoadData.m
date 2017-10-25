%% Carga de datos de barras y lineas desde archivo excel

%% Los datos del archivo Excel se especifican de la siguiente manera
    % Hoja 1: Datos de barras
    % Hoja 2: Datos de ramas prefalla
    % Hoja 3: Datos de ramas falla
    % Hoja 4: Datos de ramas postfalla
    % Hoja 5: Generadores
function [BUSDATA, LINEDATA_PRE, LINEDATA_FALLA, LINEDATA_POST, GENDATA, SIMULATIONDATA, FALLADATA] = ET_LoadData(DATAFILE)
    BUSDATA = xlsread(DATAFILE, 1);
    LINEDATA_PRE = xlsread(DATAFILE, 3);
    LINEDATA_FALLA = xlsread(DATAFILE, 4);
    LINEDATA_POST = xlsread(DATAFILE, 5);
    GENDATA = xlsread(DATAFILE, 6);
    FALLADATA = xlsread(DATAFILE, 7);
    SIMULATIONDATA = xlsread(DATAFILE, 8);
end