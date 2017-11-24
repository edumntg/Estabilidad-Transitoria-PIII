%% Carga de datos de barras y lineas desde archivo excel

%% Los datos del archivo Excel se especifican de la siguiente manera
    % Hoja 1: Datos de barras
    % Hoja 2: Datos de ramas prefalla
    % Hoja 3: Datos de ramas falla
    % Hoja 4: Datos de ramas postfalla
    % Hoja 5: Generadores
function [BUSDATA, LINEDATA_PRE, LINEDATA_FALLA, LINEDATA_POST, GENDATA, SIMULATIONDATA, FALLADATA] = ET_LoadData_DAT()
    BUSDATA = load('BUSDATA.dat', '-ascii');
    LINEDATA_PRE = load('RAMAS_PRE.dat', '-ascii');
    LINEDATA_FALLA = load('RAMAS_FALLA.dat', '-ascii');
    LINEDATA_POST = load('RAMAS_POST.dat', '-ascii');
    GENDATA = load('GENERADORES.dat', '-ascii');
    FALLADATA = load('FALLA.dat', '-ascii');
    SIMULATIONDATA = load('INTEGRACION.dat', '-ascii');
end