%% Carga de datos de barras y lineas desde archivo excel

[BDataNum, BDataStr, BDataTODO] = xlsread(DATAFILE, 'BUSDATA');
[LDataNum, LDataStr, LDataTODO] = xlsread(DATAFILE, 'LINEDATA');
[GDataNum, GDataStr, GDataTODO] = xlsread(DATAFILE, 'GENDATA');
[SDataNum, SDataStr, SDataTODO] = xlsread(DATAFILE, 'SHUNTDATA');
n = size(BDataNum, 1);      % el numero de filas en el archivo excel es igual al numero de barras
nl = size(LDataNum, 1);     % el numero de filas en el archivo excel es igual al numero de ramas
ng = size(GDataNum, 1);  % el numero de filas en el archivo excel es igual al numero de generadores
ns = size(SDataNum, 1);  % el numero de filas en el archivo excel es igual al numero de shunts
BUSDATA = zeros(n, 9);
LINEDATA = zeros(nl, 6);
GENDATA = zeros(ng, 6);
SHUNTDATA = [];
for i = 1:n
    BUSDATA(i, 1) = BDataNum(i, 1);                              % id de barra
    BUSDATA(i, 2) = BDataNum(i, 2);                              % tipo de barra
    BUSDATA(i, 3) = BDataNum(i, 3);                              % tension
    BUSDATA(i, 4) = BDataNum(i, 4);                              % angulo
    BUSDATA(i, 5) = BDataNum(i, 5);                              % Pdemandada
    BUSDATA(i, 6) = BDataNum(i, 6);                              % Qdemandada
    BUSDATA(i, 7) = BDataNum(i, 7);                              % Pgenerada
    BUSDATA(i, 8) = BDataNum(i, 8);                              % Qgenerada
    BUSDATA(i, 9) = BDataNum(i, 9);                              % Factor de participacion Ki
end

for i = 1:nl
    LINEDATA(i, 1) = LDataNum(i, 1);                             % Barra de partida
    LINEDATA(i, 2) = LDataNum(i, 2);                             % Barra de llegada
    LINEDATA(i, 3) = LDataNum(i, 3);                             % Resistencia
    LINEDATA(i, 4) = LDataNum(i, 4);                             % Reactancia
    LINEDATA(i, 5) = LDataNum(i, 5)/2;                             % Shunt
    LINEDATA(i, 6) = LDataNum(i, 7);                             % Tap
end

for i = 1:ng
    GENDATA(i, 1) = GDataNum(i, 1);                             % Barra
    GENDATA(i, 2) = GDataNum(i, 2);                             % Resistencia Interna
    GENDATA(i, 3) = GDataNum(i, 3);                             % Reactancia Transitoria
    GENDATA(i, 4) = GDataNum(i, 4);                             % Tension Interna
    GENDATA(i, 5) = GDataNum(i, 5);                             % Angulo d0
    GENDATA(i, 6) = GDataNum(i, 6);                             % H
end

if(ns > 0)
    SHUNTDATA = zeros(ns, 2);
end
for i = 1:ns
    SHUNTDATA(i, 1) = SDataNum(i, 1);                             % Barra
    SHUNTDATA(i, 2) = SDataNum(i, 2);                             % Reactancia
end