%% Eduardo Montilva 12-10089
% Script el cual tiene como funcion armar la Ybus del sistema, incluyendo
% los shunts conectados a las barras

% global n nl ns Ybus

n = size(BUSDATA, 1);       % Tamaño de la Ybus
nl = size(LINEDATA, 1);     % Nro de lineas en el sistema
ns = size(SHUNTDATA, 1);    % Nro de shunts en el sistema
Ybus = zeros(n,n);          % Se inicializa como una matriz de ceros

g = zeros(n, n);
b = zeros(n, n);

% Se agregan las impedancias de las lineas
for i = 1:nl 
    from = LINEDATA(i, 1);                      % Barra de inicio de la linea
    to = LINEDATA(i, 2);                        % Barra de fin de la linea
    zl = LINEDATA(i, 3) + 1i*LINEDATA(i, 4);    % Impedancia de la linea
    a = LINEDATA(i, 6);                         % Tap del transformador
    
    Ybus(from, from) = Ybus(from, from) + (1/zl)/(a^2);
    Ybus(to, to) = Ybus(to, to) + (1/zl)/(a^2);
    Ybus(from, to) = Ybus(from, to) - (1/zl)/a;
    Ybus(to, from) = Ybus(to, from) - (1/zl)/a;
    
    % Se agregan los shunts de las lineas
    Zshunt = 1i*LINEDATA(i, 5);
    Ybus(from, from) = Ybus(from, from) + Zshunt;
    Ybus(to, to) = Ybus(to, to) + Zshunt;
    
    g(from, to) = real(Zshunt);
    g(to, from) = real(Zshunt);
    
    b(from, to) = imag(Zshunt);
    b(to, from) = imag(Zshunt);
end

% Se agregan las impedancias de los shunts
if isempty(SHUNTDATA) == 0
    for i = 1:ns
        barra = SHUNTDATA(i, 1);                % barra donde esta conectado el shunt
        zs = 1i*SHUNTDATA(i, 2);                % impedancia del shunt
        Ybus(barra, barra) = Ybus(barra, barra) + 1/zs;
    end
end

G = real(Ybus);
B = imag(Ybus);



