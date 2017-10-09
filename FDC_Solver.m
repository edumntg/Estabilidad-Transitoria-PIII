%% Eduardo Montilva 12-10089
% Script para la solucion del flujo de carga, mediante fsolve

% global n Vabs delta G B bustype SLACK PV PQ REF Pneta Qneta Ki Ybus Pdesbalance

bustype = zeros(n, 1);                          % Usado para obtener el tipo de barra (SLACK/REF/PV/PQ)
Pneta = zeros(n, 1);                            % Potencia activa neta de cada barra
Qneta = zeros(n, 1);                            % Potencia reactiva neta de cada barra
delta = zeros(n, 1);                            % Angulo del voltaje de cada barra (especificado)
Vabs = zeros(n, 1);                             % Voltaje de cada barra (especificado)
Ki = zeros(n, 1);                               % Factor de distribucion de cada barra

Sgen = zeros(n, 1);                             
Sload = zeros(n, 1);
Vrect = zeros(n, 1);

Pgen = zeros(n, 1);
Qgen = zeros(n, 1);

Pconsig = zeros(n, 1);
Qconsig = zeros(n, 1);

Pflow = zeros(n, n);
Qflow = zeros(n, n);

Pflow_bus = zeros(n, 1);
Qflow_bus = zeros(n, 1);

Ploss = zeros(n,n);
Qloss = zeros(n,n);
Ploss_total = 0;
Qloss_total = 0;
Pload = zeros(n, 1);
Qload = zeros(n, 1);
Sshunt = zeros(n, 1);

X0 = zeros(2*n, 1);
for i = 1:n 
    
    %% Se crean variables que seran de utilidad durante el flujo de carga
    
    bustype(i) = BUSDATA(i, 2); %PV, PQ o Slack
    Vabs(i) = BUSDATA(i, 3);
    delta(i) = BUSDATA(i, 4);
    
    %% Calculamos la potencia neta por barra

    Pload(i) = -BUSDATA(i, 5);
    Qload(i) = -BUSDATA(i, 6);
    
    Pconsig(i) = BUSDATA(i, 7);
    Qconsig(i) = BUSDATA(i, 8);
    
    %% Factor de distribucion
    Ki(i) = BUSDATA(i, 9);
    
    Pneta(i) = Pconsig(i)-Pload(i);
    Qneta(i) = Qconsig(i)-Qload(i);
end

Pdesbalance = sum(Pconsig) - abs(sum(Pload));

V = Vabs;
d = delta;

for i = 1:n
    if bustype(i) == SLACK || bustype(i) == REF % incognitas: P y Q
        X0(2*i-1) = Pconsig(i);
        X0(2*i) = Qconsig(i);
    elseif bustype(i) == PV % incognitas: delta y Q
        X0(2*i-1) = d(i);
        X0(2*i) = Qconsig(i);
    elseif bustype(i) == PQ %incognitas: V y delta
        X0(2*i-1) = d(i);
        X0(2*i) = V(i);
    end
end

%% Ejecucion del fsolve (iteraciones)
options = optimset('Display','off');
X = fsolve(@(x)FDC_SolverFunc(x, SHUNTDATA, bustype, V, d, Ki, Pload, Qload, Pconsig, Qconsig, Ybus, g, b, Pdesbalance), X0, options);

%% Una vez terminadas las iteraciones, se obtienen las variables de salida y se recalculan potencias
for i = 1:n
    if bustype(i) == SLACK || bustype(i) == REF % incognitas: P y Q
        Pgen(i) = X(2*i-1);
        Qgen(i) = X(2*i);
    elseif bustype(i) == PV % incognitas: delta y Q
        d(i) = X(2*i-1);
        Qgen(i) = X(2*i);
    elseif bustype(i) == PQ %incognitas: V y delta
        d(i) = X(2*i-1);
        V(i) = X(2*i);
    end
end

%% Calculo de flujos en lineas y perdidas
for i = 1:n
    for k = 1:n
        if i ~= k
            Pflow(i,k) = (-G(i,k) + g(i,k))*V(i)^2 + V(i)*V(k)*(G(i,k)*cos(d(i) - d(k)) + B(i,k)*sin(d(i) - d(k)));
            Qflow(i,k) = (B(i,k) - b(i,k))*V(i)^2 + V(i)*V(k)*(-B(i,k)*cos(d(i) - d(k)) + G(i,k)*sin(d(i) - d(k)));
        end
    end
    Pflow_bus(i) = sum(Pflow(i, 1:size(Pflow, 2)));
    Qflow_bus(i) = sum(Qflow(i, 1:size(Qflow, 2)));
end

%% Calculo de las perdidas
for i = 1:n
    for k = 1:n
        %% Calculo de perdidas
        if i ~= k
            Ploss(i,k) = Pflow(i,k) + Pflow(k,i);
            Qloss(i,k) = Qflow(i,k) + Qflow(k,i);
            
            if k > i
                Ploss_total = Ploss_total + Ploss(i,k);
                Qloss_total = Qloss_total + Qloss(i,k);
            end
        end
    end
end

%% Calculos de la P neta
% for i = 1:n
%     Pneta(i) = 0;
%     Qneta(i) = 0;
%     for k = 1:n
%         Pneta(i) = Pneta(i) + V(i)*V(k)*(G(i,k)*cos(d(i) - d(k)) + B(i,k)*sin(d(i) - d(k)));
%         Qneta(i) = Qneta(i) + V(i)*V(k)*(-B(i,k)*cos(d(i) - d(k)) + G(i,k)*sin(d(i) - d(k)));
%     end
% end

if isempty(SHUNTDATA) ~= 1
    for i = 1:size(SHUNTDATA, 1)
        bus = SHUNTDATA(i, 1);
        z = 1i*SHUNTDATA(i, 2);
        Sshunt(bus) = conj(z)\V(bus)^2;
        Qneta(bus) = Qneta(bus) - imag(Sshunt(bus));
    end
end

for i = 1:n
    Pgen(i) = 0;
    Qgen(i) = 0;
    if bustype(i) ~= PQ
        Pgen(i) = abs(Pload(i)) + Pflow_bus(i);
        Qgen(i) = abs(Qload(i)) + Qflow_bus(i) + imag(Sshunt(i));
    end
end

for i = 1:n
    Pneta(i) = Pgen(i) - abs(Pload(i));
    Qneta(i) = Qgen(i) - abs(Qload(i)) - imag(Sshunt(i));
end

Pdesbalance = Pdesbalance + Ploss_total;

%% VARIABLES PARA GARANTIZAR EL BUEN FUNCIONAMIENTO DEL PROGRAMA
% La P de salida en cada barra debe ser igual a la P neta de la misma
fprintf('Diferencia entre Pneta y Psalida para cada barra: %s\n', mat2str(Pgen - abs(Pload) - Pflow_bus));
fprintf('Diferencia entre Qneta y Qsalida para cada barra: %s\n', mat2str(Qgen - imag(Sshunt) - abs(Qload) - Qflow_bus));
Pdesbalance_result = sum(Pgen) - abs(sum(Pload));
fprintf('Desbalance inicial en el sistema: %s\n', num2str(Pdesbalance*Sb^ShowUnits));
fprintf('Desbalance final en el sistema: %s\n\n', num2str(Pdesbalance_result*Sb^ShowUnits));