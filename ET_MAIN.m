%%  Programa generico para estabilidad transitoria
clc, clear, close all;

%%  Los datos del archivo Excel se especifican de la siguiente manera
    % Hoja 1: Datos de barras
    % Hoja 2: Datos de ramas prefalla
    % Hoja 3: Datos de ramas falla
    % Hoja 4: Datos de ramas postfalla
    % Hoja 5: Generadores

% DATAFILE = 'BUSDATA.xlsx';
DATAFILE = 'Datos9barras_confalla.xlsx';
% DATAFILE = 'BUSDATA_fallabarra.xlsx';
% DATAFILE = 'BUSDATA_fallalong.xlsx';

ShowUnits = 0;
Vb = 115;
Sb = 100;

tic

[BUSDATA, LINEDATA_PRE, LINEDATA_FALLA, LINEDATA_POST, GENDATA, SIMULATIONDATA, FALLADATA] = ET_LoadData(DATAFILE);

%%  No tocar
ti = 0;
f = SIMULATIONDATA(1);
tp = SIMULATIONDATA(2)/f;
td = SIMULATIONDATA(3)/f;
tf = SIMULATIONDATA(4)/f;
dt = SIMULATIONDATA(5)/f;

tvec = [ti tp td tf dt];

we = 2*pi*f;
%%

n = size(BUSDATA, 1);                   % el numero de filas en el archivo excel es igual al numero de barras
nl_pre = size(LINEDATA_PRE, 1);         % el numero de filas en el archivo excel es igual al numero de ramas
nl_falla = size(LINEDATA_FALLA, 1);     % el numero de filas en el archivo excel es igual al numero de ramas
nl_post = size(LINEDATA_POST, 1);       % el numero de filas en el archivo excel es igual al numero de ramas
ng = size(GENDATA, 1);                  % el numero de filas en el archivo excel es igual al numero de generadores

%% Aqui se cargaran en vectores los parametros de las maquinas como lo son resistencias, reactancias, constantes de
% tiempo y constantes de inercia
for i = 1:ng
    Ra(i, 1) = GENDATA(i, 2);
    Xd(i, 1) = GENDATA(i, 3);
    Xdp(i, 1) = GENDATA(i, 4);
    Xdpp(i, 1) = GENDATA(i, 5);
    Xq(i, 1) = GENDATA(i, 6);
    Xqp(i, 1) = GENDATA(i, 7);
    Xqpp(i, 1) = GENDATA(i, 8);
    Td0p(i, 1) = GENDATA(i, 9);
    if(Td0p(i, 1) >= 1e10)
        Td0p(i, 1) = Inf;
    end
    Tq0p(i, 1) = GENDATA(i, 10);
    if(Tq0p(i, 1) >= 1e10)
        Tq0p(i, 1) = Inf;
    end
    Td0pp(i, 1) = GENDATA(i, 11);
    if(Td0pp(i, 1) >= 1e10)
        Td0pp(i, 1) = Inf;
    end
    Tq0pp(i, 1) = GENDATA(i, 12);
    if(Tq0pp(i, 1) >= 1e10)
        Tq0pp(i, 1) = Inf;
    end
    H(i, 1) = GENDATA(i, 13);
    
    Kv(i) = GENDATA(i, 14);
    Ki(i) = GENDATA(i, 15);
    Kt(i) = GENDATA(i, 16);
      
    R(i, 1) = GENDATA(i, 17);
    
    Tv(i, 1) = GENDATA(i, 18);
    if(Tv(i, 1) >= 1e10)
        Tv(i, 1) = Inf;
    end
    
    Tt(i, 1) = GENDATA(i, 19);
    if(Tt(i, 1) >= 1e10)
        Tt(i, 1) = Inf;
    end
end

%%  Formacion de la Ybus para el FDC
[Ybus_pre, G_pre, B_pre, g_pre, b_pre] = ET_Ybus(BUSDATA, LINEDATA_PRE, n, [], [], []);
%%  Ejecucion del FDC
[V, theta, Pgen, Qgen, Pneta, Qneta, Sshunt, Pflow, Pflow_bus, ...
Qflow, Qflow_bus, Ploss, Qloss, Pload, Qload] = ET_FDC(BUSDATA, LINEDATA_PRE, G_pre, B_pre, g_pre, b_pre);

ET_PrintFDC;

%%  Una vez que se realizo el FDC, se pueden modelar las cargas como impedancia y se agregan a la Ybus
[Ybusc_pre, Gc_pre, Bc_pre, gc_pre, bc_pre] = ET_Ybus(BUSDATA, LINEDATA_PRE, n, V, theta, []);
[Ybusc_falla, Gcf, Bcf, gcf, bcf] = ET_Ybus(BUSDATA, LINEDATA_FALLA, n, V, theta, FALLADATA);
[Ybusc_post, Gcpt, Bcpt, gcpt, bcpt] = ET_Ybus(BUSDATA, LINEDATA_POST, n, V, theta, []);

%%  Estos son los equivalentes de Kron sin impedancias de generadores. En este equivalente se eliminan barras PQ
[YKron_pre, YaR_pre, YbR_pre, YcR_pre, YdR_pre] = ET_Kron(Ybusc_pre, ng);                                % Se obtiene la matriz equivalente de Kron para el sistema pre-falla
[YKron_falla, YaRf, YbRf, YcRf, YdRf] = ET_Kron(Ybusc_falla, ng);                                % Se obtiene la matriz equivalente de Kron para el sistema pre-falla
[YKron_post, YaRpt, YbRpt, YcRpt, YdRpt] = ET_Kron(Ybusc_post, ng);                                % Se obtiene la matriz equivalente de Kron para el sistema pre-falla

%%  Matrices RM
YKrm_pre = ET_YRM(YKron_pre);
YKrm_falla = ET_YRM(YKron_falla);
YKrm_post = ET_YRM(YKron_post);

%%  Calculos de los VRM e IRM
[Vrm_pre, Irm_pre] = ET_VRMIRM(V, theta, YKrm_pre, BUSDATA);

%%  Calculos de las tensiones internas ficticas, a fin de obtener los angulos delta para cada generador
[Ef_pre, d0_pre, Vt_pre, It_pre] = ET_EFICT(Vrm_pre, Irm_pre, Ra, Xq, Xd, ng);

%%  Calculo de la potencia mecanica y electrica de los generadores, pre-falla
Pmgap0 = ET_PMEC(Vrm_pre, Irm_pre, Ra, ng);

%%  Calculo de la matriz T de Park para calcular la Eexc en el momento pre falla
Mp = ET_MATRIZM(Ra, Xqp, Xdp, ng);
Mpp = ET_MATRIZM(Ra, Xqpp, Xdpp, ng);
T = ET_TPARK(d0_pre, ng);
A = inv(T)*inv(Mp)*T;
App = inv(T)*inv(Mpp)*T;

%%  Calculo de tensiones y corrientes qd
[Vqdp, Iqd, Vqp, Vdp, Iq0, Id0] = ET_VIQD(Vrm_pre, Irm_pre, T, ng);
[Eq0, Ed0, Eqp0, Edp0, Eqpp0, Edpp0, Eexc] = ET_EQDEXC(Vqp, Vdp, Iq0, Id0, Ra, Xd, Xdp, Xdpp, Xq, Xqp, Xqpp, ng);
%%  %%%%%%%%%%%%%%%%%%%%%%%%%% Se empieza con las integraciones %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Se integrara por separado: prefalla, falla y postfalla

%%  Caso pre-falla
%   Se definen los valores iniciales para caso pre-falla
w0 = zeros(ng, 1);
d0 = d0_pre;
Pegap0 = Pmgap0;
Xv0 = Pegap0;
Pc0 = Pegap0;


[w_pre, d_pre, Eqp_pre, Eqpp_pre, Edp_pre, Edpp_pre, Pmgap_pre, Xv_pre, Pc_pre, Pegap_pre, Vt_pre, theta_pre, Iq_pre, Id_pre, ...
 w_falla, d_falla, Eqp_falla, Eqpp_falla, Edp_falla, Edpp_falla, Pmgap_falla, Xv_falla, Pc_falla, Pegap_falla, Vt_falla, theta_falla, Iq_falla, Id_falla, ...
 w_post, d_post, Eqp_post, Eqpp_post, Edp_post, Edpp_post, Pmgap_post, Xv_post, Pc_post, Pegap_post, Vt_post, theta_post, Iq_post, Id_post, Pegappostest] = ET_Integracion2(w0, d0, Pegap0, Iq0, Id0, Eqp0, Eqpp0, Edp0, Edpp0, Pmgap0, Xv0, Pc0, Ra, Xq, Xqp, Xqpp, Xd, Xdp, Xdpp, Tq0p, Tq0pp, Td0p, Td0pp, Tt, Kv, Tv, Kt, Ki, R, YKrm_pre, YKrm_falla, YKrm_post, Ybusc_pre, Ybusc_falla, Ybusc_post, Mp, Mpp, Eexc, we, H, ng, n, tvec);

% [w_pre, d_pre, Pegap_pre, Eqp_pre, Eqpp_pre, Edp_pre, Edpp_pre, Vt_pre, theta_pre, Iq_pre, Id_pre, Pmgap_pre, Xv_pre, Pc_pre] = ET_Integracion(w0, d0, Pegap0, Iq0, Id0, Eqp0, Eqpp0, Edp0, Edpp0, Pmgap0, Xv0, Pc0, Ra, Xq, Xqp, Xqpp, Xd, Xdp, Xdpp, Tq0p, Tq0pp, Td0p, Td0pp, Tt, Kv, Tv, Kt, Ki, R, YKrm_pre, Ybusc_pre, Mp, Mpp, Eexc, we, H, ng, n, [ti dt tp]);
% 
% %%  Caso falla
% for i = 1:ng
%     w0(i) = w_pre(i, length(w_pre));
%     d0(i) = d_pre(i, length(d_pre));
%     Pegap0(i) = Pegap_pre(i, length(Pegap_pre));
%     Eqp0(i) = Eqp_pre(i, length(Eqp_pre));
%     Edp0(i) = Edp_pre(i, length(Edp_pre));
%     Eqpp0(i) = Eqpp_pre(i, length(Eqpp_pre));
%     Edpp0(i) = Edpp_pre(i, length(Edpp_pre));
%     Pmgap0(i) = Pmgap_pre(i, length(Pmgap_pre));
%     Xv0(i) = Xv_pre(i, length(Xv_pre));
%     Pc0(i) = Pc_pre(i, length(Pc_pre));
% end
% 
% [w_falla, d_falla, Pegap_falla, Eqp_falla, Eqpp_falla, Edp_falla, Edpp_falla, Vt_falla, theta_falla, Iq_falla, Id_falla, Pmgap_falla, Xv_falla, Pc_falla] = ET_Integracion(w0, d0, Pegap0, Iq_pre, Id_pre, Eqp0, Eqpp0, Edp0, Edpp0, Pmgap0, Xv0, Pc0, Ra, Xq, Xqp, Xqpp, Xd, Xdp, Xdpp, Tq0p, Tq0pp, Td0p, Td0pp, Tt, Kv, Tv, Kt, Ki, R, YKrm_falla, Ybusc_falla, Mp, Mpp, Eexc, we, H, ng, n, [tp dt td]);
% 
% %% Caso post-falla
% for i = 1:ng
%     w0(i) = w_falla(i, length(w_falla));
%     d0(i) = d_falla(i, length(d_falla));
%     Pegap0(i) = Pegap_falla(i, length(Pegap_falla));
%     Eqp0(i) = Eqp_falla(i, length(Eqp_falla));
%     Edp0(i) = Edp_falla(i, length(Edp_falla));
%     Eqpp0(i) = Eqpp_falla(i, length(Eqpp_falla));
%     Edpp0(i) = Edpp_falla(i, length(Edpp_falla));
%     Pmgap0(i) = Pmgap_falla(i, length(Pmgap_falla));
%     Xv0(i) = Xv_falla(i, length(Xv_falla));
%     Pc0(i) = Pc_falla(i, length(Pc_falla));
% end
% [w_post, d_post, Pegap_post, Eqp_post, Eqpp_post, Edp_post, Edpp_post, Vt_post, theta_post, Iq_post, Id_post, Pmgap_post, Xv_post, Pc_post] = ET_Integracion(w0, d0, Pegap0, Iq_falla, Id_falla, Eqp0, Eqpp0, Edp0, Edpp0, Pmgap0, Xv0, Pc0, Ra, Xq, Xqp, Xqpp, Xd, Xdp, Xdpp, Tq0p, Tq0pp, Td0p, Td0pp, Tt, Kv, Tv, Kt, Ki, R, YKrm_post, Ybusc_post, Mp, Mpp, Eexc, we, H, ng, n, [td dt tf]);

% %% Velocidades en las barras
% for i = 1:n
%     vel_pre(i, 1:size(theta_pre, 2)-1) = diff(theta_pre(i, 1:size(theta_pre, 2)))./diff(ti:dt:tp);
% end
% 
% for i = 1:n
%     vel_falla(i, 1:size(theta_falla, 2)-1) = diff(theta_falla(i, 1:size(theta_falla, 2)))./dt;
% end
% 
% for i = 1:n
%     vel_post(i, 1:size(theta_post, 2)-1) = diff(theta_post(i, 1:size(theta_post, 2)))./dt;
% end
% 
% % frecuencias
% f_pre = vel_pre./(2*pi) + f;
% f_falla = vel_falla./(2*pi) + f;
% f_post = vel_post./(2*pi) + f;
% 
% %% Droops
% R = R./f;
% deltaP = 70/Sb; % negativo = bote de carga (deltaP)
% beta = sum(1./R);
% df = -(1/beta)*deltaP;
% for i = 1:ng
%     dPm(i) = -(1/R(i))*df*Sb;
% end

toc
%% Plots
ET_Plot(w_pre, d_pre, Pegap_pre, Eqp_pre, Eqpp_pre, Edp_pre, Edpp_pre, Vt_pre, theta_pre, Pmgap_pre, Xv_pre, Pc_pre, w_falla, d_falla, Pegap_falla, Eqp_falla, Eqpp_falla, Edp_falla, Edpp_falla, Vt_falla, theta_falla, Pmgap_falla, Xv_falla, Pc_falla, w_post, d_post, Pegap_post, Eqp_post, Eqpp_post, Edp_post, Edpp_post, Vt_post, theta_post, Pmgap_post, Xv_post, Pc_post, ng, n, [ti tp td tf dt]);
