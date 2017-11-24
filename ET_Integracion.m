%% Aqui se haran las iteraciones y se resolvera la ecuacion de oscilacion

function [w_pre, d_pre, Eqp_pre, Eqpp_pre, Edp_pre, Edpp_pre, Pmgap_pre, Xv_pre, Pc_pre, Pegap_pre, Vt_pre, theta_pre, Iq_pre, Id_pre, ...
          w_falla, d_falla, Eqp_falla, Eqpp_falla, Edp_falla, Edpp_falla, Pmgap_falla, Xv_falla, Pc_falla, Pegap_falla, Vt_falla, theta_falla, Iq_falla, Id_falla, ...
          w_post, d_post, Eqp_post, Eqpp_post, Edp_post, Edpp_post, Pmgap_post, Xv_post, Pc_post, Pegap_post, Vt_post, theta_post, Iq_post, Id_post, Pegapposttest] = ET_Integracion(w0, d0, Pegap0, Iq0, Id0, Eqp0, Eqpp0, Edp0, Edpp0, Pmgap0, Xv0, Pc0, Ra, Xq, Xqp, Xqpp, Xd, Xdp, Xdpp, Tq0p, Tq0pp, Td0p, Td0pp, Tt, Kv, Tv, Kt, Ki, R, YKrm_pre, YKrm_falla, YKrm_post, Ybus_pre, Ybus_falla, Ybus_post, Mp, Mpp, Eexc, we, H, AGC, JAULA, ng, n, tvec)

    ti = tvec(1);
    tp = tvec(2);
    td = tvec(3);
    tf = tvec(4);
    dt = tvec(5);
    
    size_pre = length(ti:dt:tp);
    size_falla = length(tp:dt:td);
    size_post = length(td:dt:tf);

    options = optimset('Display', 'off');
    
    %%  Se definen las matrices para la integracion pre-falla
    d_pre = zeros(ng, size_pre);
    w_pre = zeros(ng, size_pre);
    Pegap_pre = zeros(ng, size_pre);
    Eqp_pre = zeros(ng, size_pre);
    Edp_pre = zeros(ng, size_pre);
    Eqpp_pre = zeros(ng, size_pre);
    Edpp_pre = zeros(ng, size_pre);
    Vt_pre = zeros(ng, size_pre);
    theta_pre = zeros(ng, size_pre);
    Iq_pre = zeros(ng, size_pre);
    Id_pre = zeros(ng, size_pre);
    Pmgap_pre = zeros(ng, size_pre);
    Xv_pre = zeros(ng, size_pre);
    Pc_pre = zeros(ng, size_pre);
    
    %% Se definen las matrices para la integracion falla
    d_falla = zeros(ng, size_falla);
    w_falla = zeros(ng, size_falla);
    Pegap_falla = zeros(ng, size_falla);
    Eqp_falla = zeros(ng, size_falla);
    Edp_falla = zeros(ng, size_falla);
    Eqpp_falla = zeros(ng, size_falla);
    Edpp_falla = zeros(ng, size_falla);
    Vt_falla = zeros(ng, size_falla);
    theta_falla = zeros(ng, size_falla);
    Iq_falla = zeros(ng, size_falla);
    Id_falla = zeros(ng, size_falla);
    Pmgap_falla = zeros(ng, size_falla);
    Xv_falla = zeros(ng, size_falla);
    Pc_falla = zeros(ng, size_falla);
    
    %%  Se definen las matrices para la integracion post-falla
    d_post = zeros(ng, size_post);
    w_post = zeros(ng, size_post);
    Pegap_post = zeros(ng, size_post);
    Eqp_post = zeros(ng, size_post);
    Edp_post = zeros(ng, size_post);
    Eqpp_post = zeros(ng, size_post);
    Edpp_post = zeros(ng, size_post);
    Vt_post = zeros(ng, size_post);
    theta_post = zeros(ng, size_post);
    Iq_post = zeros(ng, size_post);
    Id_post = zeros(ng, size_post);
    Pmgap_post = zeros(ng, size_post);
    Xv_post = zeros(ng, size_post);
    Pc_post = zeros(ng, size_post);
    
    %%  Estas variables se utilizara para los tres casos de integracion
    wv = zeros(ng, 1);
    dv = zeros(ng, 1);
    Pegapv = zeros(ng, 1);
    Eqpv = zeros(ng, 1);
    Edpv = zeros(ng, 1);
    Eqppv = zeros(ng, 1);
    Edppv = zeros(ng, 1);
    Iqv = zeros(ng, 1);
    Idv = zeros(ng, 1);
    Pmgapv = zeros(ng, 1);
    Xvv = zeros(ng, 1);
    Pcv = zeros(ng, 1);
    
    Eqdpn = zeros(2*ng, 1);
    Eqdppn = zeros(2*ng, 1);
    Eqd = zeros(2*ng, 1);
    
    Vtn = zeros(ng, 1);
    Itn = zeros(ng, 1);
    
    vv = 4;
    vv = vv + sum(JAULA.*2);
    vv = vv + sum(AGC.*3);

    x0 = zeros(1, vv);
    
    %%  Integracion pre-falla
    for i = 1:ng
        wv(i) = w0(i);
        dv(i) = d0(i);
        Eqpv(i) = Eqp0(i);
        Edpv(i) = Edp0(i);
        Eqppv(i) = Eqpp0(i);
        Edppv(i) = Edpp0(i);
        Pmgapv(i) = Pmgap0(i);
        Xvv(i) = Xv0(i);
        Pcv(i) = Pc0(i);
        
        Pegapv(i) = Pegap0(i);
        Iqv(i) = Iq0(i);
        Idv(i) = Id0(i);
    end
    
    k = 1;
    for t = ti:dt:tp
        if(k > 1)       % Ya no se esta en la primera integracion, por tanto los valores viejos seran los de la integracion anterior
            v = 1;
            for i = 1:ng
                x0(v) = w_pre(i, k - 1);
                x0(v+1) = d_pre(i, k - 1);
                x0(v+2) = Eqp_pre(i, k - 1);
                x0(v+3) = Edp_pre(i, k - 1);
                v = v + 4;
                if(JAULA(i) == 1)
                    x0(v) = Eqpp_pre(i, k - 1);
                    x0(v+1) = Edpp_pre(i, k - 1);
                    v = v + 2;
                end
                if(AGC(i) == 1)
                    x0(v) = Pmgap_pre(i, k - 1);
                    x0(v+1) = Xv_pre(i, k - 1);
                    x0(v+2) = Pc_pre(i, k - 1);
                    v = v + 3;
                end
            end
        else
            % Se va a realizar la primera integracion, por tanto los valores viejos seran los valores calculados como iniciales
            
            v = 1;
            for i = 1:ng
                x0(v) = wv(i);                  % w0
                x0(v+1) = dv(i);                % d0
                x0(v+2) = Eqpv(i);              % Eqp0
                x0(v+3) = Edpv(i);              % Edp0
                v = v + 4;
                if(JAULA(i) == 1)
                    x0(v) = Eqppv(i);             % Eqpp0
                    x0(v+1) = Edppv(i);             % Edpp0
                    v = v + 2;
                end
                if(AGC(i) == 1)
                    x0(v) = Pmgapv(i);            % Pmgap
                    x0(v+1) = Xvv(i);               % Xv
                    x0(v+2) = Pcv(i);               % Pc
                    v = v + 3;
                end
            end
        end
        exitflag = 0;
        iterations = 0;
        while(exitflag ~= 1)
            [x, ~, exitflag, ~] = fsolve(@(x)ET_EcuacionesFSOLVE(x, YKrm_pre, wv, dv, Eqpv, Eqppv, Edpv, Edppv, Pmgapv, Xvv, Pcv, Pegapv, Iqv, Idv, Mp, Mpp, Ra, Xq, Xqp, Xqpp, Xd, Xdp, Xdpp, Tq0p, Tq0pp, Td0p, Td0pp, Tt, Kv, Tv, Kt, Ki, R, Eexc, H, JAULA, AGC, we, ng, dt), x0, options);

            iterations = iterations + 1;
            if(iterations == 100)
                fprintf('FSOLVE no convergio en 100 iteraciones. Caso: pre-falla\n');
                pause;
            end
        end
        
        %%  Se guardan los nuevos valores
        v = 1;
        for i = 1:ng
            w_pre(i, k) = x(v);
            d_pre(i, k) = x(v+1);
            Eqp_pre(i, k) = x(v+2);
            Edp_pre(i, k) = x(v+3);
            v = v + 4;
            if(JAULA(i) == 1)
                Eqpp_pre(i, k) = x(v);
                Edpp_pre(i, k) = x(v+1);
                v = v + 2;
            end
            if(AGC(i) == 1)
                Pmgap_pre(i, k) = x(v);
                Xv_pre(i, k) = x(v+1);
                Pc_pre(i, k) = x(v+2);
                v = v + 3;
            else
                Pmgap_pre(i, k) = Pmgap0(i);
            end
            
            Eqdpn(2*i-1, 1) = Eqp_pre(i, k);
            Eqdpn(2*i, 1) = Edp_pre(i, k);
            Eqdppn(2*i-1, 1) = Eqpp_pre(i, k);
            Eqdppn(2*i, 1) = Edpp_pre(i, k);
            
            Eqd(2*i-1, 1) = Eqp_pre(i, k);
            Eqd(2*i, 1) = Edp_pre(i, k);
            if(JAULA(i) == 1)
                Eqd(2*i-1, 1) = Eqpp_pre(i, k);
                Eqd(2*i, 1) = Edpp_pre(i, k);
            end
        end
        
        Tn = ET_TPARK(d_pre(1:ng, k), ng);
        An = inv(Tn)*inv(Mp)*Tn;
%         Appn = inv(Tn)*inv(Mpp)*Tn;

%         Ermn = inv(Tn)*Eqdpn;
        Ermn = Tn\Eqd;
        Vrmn = (YKrm_pre + An)\(An*Ermn);
        Irmn = YKrm_pre*Vrmn;
        
        Iqdn = Tn*Irmn;
        
        Ign = zeros(n, 1);
        for i = 1:ng
            Iq_pre(i, k) = Iqdn(2*i-1);
            Id_pre(i, k) = Iqdn(2*i);
            
            Vtn(i) = Vrmn(2*i-1) + 1i*Vrmn(2*i);
            Itn(i) = Irmn(2*i-1) + 1i*Irmn(2*i);
            
            Ign(i) = Itn(i);
            Pegap_pre(i, k) = real(Vtn(i)*conj(Itn(i))) + Ra(i)*abs(Itn(i))^2;
        end
        
        Vbn = Ybus_pre\Ign;
        
        for i = 1:n
            Vt_pre(i, k) = abs(Vbn(i));
            theta_pre(i, k) = angle(Vbn(i));
        end
        
        for i = 1:ng
            wv(i) = w_pre(i, k);
            dv(i) = d_pre(i, k);
            Eqpv(i) = Eqp_pre(i, k);
            Edpv(i) = Edp_pre(i, k);
            if(JAULA(i) == 1)
                Eqppv(i) = Eqpp_pre(i, k);
                Edppv(i) = Edpp_pre(i, k);
            end
            if(AGC(i) == 1)
                Pmgapv(i) = Pmgap_pre(i, k);
                Xvv(i) = Xv_pre(i, k);
                Pcv(i) = Pc_pre(i, k);
            end

            Pegapv(i) = Pegap_pre(i, k);
            Iqv(i) = Iq_pre(i, k);
            Idv(i) = Id_pre(i, k);
        end
        
        k = k + 1;
    end

    %% Integracion falla
    
    k = 1;
    for t = tp:dt:td
        if(k > 1)       % Ya no se esta en la primera integracion, por tanto los valores viejos seran los de la integracion anterior
            v = 1;
            for i = 1:ng
                x0(v) = w_falla(i, k - 1);
                x0(v+1) = d_falla(i, k - 1);
                x0(v+2) = Eqp_falla(i, k - 1);
                x0(v+3) = Edp_falla(i, k - 1);
                v = v + 4;
                if(JAULA(i) == 1)
                    x0(v) = Eqpp_falla(i, k - 1);
                    x0(v+1) = Edpp_falla(i, k - 1);
                    v = v + 2;
                end
                if(AGC(i) == 1)
                    x0(v) = Pmgap_falla(i, k - 1);
                    x0(v+1) = Xv_falla(i, k - 1);
                    x0(v+2) = Pc_falla(i, k - 1);
                    v = v + 3;
                end
            end
        else
            % Se va a realizar la primera integracion, por tanto los valores viejos seran los valores calculados como iniciales
            v = 1;
            for i = 1:ng
                x0(v) = wv(i);                  % w0
                x0(v+1) = dv(i);                % d0
                x0(v+2) = Eqpv(i);              % Eqp0
                x0(v+3) = Edpv(i);              % Edp0
                v = v + 4;
                if(JAULA(i) == 1)
                    x0(v) = Eqppv(i);             % Eqpp0
                    x0(v+1) = Edppv(i);             % Edpp0
                    v = v + 2;
                end
                if(AGC(i) == 1)
                    x0(v) = Pmgapv(i);            % Pmgap
                    x0(v+1) = Xvv(i);               % Xv
                    x0(v+2) = Pcv(i);               % Pc
                    v = v + 3;
                end
            end
        end
        
        exitflag = 0;
        iterations = 0;
        while(exitflag ~= 1)
            [x, ~, exitflag, ~] = fsolve(@(x)ET_EcuacionesFSOLVE(x, YKrm_falla, wv, dv, Eqpv, Eqppv, Edpv, Edppv, Pmgapv, Xvv, Pcv, Pegapv, Iqv, Idv, Mp, Mpp, Ra, Xq, Xqp, Xqpp, Xd, Xdp, Xdpp, Tq0p, Tq0pp, Td0p, Td0pp, Tt, Kv, Tv, Kt, Ki, R, Eexc, H, JAULA, AGC, we, ng, dt), x0, options);

            iterations = iterations + 1;
            if(iterations == 100)
                fprintf('FSOLVE no convergio en 100 iteraciones. Caso: falla\n');
                pause;
            end
        end
        
        %%  Se guardan los nuevos valores
        if(t == tp)
            for i = 1:ng
                % Se esta en t0+, por tanto los valores de las variables de estado deben permanecer igual a t0-
                % Estos corresponden a los ultimos calculados en
                % pre-falla
                w_falla(i, k) = w_pre(i, size_pre);
                d_falla(i, k) = d_pre(i, size_pre);
                Eqp_falla(i, k) = Eqp_pre(i, size_pre);
                Edp_falla(i, k) = Edp_pre(i, size_pre);
                if(JAULA(i) == 1)
                    Eqpp_falla(i, k) = Eqpp_pre(i, size_pre);
                    Edpp_falla(i, k) = Edpp_pre(i, size_pre);
                end
                if(AGC(i) == 1)
                    Pmgap_falla(i, k) = Pmgap_pre(i, size_pre);
                    Xv_falla(i, k) = Xv_pre(i, size_pre);
                    Pc_falla(i, k) = Pc_pre(i, size_pre);
                else
                    Pmgap_falla(i, k) = Pmgap0(i);
                end
            end
        else
            %   Si ya se esta en t0+ +dt entonces si se guarda el valor
            %   calculado
            v = 1;
            for i = 1:ng
                w_falla(i, k) = x(v);
                d_falla(i, k) = x(v+1);
                Eqp_falla(i, k) = x(v+2);
                Edp_falla(i, k) = x(v+3);
                v = v + 4;
                if(JAULA(i) == 1)
                    Eqpp_falla(i, k) = x(v);
                    Edpp_falla(i, k) = x(v+1);
                    v = v + 2;
                end
                if(AGC(i) == 1)
                    Pmgap_falla(i, k) = x(v);
                    Xv_falla(i, k) = x(v+1);
                    Pc_falla(i, k) = x(v+2);
                    v = v + 3;
                else
                    Pmgap_falla(i, k) = Pmgap0(i);
                end
                
                Eqdpn(2*i-1, 1) = Eqp_falla(i, k);
                Eqdpn(2*i, 1) = Edp_falla(i, k);
                Eqdppn(2*i-1, 1) = Eqpp_falla(i, k);
                Eqdppn(2*i, 1) = Edpp_falla(i, k);
                
                Eqd(2*i-1, 1) = Eqp_falla(i, k);
                Eqd(2*i, 1) = Edp_falla(i, k);
                if(JAULA(i) == 1)
                    Eqd(2*i-1, 1) = Eqpp_falla(i, k);
                    Eqd(2*i, 1) = Edpp_falla(i, k);
                end
            end
        end
        
        Tn = ET_TPARK(d_falla(1:ng, k), ng);
        An = inv(Tn)*inv(Mp)*Tn;
%         Appn = inv(Tn)*inv(Mpp)*Tn;

%         Ermn = inv(Tn)*Eqdpn;
        Ermn = Tn\Eqd;
        Vrmn = (YKrm_falla + An)\(An*Ermn);
        Irmn = YKrm_falla*Vrmn;
        
        Iqdn = Tn*Irmn;
        
        Ign = zeros(n, 1);
        for i = 1:ng
            Iq_falla(i, k) = Iqdn(2*i-1);
            Id_falla(i, k) = Iqdn(2*i);
            
            Vtn(i) = Vrmn(2*i-1) + 1i*Vrmn(2*i);
            Itn(i) = Irmn(2*i-1) + 1i*Irmn(2*i);
            
            Ign(i) = Itn(i);
            
            Pegap_falla(i, k) = real(Vtn(i)*conj(Itn(i))) + Ra(i)*abs(Itn(i))^2;
        end
        
        Vbn = Ybus_falla\Ign;
        
        for i = 1:n
            Vt_falla(i, k) = abs(Vbn(i));
            theta_falla(i, k) = angle(Vbn(i));
        end
        
        %% Se almacenan valores viejos
        for i = 1:ng
            wv(i) = w_falla(i, k);
            dv(i) = d_falla(i, k);
            Eqpv(i) = Eqp_falla(i, k);
            Edpv(i) = Edp_falla(i, k);
            if(JAULA(i) == 1)
                Eqppv(i) = Eqpp_falla(i, k);
                Edppv(i) = Edpp_falla(i, k);
            end
            if(AGC(i) == 1)
                Pmgapv(i) = Pmgap_falla(i, k);
                Xvv(i) = Xv_falla(i, k);
                Pcv(i) = Pc_falla(i, k);
            end

            Pegapv(i) = Pegap_falla(i, k);
            Iqv(i) = Iq_falla(i, k);
            Idv(i) = Id_falla(i, k);
        end
        
        k = k + 1;
    end
    
    %%  Integracion post-falla
    %   Los valores viejos de la primera integracion seran los ultimos
    %   calculados en el caso falla
    
    k = 1;
    for t = td:dt:tf
        if(k > 1)       % Ya no se esta en la primera integracion, por tanto los valores viejos seran los de la integracion anterior
            v = 1;
            for i = 1:ng
                x0(v) = w_post(i, k - 1);
                x0(v+1) = d_post(i, k - 1);
                x0(v+2) = Eqp_post(i, k - 1);
                x0(v+3) = Edp_post(i, k - 1);
                v = v + 4;
                if(JAULA(i) == 1)
                    x0(v) = Eqpp_post(i, k - 1);
                    x0(v+1) = Edpp_post(i, k - 1);
                    v = v + 2;
                end
                if(AGC(i) == 1)
                    x0(v) = Pmgap_post(i, k - 1);
                    x0(v+1) = Xv_post(i, k - 1);
                    x0(v+2) = Pc_post(i, k - 1);
                    v = v + 3;
                end
            end
        else
            % Se va a realizar la primera integracion, por tanto los valores viejos seran los valores calculados como iniciales
            v = 1;
            for i = 1:ng
                x0(v) = wv(i);                  % w0
                x0(v+1) = dv(i);                % d0
                x0(v+2) = Eqpv(i);              % Eqp0
                x0(v+3) = Edpv(i);              % Edp0
                v = v + 4;
                if(JAULA(i) == 1)
                    x0(v) = Eqppv(i);             % Eqpp0
                    x0(v+1) = Edppv(i);             % Edpp0
                    v = v + 2;
                end
                if(AGC(i) == 1)
                    x0(v) = Pmgapv(i);            % Pmgap
                    x0(v+1) = Xvv(i);               % Xv
                    x0(v+2) = Pcv(i);               % Pc
                    v = v + 3;
                end
            end
        end
        
        exitflag = 0;
        iterations = 0;
        while(exitflag ~= 1)
            [x, ~, exitflag, ~] = fsolve(@(x)ET_EcuacionesFSOLVE(x, YKrm_post, wv, dv, Eqpv, Eqppv, Edpv, Edppv, Pmgapv, Xvv, Pcv, Pegapv, Iqv, Idv, Mp, Mpp, Ra, Xq, Xqp, Xqpp, Xd, Xdp, Xdpp, Tq0p, Tq0pp, Td0p, Td0pp, Tt, Kv, Tv, Kt, Ki, R, Eexc, H, JAULA, AGC, we, ng, dt), x0, options);

            iterations = iterations + 1;
            if(iterations == 100)
                fprintf('FSOLVE no convergio en 100 iteraciones. Caso: post-falla\n');
                pause;
            end
        end
        
        %%  Se guardan los nuevos valores
        if(t == td)
            for i = 1:ng
                % Se esta en t0+, por tanto los valores de las variables de estado deben permanecer igual a t0-
                % Estos corresponden a los ultimos calculados en
                % pre-falla
                w_post(i, k) = w_falla(i, size_falla);
                d_post(i, k) = d_falla(i, size_falla);
                Eqp_post(i, k) = Eqp_falla(i, size_falla);
                Edp_post(i, k) = Edp_falla(i, size_falla);
                if(JAULA(i) == 1)
                    Eqpp_post(i, k) = Eqpp_falla(i, size_falla);
                    Edpp_post(i, k) = Edpp_falla(i, size_falla);
                end
                if(AGC(i) == 1)
                    Pmgap_post(i, k) = Pmgap_falla(i, size_falla);
                    Xv_post(i, k) = Xv_falla(i, size_falla);
                    Pc_post(i, k) = Pc_falla(i, size_falla);
                else
                    Pmgap_post(i, k) = Pmgap0(i);
                end
            end
        else
            %   Si ya se esta en t0+ +dt entonces si se guarda el valor
            %   calculado
            v = 1;
            for i = 1:ng
                w_post(i, k) = x(v);
                d_post(i, k) = x(v+1);
                Eqp_post(i, k) = x(v+2);
                Edp_post(i, k) = x(v+3);
                v = v + 4;
                if(JAULA(i) == 1)
                    Eqpp_post(i, k) = x(v);
                    Edpp_post(i, k) = x(v+1);
                    v = v + 2;
                end
                if(AGC(i) == 1)
                    Pmgap_post(i, k) = x(v);
                    Xv_post(i, k) = x(v+1);
                    Pc_post(i, k) = x(v+2);
                    v = v + 3;
                else
                    Pmgap_post(i, k) = Pmgap0(i);
                end
                
                Eqdpn(2*i-1, 1) = Eqp_post(i, k);
                Eqdpn(2*i, 1) = Edp_post(i, k);
                Eqdppn(2*i-1, 1) = Eqpp_post(i, k);
                Eqdppn(2*i, 1) = Edpp_post(i, k);
                
                Eqd(2*i-1, 1) = Eqp_post(i, k);
                Eqd(2*i, 1) = Edp_post(i, k);
                if(JAULA(i) == 1)
                    Eqd(2*i-1, 1) = Eqpp_post(i, k);
                    Eqd(2*i, 1) = Edpp_post(i, k);
                end
            end
        end
        
        Tn = ET_TPARK(d_post(1:ng, k), ng);
        An = inv(Tn)*inv(Mp)*Tn;
%         Appn = inv(Tn)*inv(Mpp)*Tn;

%         Ermn = inv(Tn)*Eqdpn;
        Ermn = Tn\Eqd;
        Vrmn = (YKrm_post + An)\(An*Ermn);
        Irmn = YKrm_post*Vrmn;
        
        Iqdn = Tn*Irmn;
        
        Ign = zeros(n, 1);
        for i = 1:ng
            Iq_post(i, k) = Iqdn(2*i-1);
            Id_post(i, k) = Iqdn(2*i);
            
            Vtn(i) = Vrmn(2*i-1) + 1i*Vrmn(2*i);
            Itn(i) = Irmn(2*i-1) + 1i*Irmn(2*i);
            
            Ign(i) = Itn(i);
            
            Pegap_post(i, k) = real(Vtn(i)*conj(Itn(i))) + Ra(i)*abs(Itn(i))^2;
        end
        
        Vbn = Ybus_post\Ign;
        
        for i = 1:n
            Vt_post(i, k) = abs(Vbn(i));
            theta_post(i, k) = angle(Vbn(i));
        end
        
        %% Se almacenan valores viejos
        for i = 1:ng
            wv(i) = w_post(i, k);
            dv(i) = d_post(i, k);
            Eqpv(i) = Eqp_post(i, k);
            Edpv(i) = Edp_post(i, k);
            if(JAULA(i) == 1)
                Eqppv(i) = Eqpp_post(i, k);
                Edppv(i) = Edpp_post(i, k);
            end
            if(AGC(i) == 1)
                Pmgapv(i) = Pmgap_post(i, k);
                Xvv(i) = Xv_post(i, k);
                Pcv(i) = Pc_post(i, k);
            end
            
            Pegapv(i) = Pegap_post(i, k);
            Iqv(i) = Iq_post(i, k);
            Idv(i) = Id_post(i, k);
        end

        k = k + 1;
    end
end
