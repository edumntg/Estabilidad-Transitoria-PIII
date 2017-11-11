%% Aqui se haran las iteraciones y se resolvera la ecuacion de oscilacion

function [w, d, Pegap, Eqp, Eqpp, Edp, Edpp, Vt, theta, Iq, Id] = ET_Integracion(wo, do, Peo, Iqo, Ido, Ra, Xq, Xqp, Xqpp, Xd, Xdp, Xdpp, Eqpo, Eqppo, Edpo, Edppo, Tq0p, Tq0pp, Td0p, Td0pp, Pmgap, YKrm, Ybus, Mp, Mpp, Eexc, we, H, ng, n, tvec)

    ti = tvec(1);
    dt = tvec(2);
    tf = tvec(3);

    options = optimset('Display', 'off');

    d = zeros(ng, length(ti:dt:tf));
    w = zeros(ng, length(ti:dt:tf));
    Pegap = zeros(ng, length(ti:dt:tf));
    Eqp = zeros(ng, length(ti:dt:tf));
    Edp = zeros(ng, length(ti:dt:tf));
    Eqpp = zeros(ng, length(ti:dt:tf));
    Edpp = zeros(ng, length(ti:dt:tf));
    Vt = zeros(ng, length(ti:dt:tf));
    theta = zeros(ng, length(ti:dt:tf));
    Iq = zeros(ng, length(ti:dt:tf));
    Id = zeros(ng, length(ti:dt:tf));
    
    Pmgapv = Pmgap;
    
    for i = 1:ng
        wv(i) = wo(i);
        dv(i) = do(i);
        Pegapv(i) = Peo(i);
        Eqpv(i) = Eqpo(i);
        Edpv(i) = Edpo(i);
        Eqppv(i) = Eqppo(i);
        Edppv(i) = Edppo(i);
        Iqv(i) = Iqo(i);
        Idv(i) = Ido(i);
    end

    k = 1;
    %% Se empieza a iterar
    for t = ti:dt:tf 
        if(k > 1)
            v = 1;
            for i = 1:ng
                
                x0(v) = w(i, k - 1); % w0
                x0(v+1) = d(i, k - 1); % d0
                x0(v+2) = Eqp(i, k - 1); % Eqp0
                x0(v+3) = Edp(i, k - 1); % Edp0
                x0(v+4) = Eqpp(i, k - 1); % Eqpp0
                x0(v+5) = Edpp(i, k - 1); % Edpp0
                
                wv(i) = w(i, k - 1);
                dv(i) = d(i, k - 1);
                Pegapv(i) = Pegap(i, k - 1);
                Eqpv(i) = Eqp(i, k - 1);
                Edpv(i) = Edp(i, k - 1);
                Eqppv(i) = Eqpp(i, k - 1);
                Edppv(i) = Edpp(i, k - 1);
                Iqv(i) = Iq(i, k - 1);
                Idv(i) = Id(i, k - 1);
                
                v = v + 6;
            end
        else
            v = 1;
            for i = 1:ng
                x0(v) = wv(i); % w0
                x0(v+1) = dv(i); % d0
                x0(v+2) = Eqpv(i); % Eqp0
                x0(v+3) = Edpv(i); % Edp0
                x0(v+4) = Eqppv(i); % Eqpp0
                x0(v+5) = Edppv(i); % Edpp0
                v = v + 6;
            end
        end

        [x,fval,exitflag,output,jacobian] = fsolve(@(x)ET_EcuacionesFSOLVE(x, YKrm, wv, dv, Eqpv, Eqppv, Edpv, Edppv, Iqv, Idv, Mp, Mpp, Pegapv, Pmgapv, Ra, Xq, Xqp, Xqpp, Xd, Xdp, Xdpp, Tq0p, Tq0pp, Td0p, Td0pp, Eexc, H, we, ng, dt), x0, options);
        
        v = 1;
        for i = 1:ng
            w(i, k) = x(v);
            d(i, k) = x(v+1);
            Eqp(i, k) = x(v+2);
            Edp(i, k) = x(v+3);
            Eqpp(i, k) = x(v+4);
            Edpp(i, k) = x(v+5);
            v = v + 6;
        end
        
        Tn = ET_TPARK(d(1:ng, k), ng);
        An = inv(Tn)*inv(Mp)*Tn;
        Appn = inv(Tn)*inv(Mpp)*Tn;
        
        for i = 1:ng
            Eqdpn(2*i-1, 1) = Eqp(i, k);
            Eqdpn(2*i, 1) = Edp(i, k);
            Eqdppn(2*i-1, 1) = Eqpp(i, k);
            Eqdppn(2*i, 1) = Edpp(i, k);
        end

%         Ermn = inv(Tn)*Eqdpn;
        Ermnp = inv(Tn)*Eqdppn;
        Vrmn = inv(YKrm + Appn)*Appn*Ermnp;
        Irmn = YKrm*Vrmn;
        
        Iqdn = Tn*Irmn;
        
        for i = 1:ng
            Iq(i, k) = Iqdn(2*i-1);
            Id(i, k) = Iqdn(2*i);
        end
        
        for i = 1:ng
            Vtn(i) = Vrmn(2*i-1) + 1i*Vrmn(2*i);
            Itn(i) = Irmn(2*i-1) + 1i*Irmn(2*i);
            
            Vt(i, k) = abs(Vtn(i));
            theta(i, k) = angle(Vtn(i));
        end
        
        for i = 1:ng
            Pegap(i, k) = real(Vtn(i)*conj(Itn(i))) + Ra(i)*abs(Itn(i))^2;
        end
        
        %% vector de corrientes generadas
        Ign = zeros(n, 1);
        for i = 1:ng
            Ign(i) = Itn(i);
        end
        
        Vbn = inv(Ybus)*Ign;
        
        for i = 1:size(Ybus, 1);
            Vt(i, k) = abs(Vbn(i));
            theta(i, k) = angle(Vbn(i));
        end
        
        k = k + 1;
    end
end
