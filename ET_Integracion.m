%% Aqui se haran las iteraciones y se resolvera la ecuacion de oscilacion

function [w, d, Pe, Eqp, Edp, Vt] = ET_Integracion(wo, do, Peo, Iq, Id, Ra, Xq, Xqp, Xd, Xdp, Eqpo, Edpo, Tq0p, Td0p, Pmp, YKronRM, Mp, Eexc, w0, H, ng, tvec)

    ti = tvec(1);
    dt = tvec(2);
    tf = tvec(3);

    options = optimset('Display', 'off');

    d = zeros(ng, length(ti:dt:tf));
    w = zeros(ng, length(ti:dt:tf));
    Pe = zeros(ng, length(ti:dt:tf));
    Eqp = zeros(ng, length(ti:dt:tf));
    Edp = zeros(ng, length(ti:dt:tf));
    Vt = zeros(ng, length(ti:dt:tf));
    
    Iq = zeros(ng, length(ti:dt:tf));
    Id = zeros(ng, length(ti:dt:tf));
    
    Pmv = Pmp;
    
    Peij = zeros(ng, ng);
    
    for i = 1:ng
        wv(i) = wo(i);
        dv(i) = do(i);
        Pev(i) = Peo(i);
        Eqpv(i) = Eqpo(i);
        Edpv(i) = Edpo(i);
        Iqv(i) = Iq(i);
        Idv(i) = Id(i);
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
                
                wv(i) = w(i, k - 1);
                dv(i) = d(i, k - 1);
                Pev(i) = Pe(i, k - 1);
                Eqpv(i) = Eqp(i, k - 1);
                Edpv(i) = Edp(i, k - 1);
                Iqv(i) = Iq(i, k - 1);
                Idv(i) = Id(i, k - 1);
                
                v = v + 4;
            end
        else
            v = 1;
            for i = 1:ng
                x0(v) = w(i, k); % w0
                x0(v+1) = d(i, k); % d0
                x0(v+2) = Eqp(i, k); % Eqp0
                x0(v+3) = Edp(i, k); % Edp0
                v = v + 4;
            end
        end

        [x,fval,exitflag,output,jacobian] = fsolve(@(x)ET_EcuacionesFSOLVE(x, YKronRM, wv, dv, Eqpv, Edpv, Iqv, Idv, Mp, Pev, Pmv, Ra, Xq, Xqp, Xd, Xdp, Tq0p, Td0p, Eexc, H, w0, ng, dt), x0, options);
        
        v = 1;
        for i = 1:ng
            w(i, k) = x(v);
            d(i, k) = x(v+1);
            Eqp(i, k) = x(v+2);
            Edp(i, k) = x(v+3);
            v = v + 4;
        end
        
        Tn = ET_TPARK(d(1:ng, k), ng);
        An = inv(Tn).*Mp.*Tn;
        for i = 1:ng
            Eqdpn(2*i-1, 1) = Eqp(i, k);
            Eqdpn(2*i, 1) = Edp(i, k);
        end

        Ermn = inv(Tn)*Eqdpn;
        Vrmn = (YKronRM + An)*An*Ermn;
        Irmn = YKronRM*Vrmn;
        
        [Vqdn, Iqdn, Vqpn, Vdpn, Iqn, Idn] = ET_VIQD(Vrmn, Irmn, Tn, ng);
        
        for i = 1:ng
            Vtn(i) = Vrmn(2*i-1) + 1i*Vrmn(2*i);
            Itn(i) = Irmn(2*i-1) + 1i*Irmn(2*i);
            
            Vt(i, k) = Vtn(i);
        end
        
        for i = 1:ng
            Pe(i, k) = real(Vtn(i)*Itn(i)) + Ra(i)*abs(Itn(i))^2;
            Iq(i, k) = Iqn(i);
            Id(i, k) = Idn(i);
        end
        
        k = k + 1;
    end
end
