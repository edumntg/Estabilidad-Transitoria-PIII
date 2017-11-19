function F = ET_EcuacionesFSOLVE(x, YKrm, wv, dv, Eqpv, Eqppv, Edpv, Edppv, Iqv, Idv, Pegapv, Pmgapv, Xvv, Pcv, Mp, Mpp, Ra, Xq, Xqp, Xqpp, Xd, Xdp, Xdpp, Tq0p, Tq0pp, Td0p, Td0pp, Tt, Kv, Tv, Kt, Ki, R, Eexc, H, we, ng, dt)

    % Tendremos 4 ecuaciones por generador, en total 4*n_gen ecuaciones
    % Tambien tendremos 4 variables, las cuales seran wn, dn, eqp, edp (n de nuevo)
    Ka = 0.00;
    
    v = 1;
    for i = 1:ng
        wn(i, 1) = x(v);
        dn(i, 1) = x(v+1);
        Eqpn(i, 1) = x(v+2);
        Edpn(i, 1) = x(v+3);
        Eqppn(i, 1) = x(v+4);
        Edppn(i, 1) = x(v+5);
        Pmgapn(i, 1) = x(v+6);
        Xvn(i, 1) = x(v+7);
        Pcn(i, 1) = x(v+8);
        v = v + 9;
    end

    Tn = ET_TPARK(dn, ng);
%     An = inv(Tn)*inv(Mp)*Tn;
    Appn = inv(Tn)*inv(Mpp)*Tn;
    
    Eqdpn = zeros(2*ng, 1);
    Eqdppn = zeros(2*ng, 1);
    for i = 1:ng
        Eqdpn(2*i-1, 1) = Eqpn(i);
        Eqdpn(2*i, 1) = Edpn(i);
        Eqdppn(2*i-1, 1) = Eqppn(i);
        Eqdppn(2*i, 1) = Edppn(i);
    end

%     Ermn = inv(Tn)*Eqdpn;
    Ermn = inv(Tn)*Eqdppn;
    Vrmn = inv(YKrm + Appn)*Appn*Ermn;
    Irmn = YKrm*Vrmn;
    
    Iqdn = Tn*Irmn;
    
    v = 1;
    for i = 1:ng
        
        Iqn(i) = Iqdn(2*i-1);
        Idn(i) = Iqdn(2*i);
        
        Vtn(i) = Vrmn(2*i-1) + 1i*Vrmn(2*i);
        Itn(i) = Irmn(2*i-1) + 1i*Irmn(2*i);
        
        Pegapn(i) = real(Vtn(i)*conj(Itn(i))) + Ra(i)*abs(Itn(i))^2;
        Paceln(i) = Pmgapn(i) - Pegapn(i);
        Pamortn(i) = Ka*wn(i);

        Pacelv(i) = Pmgapv(i) - Pegapv(i);
        Pamortv(i) = Ka*wv(i);
        
        dwv(i) = (we/(2*H(i)))*(Pacelv(i) - Pamortv(i));
        ddv(i) = wv(i);
        dEqpv(i) = (1/Td0p(i))*(-Eqpv(i) + (Xd(i) - Xdp(i))*Idv(i) + Eexc(i));
        dEdpv(i) = (1/Tq0p(i))*(-Edpv(i) - (Xq(i) - Xqp(i))*Iqv(i));
        dEqppv(i) = (1/Td0pp(i))*(-Eqppv(i) + (Xdp(i) - Xdpp(i))*Idv(i) + Eqpv(i));
        dEdppv(i) = (1/Tq0pp(i))*(-Edppv(i) - (Xqp(i) - Xqpp(i))*Iqv(i) + Edpv(i));
        dPmgapv(i) = (1/Tt(i))*(Kv(i)*Xvv(i) - Pmgapv(i));
        
        Xpv(i) = -wv(i)/(2*pi*R(i));
        
        dXvv(i) = (1/Tv(i))*(Kt(i)*(Pcv(i) + Xpv(i)) - Xvv(i));
        dPcv(i) = -Ki(i)*wv(i)/(2*pi);
        
        dwn(i) = (we/(2*H(i)))*(Paceln(i) - Pamortn(i));
        ddn(i) = wn(i);
        dEqpn(i) = (1/Td0p(i))*(-Eqpn(i) + (Xd(i) - Xdp(i))*Idn(i) + Eexc(i));
        dEdpn(i) = (1/Tq0p(i))*(-Edpn(i) - (Xq(i) - Xqp(i))*Iqn(i));
        dEqppn(i) = (1/Td0pp(i))*(-Eqppn(i) + (Xdp(i) - Xdpp(i))*Idn(i) + Eqpn(i));
        dEdppn(i) = (1/Tq0pp(i))*(-Edppn(i) - (Xqp(i) - Xqpp(i))*Iqn(i) + Edpn(i));
        dPmgapn(i) = (1/Tt(i))*(Kv(i)*Xvn(i) - Pmgapn(i));
        
        Xpn(i) = -wn(i)/(2*pi*R(i));
        
        dXvn(i) = (1/Tv(i))*(Kt(i)*(Pcn(i) + Xpn(i)) - Xvn(i));
        dPcn(i) = -Ki(i)*wn(i)/(2*pi);
        
        F(v) = -wn(i) + wv(i) + (dt/2)*(dwn(i) + dwv(i));
        F(v+1) = -dn(i) + dv(i) + (dt/2)*(ddn(i) + ddv(i));
        F(v+2) = -Eqpn(i) + Eqpv(i) + (dt/2)*(dEqpn(i) + dEqpv(i));
        F(v+3) = -Edpn(i) + Edpv(i) + (dt/2)*(dEdpn(i) + dEdpv(i));
        F(v+4) = -Eqppn(i) + Eqppv(i) + (dt/2)*(dEqppn(i) + dEqppv(i));
        F(v+5) = -Edppn(i) + Edppv(i) + (dt/2)*(dEdppn(i) + dEdppv(i));
        F(v+6) = -Pmgapn(i) + Pmgapv(i) + (dt/2)*(dPmgapn(i) + dPmgapv(i));
        F(v+7) = -Xvn(i) + Xvv(i) + (dt/2)*(dXvn(i) + dXvv(i));
        F(v+8) = -Pcn(i) + Pcv(i) + (dt/2)*(dPcn(i) + dPcv(i));
        v = v + 9;
    end
end