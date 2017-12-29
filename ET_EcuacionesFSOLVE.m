function F = ET_EcuacionesFSOLVE(x, YKrm, wv, dv, Eqpv, Eqppv, Edpv, Edppv, Pmgapv, Xvv, Pcv, Vmedv, Vviv, Vav, Vwv, Vpc1v, Vpc2v, Vexcv, Pegapv, Iqv, Idv, Vtv, Vref, ...
                                Mp, Mpp, Ra, Xq, Xqp, Xqpp, Xd, Xdp, Xdpp, Tq0p, Tq0pp, Td0p, Td0pp, ...
                                Tt, Kv, Tv, Kt, Ki, R, ...
                                Kmed, Kexc, Ka, Tmed, Texc, Ta, Kd, Kp, Kvi, ...
                                Kest, Tw, T1, T2, T3, T4, ...
                                Eexc, H, JAULA, AGC, AVR, PSS, we, ng, dt)

    % Tendremos 9 ecuaciones por generador, en total 4*n_gen ecuaciones
    % Tambien tendremos 9 variables, las cuales seran wn, dn, eqp, edp (n de nuevo)
    D = 0.0;
    
    wn = zeros(ng, 1);
    dn = zeros(ng, 1);
    Eqpn = zeros(ng, 1);
    Edpn = zeros(ng, 1);
    Eqppn = zeros(ng, 1);
    Edppn = zeros(ng, 1);
    Pmgapn = zeros(ng, 1);
    Xvn = zeros(ng, 1);
    Pcn = zeros(ng, 1);
    Vmedn = zeros(ng, 1);
    Vvin = zeros(ng, 1);
    Van = zeros(ng, 1);
    Vexcn = zeros(ng, 1);
    Vwn = zeros(ng, 1);
    Vpc1n = zeros(ng, 1);
    Vpc2n = zeros(ng, 1);
    
    Iqn = zeros(ng, 1);
    Idn = zeros(ng, 1);
    Vtn = zeros(ng, 1);
    Itn = zeros(ng, 1);
    Pegapn = zeros(ng, 1);
    Paceln = zeros(ng, 1);
    Pamortn = zeros(ng, 1);
    Pacelv = zeros(ng, 1);
    Pamortv = zeros(ng, 1);
    dwv = zeros(ng, 1);
    ddv = zeros(ng, 1);
    dEqpv = zeros(ng, 1);
    dEdpv = zeros(ng, 1);
    dEqppv = zeros(ng, 1);
    dEdppv = zeros(ng, 1);
    dPmgapv = zeros(ng, 1);
    Xpv = zeros(ng, 1);
    dXvv = zeros(ng, 1);
    dPcv = zeros(ng, 1);
    dVmedv = zeros(ng, 1);
    dVviv = zeros(ng, 1);
    dVav = zeros(ng, 1);
    dVexcv = zeros(ng, 1);
    dVwv = zeros(ng, 1);
    dVpc1v = zeros(ng, 1);
    dVpc2v = zeros(ng, 1);
    
    dwn = zeros(ng, 1);
    ddn = zeros(ng, 1);
    dEqpn = zeros(ng, 1);
    dEdpn = zeros(ng, 1);
    dEqppn = zeros(ng, 1);
    dEdppn = zeros(ng, 1);
    dPmgapn = zeros(ng, 1);
    Xpn = zeros(ng, 1);
    dXvn = zeros(ng, 1);
    dPcn = zeros(ng, 1);
    dVmedn = zeros(ng, 1);
    dVvin = zeros(ng, 1);
    dVan = zeros(ng, 1);
    dVexcn = zeros(ng, 1);
    dVwn = zeros(ng, 1);
    dVpc1n = zeros(ng, 1);
    dVpc2n = zeros(ng, 1);
    
    vv = 4;
    vv = vv + sum(JAULA.*2);
    vv = vv + sum(AGC.*3);
    vv = vv + sum(AVR.*4);
    vv = vv + sum(PSS.*3);
    
    F = zeros(vv, 1);
        
    v = 1;
    for i = 1:ng
        wn(i) = x(v);
        dn(i) = x(v+1);
        Eqpn(i) = x(v+2);
        Edpn(i) = x(v+3);
        v = v + 4;
        if(JAULA(i) == 1)
            Eqppn(i) = x(v);
            Edppn(i) = x(v+1);
            v = v + 2;
        end
        if(AGC(i) == 1)
            Pmgapn(i) = x(v);
            Xvn(i) = x(v+1);
            Pcn(i) = x(v+2);
            v = v + 3;
        else
            Pmgapn(i) = Pmgapv(i);
        end
        if(AVR(i) == 1)
            Vmedn(i) = x(v);
            Vvin(i) = x(v+1);
            Van(i) = x(v+2);
            Vexcn(i) = x(v+3);
            v = v + 4;
        else
            Vexcn(i) = Eexc(i);
        end
        if(PSS(i) == 1);
            Vwn(i) = x(v);
            Vpc1n(i) = x(v+1);
            Vpc2n(i) = x(v+2);
            v = v + 3;
        end
    end

    Tn = ET_TPARK(dn, ng);
    An = Tn\(Mp\Tn);
%     Appn = Tn\(Mpp\Tn);
    
    Eqdpn = zeros(2*ng, 1);
    Eqdppn = zeros(2*ng, 1);
    Eqd = zeros(2*ng, 1);
    for i = 1:ng
        Eqdpn(2*i-1, 1) = Eqpn(i);
        Eqdpn(2*i, 1) = Edpn(i);
        Eqdppn(2*i-1, 1) = Eqppn(i);
        Eqdppn(2*i, 1) = Edppn(i);
        
        Eqd(2*i-1, 1) = Eqpn(i);
        Eqd(2*i, 1) = Edpn(i);
        if(JAULA(i) == 1)
            Eqd(2*i-1, 1) = Eqppn(i);
            Eqd(2*i, 1) = Edppn(i);
        end
        
    end

%     Ermn = inv(Tn)*Eqdpn;
    Ermn = Tn\Eqd;
    Vrmn = (YKrm + An)\(An*Ermn);
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
        Pamortn(i) = D*wn(i);

        Pacelv(i) = Pmgapv(i) - Pegapv(i);
        Pamortv(i) = D*wv(i);
        
        dwv(i) = (1/(2*H(i)))*(Pacelv(i) - Pamortv(i));
        ddv(i) = wv(i)*we;
        
        dVmedv(i) = (1/Tmed(i))*(Kmed(i)*abs(Vtv(i)) - Vmedv(i));
        if(PSS(i) == 1)
            ev(i) = Vref(i) - Vmedv(i) + Vpc2v(i);
        else
            ev(i) = Vref(i) - Vmedv(i);
        end
        dVviv(i) = Kvi(i)*ev(i);
        
        %% Ecuaciones del PSS
        dVwv(i) = Kest(i)*dwv(i) - Vwv(i)/Tw(i);
        dVpc1v(i) = (1/T2(i))*(Vwv(i) + T1(i)*dVwv(i) - Vpc1v(i));
        dVpc2v(i) = (1/T4(i))*(Vpc1v(i) + T3(i)*dVpc1v(i) - Vpc2v(i));
                
        Vvpv(i) = Kp(i)*ev(i);

        if(PSS(i) == 1)
            Vvdv(i) = Kd(i)*(-dVmedv(i) + dVpc2v(i));
        else
            Vvdv(i) = -Kd(i)*dVmedv(i);
        end
        Vxv(i) = Vvpv(i) + Vvdv(i) + Vviv(i);
        dVav(i) = (1/Ta(i))*(Ka(i)*Vxv(i) - Vav(i));
        
        dVexcv(i) = (1/Texc(i))*(Kexc(i)*Vav(i) - Vexcv(i));
        
        dEqpv(i) = (1/Td0p(i))*(-Eqpv(i) + (Xd(i) - Xdp(i))*Idv(i) + Vexcv(i));
        dEdpv(i) = (1/Tq0p(i))*(-Edpv(i) - (Xq(i) - Xqp(i))*Iqv(i));
        dEqppv(i) = (1/Td0pp(i))*(-Eqppv(i) + (Xdp(i) - Xdpp(i))*Idv(i) + Eqpv(i));
        dEdppv(i) = (1/Tq0pp(i))*(-Edppv(i) - (Xqp(i) - Xqpp(i))*Iqv(i) + Edpv(i));
        dPmgapv(i) = (1/Tt(i))*(Kt(i)*Xvv(i) - Pmgapv(i));
        
        Xpv(i) = -wv(i)/R(i);
        
        dXvv(i) = (1/Tv(i))*(Kv(i)*(Pcv(i) + Xpv(i)) - Xvv(i));
        dPcv(i) = -Ki(i)*wv(i);
        
        dwn(i) = (1/(2*H(i)))*(Paceln(i) - Pamortn(i));
        ddn(i) = wn(i)*we;
        
        dVmedn(i) = (1/Tmed(i))*(Kmed(i)*abs(Vtn(i)) - Vmedn(i));
        if(PSS(i) == 1)
            en(i) = Vref(i) - Vmedn(i) + Vpc2n(i);
        else
            en(i) = Vref(i) - Vmedn(i);
        end
        dVvin(i) = Kvi(i)*en(i);

        dVwn(i) = Kest(i)*dwn(i) - Vwn(i)/Tw(i);
        dVpc1n(i) = (1/T2(i))*(Vwn(i) + T1(i)*dVwn(i) - Vpc1n(i));
        dVpc2n(i) = (1/T4(i))*(Vpc1n(i) + T3(i)*dVpc1n(i) - Vpc2n(i));
        
        Vvpn(i) = Kp(i)*en(i);
        if(PSS(i) == 1)
            Vvdn(i) = Kd(i)*(-dVmedn(i) + dVpc2n(i));
        else
            Vvdn(i) = -Kd(i)*dVmedn(i);
        end
        Vxn(i) = Vvpn(i) + Vvdn(i) + Vvin(i);
        dVan(i) = (1/Ta(i))*(Ka(i)*Vxn(i) - Van(i));
        
        dVexcn(i) = (1/Texc(i))*(Kexc(i)*Van(i) - Vexcn(i));        
        
        dEqpn(i) = (1/Td0p(i))*(-Eqpn(i) + (Xd(i) - Xdp(i))*Idn(i) + Vexcn(i));
        dEdpn(i) = (1/Tq0p(i))*(-Edpn(i) - (Xq(i) - Xqp(i))*Iqn(i));
        dEqppn(i) = (1/Td0pp(i))*(-Eqppn(i) + (Xdp(i) - Xdpp(i))*Idn(i) + Eqpn(i));
        dEdppn(i) = (1/Tq0pp(i))*(-Edppn(i) - (Xqp(i) - Xqpp(i))*Iqn(i) + Edpn(i));
        dPmgapn(i) = (1/Tt(i))*(Kt(i)*Xvn(i) - Pmgapn(i));
        
        Xpn(i) = -wn(i)/R(i);
        
        dXvn(i) = (1/Tv(i))*(Kv(i)*(Pcn(i) + Xpn(i)) - Xvn(i));
        dPcn(i) = -Ki(i)*wn(i);
               
        F(v) = -wn(i) + wv(i) + (dt/2)*(dwn(i) + dwv(i));
        F(v+1) = -dn(i) + dv(i) + (dt/2)*(ddn(i) + ddv(i));
        F(v+2) = -Eqpn(i) + Eqpv(i) + (dt/2)*(dEqpn(i) + dEqpv(i));
        F(v+3) = -Edpn(i) + Edpv(i) + (dt/2)*(dEdpn(i) + dEdpv(i));
        v = v + 4;
        if(JAULA(i) == 1)
            F(v) = -Eqppn(i) + Eqppv(i) + (dt/2)*(dEqppn(i) + dEqppv(i));
            F(v+1) = -Edppn(i) + Edppv(i) + (dt/2)*(dEdppn(i) + dEdppv(i));
            v = v + 2;
        end
        if(AGC(i) == 1)
            F(v) = -Pmgapn(i) + Pmgapv(i) + (dt/2)*(dPmgapn(i) + dPmgapv(i));
            F(v+1) = -Xvn(i) + Xvv(i) + (dt/2)*(dXvn(i) + dXvv(i));
            F(v+2) = -Pcn(i) + Pcv(i) + (dt/2)*(dPcn(i) + dPcv(i));
            v = v + 3;
        end
        if(AVR(i) == 1)
            F(v) = -Vmedn(i) + Vmedv(i) + (dt/2)*(dVmedn(i) + dVmedv(i));
            F(v+1) = -Vvin(i) + Vviv(i) + (dt/2)*(dVvin(i) + dVviv(i));
            F(v+2) = -Van(i) + Vav(i) + (dt/2)*(dVan(i) + dVav(i));
            F(v+3) = -Vexcn(i) + Vexcv(i) + (dt/2)*(dVexcn(i) + dVexcv(i));
            v = v + 4;
        end
        if(PSS(i) == 1)
            F(v) = -Vwn(i) + Vwv(i) + (dt/2)*(dVwn(i) + dVwv(i));
            F(v+1) = -Vpc1n(i) + Vpc1v(i) + (dt/2)*(dVpc1n(i) + dVpc1v(i));
            F(v+2) = -Vpc2n(i) + Vpc2v(i) + (dt/2)*(dVpc2n(i) + dVpc2v(i));
            v = v + 3;
        end
    end
end