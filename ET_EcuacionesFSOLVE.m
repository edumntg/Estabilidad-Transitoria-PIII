function F = ET_EcuacionesFSOLVE(x, k, ng, dt, YKron, YShuntKron, Pmn, w, d, Pmv, Pev, E, H, w0)

    % Tendremos dos ecuaciones por generador, en total 2*n_gen ecuaciones
    % Tambien tendremos 2 variables, las cuales seran wn y dn (n de nuevo)
    
    G = real(YKron);
    B = imag(YKron);
    g = real(YShuntKron);
%     b = imag(YShuntKron);
    
    Ka = 0.0;
    
    for i = 1:ng
        wn(i) = x(2*i - 1);
        dn(i) = x(2*i);
    end
    
    for i = 1:ng
        for j = 1:ng
            if(i ~= j)
                dni = dn(i);
                dnj = dn(j);
                Peij(i, j) = (-G(i,j) + g(i,j))*E(i)^2 + E(i)*E(j)*(G(i,j)*cos(dni - dnj) + B(i,j)*sin(dni - dnj));
            end
        end
    end

    for i = 1:ng
        
        Pen(i) = sum(Peij(i, 1:ng));
        
        wvi = w(i, length(w));
        dvi = d(i, length(d));
        Pevi = Pev(i, length(Pev));
        Pmvi = Pmv(i);
        
        wni = wn(i);
        dni = dn(i);
        Pmni = Pmn(i);
        Peni = Pen(i);
        
        Hi = H(i);
        
        F(2*i - 1) = wni - wvi - (dt/2)*(w0/(2*Hi))*(Pmni - Peni - Ka*wni + Pmvi - Pevi - Ka*wvi);
        F(2*i) = dni - dvi - (dt/2)*(wni + wvi);
    end
end