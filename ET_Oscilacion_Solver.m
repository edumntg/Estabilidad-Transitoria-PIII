function F = ET_Oscilacion_Solver(x, k, n_gen, t, td, dt, YKronPre, YKronFalla, YKronPost, YShuntKronPre, YShuntKronFalla, YShuntKronPost, Pmn, w, d, Pmv, Pev, Em, H, w0)

    % Tendremos dos ecuaciones por generador, en total 2*n_gen ecuaciones
    % Tambien tendremos 2 variables, las cuales seran wn y dn (n de nuevo)
    
    Gp = real(YKronPre);
    Bp = imag(YKronPre);
    gp = real(YShuntKronPre);
    
    Gf = real(YKronFalla);
    Bf = imag(YKronFalla);
    gf = real(YShuntKronFalla);
    
    Gpt = real(YKronPost);
    Bpt = imag(YKronPost);
    gpt = real(YShuntKronPost);
    
    Ka = 0.0;
    
    for i = 1:n_gen
        wn(i) = x(2*i - 1);
        dn(i) = x(2*i);
    end
    
    for i = 1:n_gen
        for j = 1:n_gen
            if(i ~= j)
                
                dni = dn(i);
                dnj = dn(j);
                
                if(t < td)      % Esta en condicion de falla
                    Peij(i, j) = (-Gf(i,j) + gf(i,j))*Em(i)^2 + Em(i)*Em(j)*(Gf(i,j)*cos(dni - dnj) + Bf(i,j)*sin(dni - dnj));
                else            % Se ha despejado la falla
                    Peij(i, j) = (-Gp(i,j) + gp(i,j))*Em(i)^2 + Em(i)*Em(j)*(Gp(i,j)*cos(dni - dnj) + Bp(i,j)*sin(dni - dnj));
                end
            end
        end
    end

    for i = 1:n_gen
        
        Pen(i) = sum(Peij(i, 1:n_gen));
        
        wni = wn(i);
        wvi = w(i, k - 1);
        
        dni = dn(i);
        dvi = d(i, k - 1);
        
        Pmni = Pmn(i);
        Pmvi = Pmv(i);
        
        Peni = Pen(i);
        Pevi = Pev(i, k - 1);
        
        Hi = H(i);
        
        F(2*i - 1) = wni - wvi - (dt/2)*(w0/(2*Hi))*(Pmni - Peni - Ka*wni + Pmvi - Pevi - Ka*wvi);
        F(2*i) = dni - dvi - (dt/2)*(wni + wvi);
    end
end