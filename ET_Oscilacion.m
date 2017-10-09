%% Aqui se haran las iteraciones y se resolvera la ecuacion de oscilacion

% Gp = real(YKronPre);
% Bp = imag(YKronPre);
% gp = real(YShuntKronPre);
% 
% Gf = real(YKronFalla);
% Bf = imag(YKronFalla);
% gf = real(YShuntKronFalla);
% 
% Gpt = real(YKronPost);
% Bpt = imag(YKronPost);
% gpt = real(YShuntKronPost);

k = 1;
for t = ti:dt:tp-dt
    % En esta parte aun no ocurre la falla, por tanto todas las variables
    % permanecen constantes
    for g = 1:n_gen
        w(g, k) = 0;
        d(g, k) = d0(g);
        Pe(g, k) = Pe0(g);
    end
    k = k + 1;
end


options = optimset('Display', 'off');

%% Situacion de Falla
for t = tp:dt:td-dt 
    for g = 1:n_gen
        x0(2*g - 1) = w(g, k - 1); % w0
        x0(2*g) = d(g, k - 1); % d0
    end
    
    x = fsolve(@(x)ET_Oscilacion_Solver(x, k, n_gen, t, td, dt, YKronPre, YKronFalla, YKronPost, YShuntKronPre, YShuntKronFalla, YShuntKronPost, Pmfalla, w, d, Pmfalla, Pe, Em, H, w0), x0, options);
    for g = 1:n_gen
        w(g, k) = x(2*g - 1);
        d(g, k) = x(2*g);
    end
    
    for i = 1:n_gen
        for j = 1:n_gen
            if (i ~= j)
                di = d(i, k);
                dj = d(j, k);

                Peij(i, j) = (-Gf(i, j) + gf(i, j))*Em(i)^2 + Em(i)*Em(j)*(Gf(i,j)*cos(di - dj) + Bf(i,j)*sin(di - dj));
            end
        end
        Pe(i, k) = sum(Peij(i, 1:n_gen));
    end
    k = k + 1;
end

%% Luego del despeje
for t = td:dt:tf
    for g = 1:n_gen
        x0(2*g - 1) = w(g, k - 1); % w0
        x0(2*g) = d(g, k - 1); % d0
    end
    
    x = fsolve(@(x)ET_Oscilacion_Solver(x, k, n_gen, t, td, dt, YKronPre, YKronFalla, YKronPost, YShuntKronPre, YShuntKronFalla, YShuntKronPost, Pmpost, w, d, Pmpost, Pe, Em, H, w0), x0, options);
    for g = 1:n_gen
        w(g, k) = x(2*g - 1);
        d(g, k) = x(2*g);
    end
    
    for i = 1:n_gen
        for j = 1:n_gen
            if (i ~= j)
                di = d(i, k);
                dj = d(j, k);

                Peij(i, j) = (-Gpt(i, j) + gpt(i, j))*Em(i)^2 + Em(i)*Em(j)*(Gpt(i,j)*cos(di - dj) + Bpt(i,j)*sin(di - dj));
            end
        end
        Pe(i, k) = sum(Peij(i, 1:n_gen));
    end
    k = k + 1;
end
