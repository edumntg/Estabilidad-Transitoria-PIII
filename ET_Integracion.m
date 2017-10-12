%% Aqui se haran las iteraciones y se resolvera la ecuacion de oscilacion

function [w, d, Pe] = ET_Integracion(H, E, d0, w0, Pe0, Pm, YKron, YShuntKron, f, tvec)

    ti = tvec(1);
    dt = tvec(2);
    tf = tvec(3);

    ng = size(YKron, 1);

    G = real(YKron);
    B = imag(YKron);
    g = real(YShuntKron);
%     b = imag(YShuntKron);

    options = optimset('Display', 'off');

    d = zeros(ng, length(ti:dt:tf));
    w = zeros(ng, length(ti:dt:tf));
    Pe = zeros(ng, length(ti:dt:tf));
    Peij = zeros(ng, ng);
    
    for gi = 1:ng
        w(gi, 1) = w0(gi);
        d(gi, 1) = d0(gi);
        Pe(gi, 1) = Pe0(gi);
    end

    k = 1;
    %% Se empieza a iterar
    for t = ti:dt:tf 
        for gi = 1:ng
            x0(2*gi - 1) = w(gi, k); % w0
            x0(2*gi) = d(gi, k); % d0
        end

        [x,fval,exitflag,output,jacobian] = fsolve(@(x)ET_EcuacionesFSOLVE(x, k, ng, dt, YKron, YShuntKron, Pm, w, d, Pm, Pe, E, H, 2*pi*f), x0, options);

        for gi = 1:ng
            w(gi, k) = x(2*gi - 1);
            d(gi, k) = x(2*gi);
        end

        for i = 1:ng
            for j = 1:ng
                if (i ~= j)
                    di = d(i, k);
                    dj = d(j, k);

                    Peij(i, j) = (-G(i, j) + g(i, j))*E(i)^2 + E(i)*E(j)*(G(i,j)*cos(di - dj) + B(i,j)*sin(di - dj));
                end
            end
            Pe(i, k) = sum(Peij(i, 1:ng));
        end
        k = k + 1;
    end
end