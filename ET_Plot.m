function ET_Plot(w_pre, d_pre, Pe_pre, Eqp_pre, Edp_pre, Vt_pre, w_falla, d_falla, Pe_falla, Eqp_falla, Edp_falla, Vt_falla, w_post, d_post, Pe_post, Eqp_post, Edp_post, Vt_post, ng, tvec)

    f = 60;
    
    ti = tvec(1)*f;
    tp = tvec(2)*f;
    td = tvec(3)*f;
    tf = tvec(4)*f;
    dt = tvec(5)*f;

    plegend = cell(length(ng), 1);
    for i = 1:ng
        plegend{i} = strcat('Generador ', num2str(i));
    end

    figure(1), hold on
    for i = 1:ng
        plot(ti:dt:tp, Pe_pre, tp:dt:td, Pe_falla, td:dt:tf, Pe_post), grid minor
        line([tp tp], [Pe_falla(i, 1) Pe_pre(i, length(Pe_pre))]);
        line([td td], [Pe_post(i, 1) Pe_falla(i, length(Pe_falla))]);
    end
    title('Potencia electrica');
    legend(plegend);
    hold off;
    
    k = 1;
    for i = 1:ng
        kp = 1;
        kf = 1;
        kpt = 1;
        k = 1;
        for t = ti:dt:tp-dt
            wtot(i, k) = w_pre(i, kp);
            dtot(i, k) = d_pre(i, kp);
            Eqptot(i, k) = Eqp_pre(i, kp);
            Edptot(i, k) = Edp_pre(i, kp);
            kp = kp + 1;
            k = k + 1;
        end
        for t = tp:dt:td-dt
            wtot(i, k) = w_falla(i, kf);
            dtot(i, k) = d_falla(i, kf);
            Eqptot(i, k) = Eqp_falla(i, kf);
            Edptot(i, k) = Edp_falla(i, kf);
            k = k + 1;
            kf = kf + 1;
        end
        for t = td:dt:tf
            wtot(i, k) = w_post(i, kpt);
            dtot(i, k) = d_post(i, kpt);
            Eqptot(i, k) = Eqp_post(i, kpt);
            Edptot(i, k) = Edp_post(i, kpt);
            k = k + 1;
            kpt = kpt + 1;
        end
    end
    
    figure(2)
    for i = 1:ng
        plot(ti:dt:tf, wtot), grid minor
    end
    title('Desviacion de velocidad');
    legend(plegend);
    
    figure(3)
    for i = 1:ng
        plot(ti:dt:tf, dtot), grid minor
    end
    title('Angulo delta');
    legend(plegend);
    
    figure(4)
    for i = 1:ng
        plot(ti:dt:tf, Eqptot), grid minor
    end
    title('Eqp');
    legend(plegend);
    
    figure(5)
    for i = 1:ng
        plot(ti:dt:tf, Edptot), grid minor
    end
    title('Edp');
    legend(plegend);
    
    figure(6), hold on
    for i = 1:ng
        plot(ti:dt:tp, Vt_pre, tp:dt:td, Vt_falla, td:dt:tf, Vt_post), grid minor
        line([tp tp], [Vt_falla(i, 1) Vt_pre(i, length(Vt_pre))]);
        line([td td], [Vt_post(i, 1) Vt_falla(i, length(Vt_falla))]);
    end
    title('Voltajes en barras');
%     legend(plegend);
    hold off;
    
    
end