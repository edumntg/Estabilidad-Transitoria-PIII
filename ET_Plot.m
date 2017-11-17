function ET_Plot(w_pre, d_pre, Pe_pre, Eqp_pre, Eqpp_pre, Edp_pre, Edpp_pre, Vt_pre, theta_pre, Pmgap_pre, Xv_pre, Pc_pre, w_falla, d_falla, Pe_falla, Eqp_falla, Eqpp_falla, Edp_falla, Edpp_falla, Vt_falla, theta_falla, Pmgap_falla, Xv_falla, Pc_falla, w_post, d_post, Pe_post, Eqp_post, Eqpp_post, Edp_post, Edpp_post, Vt_post, theta_post, Pmgap_post, Xv_post, Pc_post, ng, n, tvec)

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
    
    figure(2), hold on
    for i = 1:ng
        plot(ti:dt:tp, Eqp_pre, tp:dt:td, Eqp_falla, td:dt:tf, Eqp_post), grid minor
        line([tp tp], [Eqp_falla(i, 1) Eqp_pre(i, length(Eqp_pre))]);
        line([td td], [Eqp_post(i, 1) Eqp_falla(i, length(Eqp_falla))]);
    end
    title('Eqp');
    legend(plegend);
    hold off;
    
    figure(3), hold on
    for i = 1:ng
        plot(ti:dt:tp, Edp_pre, tp:dt:td, Edp_falla, td:dt:tf, Edp_post), grid minor
        line([tp tp], [Edp_falla(i, 1) Edp_pre(i, size(Edp_pre, 2))]);
        line([td td], [Edp_post(i, 1) Edp_falla(i, size(Edp_falla, 2))]);
    end
    title('Edp');
    legend(plegend);
    hold off;
    
    figure(4), hold on
    for i = 1:n
        plot(ti:dt:tp, Vt_pre, tp:dt:td, Vt_falla, td:dt:tf, Vt_post), grid minor
        line([tp tp], [Vt_falla(i, 1) Vt_pre(i, size(Vt_pre, 2))]);
        line([td td], [Vt_post(i, 1) Vt_falla(i, size(Vt_falla, 2))]);
    end
    title('Voltajes en barras');
    %     legend(plegend);
    hold off;
    
    figure(5), hold on
    for i = 1:n
        plot(ti:dt:tp, theta_pre, tp:dt:td, theta_falla, td:dt:tf, theta_post), grid minor
        line([tp tp], [theta_falla(i, 1) theta_pre(i, size(theta_pre, 2))]);
        line([td td], [theta_post(i, 1) theta_falla(i, size(theta_falla, 2))]);
    end
    title('Angulos en barras');
%     legend(plegend);
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
            Eqpptot(i, k) = Eqpp_pre(i, kp);
            Edpptot(i, k) = Edpp_pre(i, kp);
            Pmgaptot(i, k) = Pmgap_pre(i, kp);
            Xvtot(i, k) = Xv_pre(i, kp);
            Pctot(i, k) = Pc_pre(i, kp);
            kp = kp + 1;
            k = k + 1;
        end
        for t = tp:dt:td-dt
            wtot(i, k) = w_falla(i, kf);
            dtot(i, k) = d_falla(i, kf);
            Eqpptot(i, k) = Eqpp_falla(i, kf);
            Edpptot(i, k) = Edpp_falla(i, kf);
            Pmgaptot(i, k) = Pmgap_falla(i, kf);
            Xvtot(i, k) = Xv_falla(i, kf);
            Pctot(i, k) = Pc_falla(i, kf);
            k = k + 1;
            kf = kf + 1;
        end
        for t = td:dt:tf
            wtot(i, k) = w_post(i, kpt);
            dtot(i, k) = d_post(i, kpt);
            Eqpptot(i, k) = Eqpp_post(i, kpt);
            Edpptot(i, k) = Edpp_post(i, kpt);
            Pmgaptot(i, k) = Pmgap_post(i, kpt);
            Xvtot(i, k) = Xv_post(i, kpt);
            Pctot(i, k) = Pc_post(i, kpt);
            k = k + 1;
            kpt = kpt + 1;
        end
    end
    
    figure(6)
    for i = 1:ng
        plot(ti:dt:tf, wtot), grid minor
    end
    title('Desviacion de velocidad');
    legend(plegend);
    
    figure(7)
    for i = 1:ng
        plot(ti:dt:tf, dtot), grid minor
    end
    title('Angulo delta');
    legend(plegend);
    
    figure(8)
    for i = 1:ng
        plot(ti:dt:tf, Eqpptot), grid minor
    end
    title('Eqpp');
    legend(plegend);
    
    figure(9)
    for i = 1:ng
        plot(ti:dt:tf, Edpptot), grid minor
    end
    title('Edpp');
    legend(plegend);
    
    figure(10)
    for i = 1:ng
        plot(ti:dt:tf, Pmgaptot), grid minor
    end
    title('Pmgap');
    legend(plegend);
    
    figure(11)
    for i = 1:ng
        plot(ti:dt:tf, Xvtot), grid minor
    end
    title('Xv');
    legend(plegend);
    
    figure(12)
    for i = 1:ng
        plot(ti:dt:tf, Pctot), grid minor
    end
    title('Pc');
    legend(plegend);
 
end