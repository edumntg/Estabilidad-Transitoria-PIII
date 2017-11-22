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
        plot(ti:dt:tp, Pe_pre, tp:dt:td, Pe_falla, td:dt:tf, Pe_post, 'linewidth', 1.5), grid minor
        line([tp tp], [Pe_falla(i, 1) Pe_pre(i, length(Pe_pre))]);
        line([td td], [Pe_post(i, 1) Pe_falla(i, length(Pe_falla))]);
    end
    title('Potencia electrica');
    legend(plegend);
    hold off;
    
    figure(2), hold on
    for i = 1:ng
        plot(ti:dt:tp, Eqp_pre, tp:dt:td, Eqp_falla, td:dt:tf, Eqp_post, 'linewidth', 1.5), grid minor
        line([tp tp], [Eqp_falla(i, 1) Eqp_pre(i, length(Eqp_pre))]);
        line([td td], [Eqp_post(i, 1) Eqp_falla(i, length(Eqp_falla))]);
    end
    title('Eqp');
    legend(plegend);
    hold off;
    
    figure(3), hold on
    for i = 1:ng
        plot(ti:dt:tp, Edp_pre, tp:dt:td, Edp_falla, td:dt:tf, Edp_post, 'linewidth', 1.5), grid minor
        line([tp tp], [Edp_falla(i, 1) Edp_pre(i, size(Edp_pre, 2))]);
        line([td td], [Edp_post(i, 1) Edp_falla(i, size(Edp_falla, 2))]);
    end
    title('Edp');
    legend(plegend);
    hold off;
    
    figure(4), hold on
    for i = 1:n
        plot(ti:dt:tp, Vt_pre, tp:dt:td, Vt_falla, td:dt:tf, Vt_post, 'linewidth', 1.5), grid minor
        line([tp tp], [Vt_falla(i, 1) Vt_pre(i, size(Vt_pre, 2))]);
        line([td td], [Vt_post(i, 1) Vt_falla(i, size(Vt_falla, 2))]);
    end
    title('Voltajes en barras');
    %     legend(plegend);
    hold off;
    
    figure(5), hold on
    for i = 1:n
        plot(ti:dt:tp, theta_pre, tp:dt:td, theta_falla, td:dt:tf, theta_post, 'linewidth', 1.5), grid minor
        line([tp tp], [theta_falla(i, 1) theta_pre(i, size(theta_pre, 2))]);
        line([td td], [theta_post(i, 1) theta_falla(i, size(theta_falla, 2))]);
    end
    title('Angulos en barras');
%     legend(plegend);
    hold off;
    
    figure(6)
    for i = 1:ng
        plot(ti:dt:tp, w_pre, tp:dt:td, w_falla, td:dt:tf, w_post, 'linewidth', 1.5), grid minor
    end
    title('Desviacion de velocidad');
    legend(plegend);
    
    figure(7)
    for i = 1:ng
        plot(ti:dt:tp, d_pre, tp:dt:td, d_falla, td:dt:tf, d_post, 'linewidth', 1.5), grid minor
    end
    title('Angulo delta');
    legend(plegend);
    
    figure(8)
    for i = 1:ng
        plot(ti:dt:tp, Eqpp_pre, tp:dt:td, Eqpp_falla, td:dt:tf, Eqpp_post, 'linewidth', 1.5), grid minor
    end
    title('Eqpp');
    legend(plegend);
    
    figure(9)
    for i = 1:ng
        plot(ti:dt:tp, Edpp_pre, tp:dt:td, Edpp_falla, td:dt:tf, Edpp_post, 'linewidth', 1.5), grid minor
    end
    title('Edpp');
    legend(plegend);
    
    figure(10)
    for i = 1:ng
        plot(ti:dt:tp, Pmgap_pre, tp:dt:td, Pmgap_falla, td:dt:tf, Pmgap_post, 'linewidth', 1.5), grid minor
    end
    title('Pmgap');
    legend(plegend);
    
    figure(11)
    for i = 1:ng
        plot(ti:dt:tp, Xv_pre, tp:dt:td, Xv_falla, td:dt:tf, Xv_post, 'linewidth', 1.5), grid minor
    end
    title('Xv');
    legend(plegend);
    
    figure(12)
    for i = 1:ng
        plot(ti:dt:tp, Pc_pre, tp:dt:td, Pc_falla, td:dt:tf, Pc_post, 'linewidth', 1.5), grid minor
    end
    title('Pc');
    legend(plegend);
 
end