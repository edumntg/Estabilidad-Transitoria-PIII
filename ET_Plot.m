function ET_Plot(wp, dp, Pep, wf, df, Pef, wpt, dpt, Pept, tvec)

    ng = size(dp, 1);
    ti = tvec(1);
    tp = tvec(2);
    td = tvec(3);
    tf = tvec(4);
    dt = tvec(5);

    plegend = cell(length(ng), 1);
    for i = 1:ng
        plegend{i} = strcat('Generador ', num2str(i));
    end

    figure(1), hold on
    for i = 1:ng
        plot(ti:dt:tp, Pep, tp:dt:td, Pef, td:dt:tf, Pept), grid minor
        line([tp tp], [Pef(i, 1) Pep(i, length(Pep))]);
        line([td td], [Pept(i, 1) Pef(i, length(Pef))]);
    end
    title('Potencia electrica');
    legend(plegend);
    hold off;

%     figure(2)
%     for i = 1:ng
%         plot(ti:dt:tp, wp, tp:dt:td, wf, td:dt:tf, wpt), grid minor
%     end
%     title('Desviacion de velocidad');
%     legend(plegend);
%     
%     figure(3)
%     for i = 1:ng
%         plot(ti:dt:tp, dp, tp:dt:td, df, td:dt:tf, dpt), grid minor
%     end
%     title('Angulo delta');
%     legend(plegend);

%     wtot = zeros(ng, length(ti:dt:tf-dt-dt));
    k = 1;
    for g = 1:ng
        kp = 1;
        kf = 1;
        kpt = 1;
        k = 1;
        for t = ti:dt:tp-dt
            wtot(g, k) = wp(g, kp);
            dtot(g, k) = dp(g, kp);
            kp = kp + 1;
            k = k + 1;
        end
        for t = tp:dt:td-dt
            wtot(g, k) = wf(g, kf);
            dtot(g, k) = df(g, kf);
            k = k + 1;
            kf = kf + 1;
        end
        for t = td:dt:tf
            wtot(g, k) = wpt(g, kpt);
            dtot(g, k) = dpt(g, kpt);
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
    
end