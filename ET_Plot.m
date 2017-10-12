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

    figure(2)
    for i = 1:ng
        plot(ti:dt:tp, wp, tp:dt:td, wf, td:dt:tf, wpt), grid minor
    end
    title('Desviacion de velocidad');
    legend(plegend);
    
    figure(3)
    for i = 1:ng
        plot(ti:dt:tp, dp, tp:dt:td, df, td:dt:tf, dpt), grid minor
    end
    title('Angulo delta');
    legend(plegend);
%     line([tp tp], ylim, 'Color',[1 0 0]);
%     line([td td], ylim, 'Color', [0 1 0]);


%     figure(2)
%     for i = 1:ng
%         plot(t, d(i, :)), hold on
%     end
%     line([tp tp], ylim, 'Color',[1 0 0]), line([td td], ylim, 'Color', [0 1 0]), grid minor, title('Angulo de carga');
%     legend(plegend), hold off;
% 
%     figure(3)
%     for i = 1:ng
%         plot(t, Pe(i, :)), hold on
%     end
%     line([tp tp], ylim, 'Color',[1 0 0]), line([td td], ylim, 'Color', [0 1 0]), grid minor, title('Potencia electrica');
%     legend(plegend), hold off;
end