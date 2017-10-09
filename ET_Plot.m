

t = ti:dt:tf;

plegend = cell(length(n_gen), 1);
for i = 1:n_gen
    plegend{i} = strcat('Generador ', num2str(i));
end

figure(1)
for i = 1:n_gen
    plot(t, w(i, :)), hold on
end
line([tp tp], ylim, 'Color',[1 0 0]), line([td td], ylim, 'Color', [0 1 0]), grid minor, title('Desviacion de velocidad');
legend(plegend), hold off;

figure(2)
for i = 1:n_gen
    plot(t, d(i, :)), hold on
end
line([tp tp], ylim, 'Color',[1 0 0]), line([td td], ylim, 'Color', [0 1 0]), grid minor, title('Angulo de carga');
legend(plegend), hold off;

figure(3)
for i = 1:n_gen
    plot(t, Pe(i, :)), hold on
end
line([tp tp], ylim, 'Color',[1 0 0]), line([td td], ylim, 'Color', [0 1 0]), grid minor, title('Potencia electrica');
legend(plegend), hold off;
