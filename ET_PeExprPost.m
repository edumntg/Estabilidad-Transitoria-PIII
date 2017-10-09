%% Este script calculara las expresiones de la ecuacion de potencia para los generadores

% dv = sym('dv', [1 n_gen]);

Gpt = real(YKronPost);
Bpt = imag(YKronPost);
gpt = real(YShuntKronPost);
bpt = imag(YShuntKronPost);

% PeEqPost = sym(zeros(n_gen, n_gen));
% for i = 1:n_gen
%     for j = 1:n_gen
%         if(i ~= j)
%             PeEqPost(i, j) = (-Gpt(i,j) + gpt(i,j))*Em(i)^2 + Em(i)*Em(j)*(Gpt(i,j)*cos(dv(i)-dv(j)) + Bpt(i,j)*sin(dv(i)-dv(j)));
%         end
%     end
% end

%% Aprovechamos y calculamos las potencias pre-falla que deben dar igual a lo que da el flujo de carga
Pe0Post = zeros(n_gen, n_gen);
for i = 1:n_gen
    for j = 1:n_gen
        if(i ~= j)
            Pe0Post(i, j) = (-Gpt(i,j) + gpt(i,j))*Em(i)^2 + Em(i)*Em(j)*(Gpt(i,j)*cos(d0(i)-d0(j)) + Bpt(i,j)*sin(d0(i)-d0(j)));
        end
    end
end