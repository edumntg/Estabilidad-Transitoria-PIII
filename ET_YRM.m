% Eduardo Montilva 12-10089

% Programa para realizar la reduccion de Kron

function [YRM, YaRM, YbRM, YcRM, YdRM] = ET_YRM(Ykron)
    %% Debemos reducir el sistema para eliminar las barras PQ.
    
    Gk = real(Ykron);
    Bk = imag(Ykron);

    
    nk = size(Ykron, 1);
    YaRM = zeros(nk, nk);
    YbRM = zeros(nk, nk);
    %% Ahora esta matriz la convertimos a RM
    for i = 1:nk
        for j = 1:nk
            YaRM(i, j) = Gk(i, j);
            YbRM(i, j) = -Bk(i, j);
        end
    end
    YdRM = YaRM;
    YcRM = -YbRM;

    YRM = [YaRM YbRM; YcRM YdRM];
end