% Eduardo Montilva 12-10089

% Programa para realizar la reduccion de Kron

function YRM = ET_YRM(Ykron)
    %% Debemos reducir el sistema para eliminar las barras PQ.
    
    Gk = real(Ykron);
    Bk = imag(Ykron);

    
    nk = size(Ykron, 1);
    YRM = zeros(2*nk, 2*nk);
    
    %% Ahora esta matriz la convertimos a RM
    for i = 1:nk
        for j = 1:nk
            if(i == j)
                YRM(2*i-1, 2*i-1) = Gk(i,i);
                YRM(2*i, 2*i) = Gk(i,i);
                
                YRM(2*i-1, 2*i) = -Bk(i,i);
                YRM(2*i, 2*i-1) = Bk(i,i);
            end
            
            if(i ~= j)
                YRM(2*i-1, 2*j-1) = Gk(i,j);
                YRM(2*i, 2*j) = Gk(i,j);
                
                YRM(2*i-1, 2*j) = -Bk(i,j);
                YRM(2*i, 2*j-1) = Bk(i,j);
            end
        end
    end
end