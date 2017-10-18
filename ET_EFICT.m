function [Ef, d] = ET_EFICT(VRM, IRM, GENDATA)
    
    %% En este punto, los vectores VRM e IRM solo contienen tensiones y corrientes de las barras SLACK y PV
     % Por tanto ya no hace falta realizar agregar una condicion para
     % excluir barras PQ
     
     % Se calculan los Vt e It
     n = size(GENDATA, 1);
     
     for i = 1:n
         Ra = GENDATA(i, 7);
         Xd = GENDATA(i, 8);
         Xq = GENDATA(i, 9);
         
         Vt(i) = VRM(2*i-1) + 1i*VRM(2*i);
         It(i) = IRM(2*i-1) + 1i*IRM(2*i);
         
         Ef(i) = Vt(i) + (Ra + 1i*Xq)*It(i);
         d(i) = angle(Ef(i));
     end
end