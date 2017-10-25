function [Ef, d, Vt, It] = ET_EFICT(VRM, IRM, Ra, Xq, Xd, ng)
    
    %% En este punto, los vectores VRM e IRM solo contienen tensiones y corrientes de las barras SLACK y PV
     % Por tanto ya no hace falta realizar agregar una condicion para
     % excluir barras PQ
     
     % Se calculan los Vt e It
     
     for i = 1:ng         
         Vt(i) = VRM(2*i-1) + 1i*VRM(2*i);
         It(i) = IRM(2*i-1) + 1i*IRM(2*i);
         
         Ef(i) = Vt(i) + (Ra(i) + 1i*Xq(i))*It(i);
         d(i) = angle(Ef(i));
     end
end