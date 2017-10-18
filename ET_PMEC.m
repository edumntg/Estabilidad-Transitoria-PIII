function Pm = ET_PMEC(GENDATA, VRM, IRM)
    
    n = size(GENDATA, 1);
     
     for i = 1:n
         Ra = GENDATA(i, 7);
         
         Vt(i) = VRM(2*i-1) + 1i*VRM(2*i);
         It(i) = IRM(2*i-1) + 1i*IRM(2*i);
         
         Pm(i) = real(Vt(i)*conj(It(i)) + Ra*abs(It(i))^2);
     end
end