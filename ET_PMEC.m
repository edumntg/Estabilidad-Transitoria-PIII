function Pm = ET_PMEC(VRM, IRM, Ra, ng)
     
    for i = 1:ng
        Vt(i) = VRM(2*i-1) + 1i*VRM(2*i);
        It(i) = IRM(2*i-1) + 1i*IRM(2*i);

        Pm(i) = real(Vt(i)*conj(It(i))) + Ra(i)*abs(It(i))^2;
    end
end