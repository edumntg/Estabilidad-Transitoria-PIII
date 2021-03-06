function [Eq, Ed, Eqp, Edp, Eqpp, Edpp, Eexc] = ET_EQDEXC(Vq, Vd, Iq, Id, Ra, Xd, Xdp, Xdpp, Xq, Xqp, Xqpp, ng)

    for i = 1:ng
        Eq(i) = Vq(i) + Ra(i)*Iq(i) - Xd(i)*Id(i);
        Eqp(i) = Vq(i) + Ra(i)*Iq(i) - Xdp(i)*Id(i);
        Eqpp(i) = Vq(i) + Ra(i)*Iq(i) - Xdpp(i)*Id(i);
        
        Ed(i) = Vd(i) + Ra(i)*Id(i) + Xq(i)*Iq(i);
        Edp(i) = Vd(i) + Ra(i)*Id(i) + Xqp(i)*Iq(i);
        Edpp(i) = Vd(i) + Ra(i)*Id(i) + Xqpp(i)*Iq(i);
        
        Eexc(i) = abs(Eq(i) + 1i*Ed(i));
    end
end