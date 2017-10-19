function [Eq, Ed] = ET_EQD(Vq, Vd, Iq, Id, Ra, Xd, Xq, ng)

    for i = 1:ng
        Eq(i) = Vq(i) + Ra(i)*Iq(i) - Xd(i)*Id(i);
        Ed(i) = Vd(i) + Ra(i)*Id(i) - Xq(i)*Iq(i);
    end
end