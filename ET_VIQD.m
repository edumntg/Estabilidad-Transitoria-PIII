function [Vqd, Iqd, Vq, Vd, Iq, Id] = ET_VIQD(VRM, IRM, T, ng)

    Vqd = T*VRM;
    Iqd = T*IRM;
    for i = 1:ng
        Vq(i) = Vqd(2*i-1);
        Vd(i) = Vqd(2*i);
        Iq(i) = Iqd(2*i-1);
        Id(i) = Iqd(2*i);
    end
end