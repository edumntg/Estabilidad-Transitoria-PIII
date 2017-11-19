function [Vqd, Iqd, Vq, Vd, Iq, Id] = ET_VIQD(VRM, IRM, T, ng)

    Vqd = T*VRM;
    Iqd = T*IRM;

    Vq = zeros(ng, 1);
    Vd = zeros(ng, 1);
    Iq = zeros(ng, 1);
    Id = zeros(ng, 1);
    for i = 1:ng
        Vq(i) = Vqd(2*i-1);
        Vd(i) = Vqd(2*i);
        Iq(i) = Iqd(2*i-1);
        Id(i) = Iqd(2*i);
    end
end