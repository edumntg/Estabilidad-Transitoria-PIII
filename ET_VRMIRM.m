function [VRM, IRM] = ET_VRMIRM(V, theta, YRM)
    
    nv = size(V, 1);
    for i = 1:nv
        Vrect(i) = V(i)*(cos(theta(i)) + 1i*sin(theta(i)));
        VRM(2*i-1) = real(Vrect(i));
        VRM(2*i) = imag(Vrect(i));
    end
    
    IRM = YRM*VRM;
    pause;
    
end