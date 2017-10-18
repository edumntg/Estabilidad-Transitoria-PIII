function [VRM, IRM] = ET_VRMIRM(V, theta, YRM, BUSDATA)
    
    %% Solo se toma en cuenta barras SLACK y PV?
    nv = size(V, 1);
    for i = 1:nv
        if(BUSDATA(i, 2) ~= 0) % NO ES PQ
            Vrect(i) = V(i)*(cos(theta(i)) + 1i*sin(theta(i)));
            VRM(2*i-1) = real(Vrect(i));
            VRM(2*i) = imag(Vrect(i));
        end
    end
    VRM = VRM'; % Se traspone ya que matlab crea el vector como un vector de columnas
    IRM = YRM*VRM;
end