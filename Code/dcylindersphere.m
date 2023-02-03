function [Distance,P0,P1] = dcylindersphere(Ac, rc, C, rs)

Dc = diff(Ac,1,2);
t = Dc'*(C-Ac(:,1));
Dc2 = Dc'*Dc;
if and(t >=0,t <= Dc2)
    CA1 = C-Ac(:,1);
    t2 = t*t;
    t2 = t2/Dc2;
    D2 =  CA1'*CA1 - t2;
    Distance = ((max(D2,0))^0.5) - (rc+rs);
    P0 = Ac(:,1) + ((t2/Dc2)^0.5)*Dc;
    dP = C-P0;
    P0 = P0+dP*(rc/(dP(1)*dP(1)+dP(2)*dP(2)+dP(3)*dP(3))^0.5);
else
    Q = null(Dc');
    if t < 0
        P0 = Ac(:,1);
        Ac(:,1);
    else
        P0 = Ac(:,2);
        t = t-Dc2;
    end
    CA = (C-P0);
    t2 = t*t/Dc2;
    Cp = Q'*CA;
    Cp2 = Cp'*Cp;
    if Cp2 < rc*rc
        Pb = Cp;
        dCcap2 = 0;
    else
        Pb = Cp*(rc/((Cp2)^0.5));
        dCcap2 = (((Cp2)^0.5)-rc)^2;
    end
    P0 = P0 + Q*Pb;
    Distance = ((dCcap2 + t2)^0.5) - rs;
end

dP = P0-C;
P1 = C + dP*(rs/(dP(1)*dP(1)+dP(2)*dP(2)+dP(3)*dP(3))^0.5);

end
