function [SS,Y,WY,nb] = funcL2WMON(F,WF,Indec,L1,LN)
if Indec==-1; F=-F;end
nb=0;
J = L1;
IS = 1;
Z(IS) = F(J);
WZ(IS) = WF(J);
IW(IS) = 1;
SS(J) = 0;
while (J ~= LN)
    J = J + 1;
    IS = IS + 1;
    Z(IS) = F(J);
    WZ(IS) = WF(J);
    IW(IS) = 1;
    SS(J) = SS(J-1);
    while ((IS ~=1 ) &&(Z(IS - 1)>=Z(IS)))
        SS(J) = SS(J) + (Z(IS) - Z(IS - 1))^2*(WZ(IS)*WZ(IS - 1)/(WZ(IS) + WZ(IS - 1)));
        IS = IS - 1;
        Z(IS) = (WZ(IS)*Z(IS) + WZ(IS + 1)*Z(IS + 1))/(WZ(IS) + WZ(IS+1));
        WZ(IS) = WZ(IS) + WZ(IS + 1);
        IW(IS) = IW(IS) + IW(IS + 1);
    end
end
J = L1;
for I = 1:IS
    JJ = J + IW(I) - 1;
    while (J<=JJ)
        Y(J) = Z(I);
        WY(J) = WZ(I);
        J = J + 1;
    end
    nb=nb+1;
end
if Indec==-1; Y=-Y;end
end