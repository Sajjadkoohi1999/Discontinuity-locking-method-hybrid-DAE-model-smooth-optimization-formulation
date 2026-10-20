function res=model(t,x,xdot) 
global u1 u2 u3 u4 i P y14 y24 Cvo Cvn Cvc Cv4 T V R kb Po Pn Pc P4 Ck y11 y23 
%% Variables
nt=x(1); n1=x(2); n2=x(3); n3=x(4); F1=x(5); F2=x(6); F3=x(7); F4=x(8); y34=x(9); y14=x(10) ;y24=x(11) ;P=x(12); r1=x(13);
r2=x(14);  r3=x(15);
dnt=xdot(1); dn1=xdot(2); dn2=xdot(3); 
%% Equations
res=[y14-n1/nt;y24-n2/nt;y34-1+y14+y24;P-nt*T*R/V; nt-n1-n2-n3;
    F1-median([u1(i)*Cvc*sqrt((Pc+P)/2)*((Pc-P)/(sqrt(abs(Pc-P)+kb*Pc))),0,u1(i)*Ck*Cvc*Pc*0.85/sqrt(2)]);
    F2-median([u2(i)*Cvo*sqrt((Po+P)/2)*((Po-P)/(sqrt(abs(Po-P)+kb*Po))),0,u2(i)*Ck*Cvo*Po*0.85/sqrt(2)]);
    F3-median([u3(i)*Cvn*sqrt((Pn+P)/2)*((Pn-P)/(sqrt(abs(Pn-P)+kb*Pn))),0,u3(i)*Ck*Cvn*Pn*0.85/sqrt(2)]);
    F4-median([u4(i)*Cv4*sqrt((P+P4)/2)*((P-P4)/(sqrt(abs(P-P4)+kb*P))),0,u4(i)*Ck*Cv4*P*0.85/sqrt(2)]);
    dnt-F1-F2-F3+F4;dn1-F1*y11+F4*y14;dn2-F2*y23+F4*y24;r1-u2(i);r2-u3(i);r3-u4(i)];
end