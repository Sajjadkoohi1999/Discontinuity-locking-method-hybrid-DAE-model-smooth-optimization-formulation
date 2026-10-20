function y=initialize(x,j) 
global u1 u2 u3 u4 i P y14 y24 Cvo Cvn Cvc Cv4 T V R kb Po Pn Pc P4 Ck y11 y23
%% Parameters
Cvo=8; % Valve coefficent of air stream
Cvn=8; % Valve coefficent of N2 stream
Cvc=8; % Valve coefficent of fuel stream
Cv4=8; % Valve coefficent of outlet stream
T=300; % Gas temperature in the vessel
V=3;   % Vessel volume
R=8.314E-5; % Gas constant
kb=1e-5; % Regularization constants
Po=12; % Air pressure
Pn=7;  % N2 pressure
Pc=10; % Fuel pressure
P4=2;  % Outlet pressure
Ck=(0.47*(1.53)^0.5)/(0.85*(sqrt(0.47+kb))); % Regularization constants
y11=1; % Mole fraction of CH4 in fuel stream
y32=1; % Mole fraction of N2 in inert stream
y23=0.29; % Mole fraction of O2 in air stream
y33=0.71; % Mole fraction of N2 in air stream
%% Variable
nt=x(1); % Total mole numbers in vessel
n1=x(2); % Mole numbers of CH4 in vessel
n2=x(3); % Mole numbers of O2 in vessel
n3=x(4); % Mole numbers of N2 in vessel
F1=x(5); % Flow rate of fuel stream
F2=x(6); % Flow rate of air stream
F3=x(7); % Flow rate of N2 stream
F4=x(8); % Flow rate of outlet stream
y34=x(9);% Mole fraction of N2 in vessel
ntdot=x(10); n1dot=x(11); n2dot=x(12); n3dot=x(13);
%% Algebric Equations
y(4)=nt-(P*V)/(R*T);
y(1)=n1-y14*nt;
y(2)=n2-y24*nt;
y(5)=n3-nt+n1+n2;
%% If we considered the check valve/choked flow features from the valve equations
if j==0
y(6)=F1-median([u1(i)*Cvc*sqrt((Pc+P)/2)*((Pc-P)/(sqrt(abs(Pc-P)+kb*Pc))),0,u1(i)*Ck*Cvc*Pc*0.85/sqrt(2)]);
y(7)=F2-median([u2(i)*Cvo*sqrt((Po+P)/2)*((Po-P)/(sqrt(abs(Po-P)+kb*Po))),0,u2(i)*Ck*Cvo*Po*0.85/sqrt(2)]);
y(8)=F3-median([u3(i)*Cvn*sqrt((Pn+P)/2)*((Pn-P)/(sqrt(abs(Pn-P)+kb*Pn))),0,u3(i)*Ck*Cvn*Pn*0.85/sqrt(2)]);
y(9)=F4-median([u4(i)*Cv4*sqrt((P+P4)/2)*((P-P4)/(sqrt(abs(P-P4)+kb*P))),0,u4(i)*Ck*Cv4*P*0.85/sqrt(2)]);
end
%% If we remove the check valve/choked flow features from the valve equations
if j==1
y(6)=F1-u1(i)*Cvc*sqrt((Pc+P)/2)*((Pc-P)/(sqrt(abs(Pc-P)+kb*Pc)));
y(7)=F2-u2(i)*Cvo*sqrt((Po+P)/2)*((Po-P)/(sqrt(abs(Po-P)+kb*Po)));
y(8)=F3-u3(i)*Cvn*sqrt((Pn+P)/2)*((Pn-P)/(sqrt(abs(Pn-P)+kb*Pn)));
y(9)=F4-u4(i)*Cv4*sqrt((P+P4)/2)*((P-P4)/(sqrt(abs(P-P4)+kb*P)));  
end
y(3)=1-y34-y24-y14;
%% ODE Equations
y(10)=ntdot-F1-F2-F3+F4;
y(11)=n1dot-F1*y11+F4*y14;
y(12)=n2dot-F2*y23+F4*y24;
y(13)=n3dot-F3*y32-F2*y33+F4*y34;
end