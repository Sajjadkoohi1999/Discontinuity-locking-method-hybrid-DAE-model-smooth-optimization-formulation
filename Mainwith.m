clear all
close all
clc
global u1 u2 u3 u4 i P y14 y24 
% Initial condition is as follows:
y14=1; % Mole fraction of CH4 in vessel
y24=0; % Mole fraction of O2 in vessel
P=7.5; % Pressure in the vessel
i=1;
% Valve opening trajectories
u1=zeros(1,6); % Fuel gas valve opening %
u2=[0.0875 0.055 0.2 0.048 0.007 0.11]; % Air valve opening %
u3=[0 0.95 0 0 0 0.97]; % N2 valve opening %
u4=[0.075 0.1875 0 0.9888 0.49 1]; % CH4 valve opening %
o=6; %Evenly spaced control intervals
TT=[0,300]; % Simulation time
Time = TT(1);
c=ones(13,1); % Initial guesses for fsolve
[d,z,exitflag]=fsolve(@ (x)initialize (x,0),c);
a=[d(1:9);y14;y24;P;u2(i);u3(i);u4(i)]; % Initial conditions for ode15i
b=[d(10:13);zeros(11,1)]; % Initial conditions for ode15i
Tspan=[Time (TT(end)*i)/o];
Yplot=[];
Tplot=[];
while (Time<TT(end))
Tspan=[Time (TT(end)*i)/o];
x0= a;
xdot0=b;
[t,y]=ode15i(@model,Tspan,x0,xdot0);
Time=t(end);
Tplot=[Tplot;t];
Yplot=[Yplot;y];
y14=y(end,10);
y24=y(end,11);
P=y(end,12);
if i==6 
    break
end
i=i+1;
c=[y(end,1:9)';zeros(4,1)];
% Re-initialization
[d,z,exitflag]=fsolve(@ (x)initialize (x,0),c);
a=[d(1:9);y14;y24;P;u2(i);u3(i);u4(i)];
b=[d(10:13);zeros(11,1)];
end
%% plot
subplot(2,2,1)
plot(Tplot,Yplot(:,13)*100,'b-.',Tplot,Yplot(:,14)*100,'r-',Tplot,Yplot(:,15)*100,'k--')
xlabel('Time (s)')
axis( [0, 300, 0,100])
xticks([0 50 100 150 200 250 300])
xticklabels({'0','50','100','150','200','250','300'})
ylabel('Valve op%')
legend('Air','N2','Outlet')
subplot(2,2,2)
plot(Tplot,Yplot(:,12),'b')
axis( [0, 300, 1,14])
yticks([0:2:14])
yticklabels({'0','2','4','6','8','10','12','14'})
xlabel('Time (s)')
ylabel('Pressure (bar)')
legend('Pressure ')
subplot(2,2,3)
plot(Tplot,Yplot(:,6),'b-.',Tplot,Yplot(:,7),'r-',Tplot,Yplot(:,8),'k--')
xlabel('Time (s)')
ylabel('Flow rate (mole/s)')
legend('Air','N2','Outlet')
axis( [0, 300, -20,60] )
yticks([-20 -10  0 10 20 30 40 50 60])
yticklabels({'-20','-10','0','10','20','30','40','50','60'})
subplot(2,2,4)
plot(Tplot,Yplot(:,10),'r-',Tplot,Yplot(:,11),'b-.',Tplot,Yplot(:,9),'k--')
xlabel('Time (s)')
ylabel('Mole fraction')
legend('CH4','O2','N2')
axis( [0, 300, -0.02, 1] )
yticks([0:0.1:1])
yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'})
%% Ternary plot
figure
A= Yplot(:,10);
B = Yplot(:,11);
C=1-A-B;
ternplot(A,B,C)
ternlabel('y CH4', 'y O2', 'y N2')
hold on 
oi=[0.07 0.97];
qq=0.03;
ternplot(qq,oi)
hold on
oi=[0.07:0.005:.37];
qq=2*oi-0.11;
ternplot(qq,oi)
legend('Chanchoveer path','Lower Approx. Flam. Zone','Upper Approx. Flam. Zone')
