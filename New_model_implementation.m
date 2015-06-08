clear all; close all; clc

%%%Initial values
BPS0=1000e-9*6e23;
BP0=0;
BS0=0;
IPS0=0;
IP0=0;
IS0=0;
EPS0=0;
EP0=0;
ES0=0;
CPS0=0;
CP0=0;
CS0=0;
CSC0=0;

X0=[BPS0; BP0; BS0; IPS0; IP0; IS0; EPS0; EP0; ES0; CPS0; CP0; CS0; CSC0];

tottime = 100;
day = 24;
tspan=[0,day*tottime];

k1=600*60;              %Literature
k2=3.82e-7*3600;        %Literature
k3=1e-2*1e-3;           %Barlett - ktransblood*partition
k4=0.09*60;             %Literature
k5=8.7e-2;              %Literature
k6=1.6e-5*3600;         %Literature
k7=5e-1;                %Barlett - kdegenda
k8=0;                   %Assume endosmal escape is much faster than degradation
k9=1e-2;                %Toggle
k10=2.9e-2;             %Bartlet - kdeginna
k11=4.4e-5*3600;        %Literature (weird N-terminus rule)
kon=10000;              %Toggle
koff=1;                 %Toggle
Ve = 2e-4; % 1 * 10^-5; Barlett
Vi = 4e-12;             %Bartlett
Vp = 1.5e-3;            %Bartlett
kcleave=2;              %Literature

p=[k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, kon, koff, Ve, Vi, Vp, kcleave];

options=[];
[T,X] = ode15s(@New_model_equations,tspan,X0,options,p);

figure (1)
plot(T,X(:,1),T,X(:,4),T,X(:,7),T,X(:,10));
xlabel('Time (h)'); ylabel('# of molecules per L');
legend('Bloodstream SBP-siRNA','ECM SBP-siRNA','Endosomal SBP-siRNA','Cytoplasmic SBP-siRNA');

figure (2)
plot(T,X(:,4));
xlabel('Time (h)'); ylabel('# of molecules per L');
legend('ECM SBP-siRNA');

figure (3)
plot(T,X(:,7));
xlabel('Time (h)'); ylabel('# of molecules per L');
legend('Endosomal SBP-siRNA');

figure (4)
plot(T,X(:,12));

figure (5)
plot(T,X(:,13));


