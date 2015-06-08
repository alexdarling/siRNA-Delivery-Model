clear all; close all; clc

%%%Initial values
kon=1e-2;              %Toggle 
koff=1e-4;             %Toggle
BPS0=1e-9*6e23;
BP0=sqrt(BPS0/6e23 * koff/kon)*6e23;
BS0=BP0;
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

tottime = 24;
day = 7;
tspan=[0,day*tottime];

k1=600*60;              %Literature x
k2=3.82e-7*3600;        %Literature x
k3=1e-2*1e-3/1000;      %Bartlett - ktransblood*partition
k4=0.09*60;             %Literature x
k5=8.7e-2;              %Bartlett x
k6=1.6e-5*3600;         %Literature x
k7=5e-1;                %Bartlett - kdegenda
k8=0;                   %Assume endosmal escape is much faster than degradation
k9=1;                %Toggle
k10=2.9e-2;             %Bartlett - kdeginna x
k11=4.4e-5*3600;        %Literature (weird N-terminus rule) x
Ve = 2e-4; % 1 * 10^-5; %Bartlett x
Vi = 4e-11;             %Bartlett x [[11 or 12?]]
Vp = 1.5e-3;            %Bartlett x
kcleave=2;              %Literature

p=[k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, kon, koff, Ve, Vi, Vp, kcleave];
options = [];
[T,X] = ode15s(@New_model_equations,tspan,X0,options,p);

figure (1)
title('Blood Stream')
subplot(1,3,1)
plot(T,X(:,1))
subplot(1,3,2)
plot(T(1:90),X(1:90,2))
subplot(1,3,3)
plot(T(1:90),X(1:90,3))

figure (2)
title('ECM')
subplot(1,3,1)
plot(T,X(:,4))
subplot(1,3,2)
plot(T,X(:,5))
subplot(1,3,3)
plot(T,X(:,6))

figure (3)
title('Endosome')
subplot(1,3,1)
plot(T,X(:,7))
subplot(1,3,2)
plot(T,X(:,8))
subplot(1,3,3)
plot(T,X(:,9))

figure (4)
title('Cell Cytoplasm')
subplot(1,3,1)
plot(T,X(:,10))
subplot(1,3,2)
plot(T,X(:,11))
subplot(1,3,3)
plot(T,X(:,12))

figure (5)
plot(T, X(:,13))
title('Cleaved SiRNA')


