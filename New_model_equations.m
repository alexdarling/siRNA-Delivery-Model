function dX = New_model_equations(t,X,p)

%%%Declare variables

%Blood
BPS=X(1);
BP=X(2);
BS=X(3);

%Interstitial
IPS=X(4);
IP=X(5);
IS=X(6);

%Endosome
EPS=X(7);
EP=X(8);
ES=X(9);

%Cytoplasm
CPS=X(10);
CP=X(11);
CS=X(12);

CSC=X(13); %cytoplasmic sirna cleaved

%%%Declare parameters
k1=p(1);                 %siRNA clearance from blood
k2=p(2);                 %protein degradation in blood
k3=p(3);                 %protein transport out of blood
k4=p(4);                 %siRNA degradation in ECM
k5=p(5);                 %protein degradation in ECM
k6=p(6);                 %endocytosis rate
k7=p(7);                 %siRNA degradation rate in endosome
k8=p(8);                 %protein degradation rate in endosome
k9=p(9);                 %endosomal escape rate
k10=p(10);               %siRNA degradation rate in cytoplasm
k11=p(11);               %protein degradation rate cytoplasm
kon=p(12);
koff=p(13);
Ve=p(14);
Vi=p(15);
Vp=p(16);
kcleave=p(17);           %rate of Dicer cleavage

%%%Equations
dBPS = kon*BP*BS -koff*BPS -k2*BPS -k3*BPS;
dBP = -kon*BP*BS +koff*BPS -k2*BP -k3*BP; 
dBS = -kon*BP*BS +koff*BPS -k1*BS;

dIPS = kon*IP*IS -koff*IPS -k5*IPS -k6*IPS +k3*BPS*Vp/Ve;
dIP = -kon*IP*IS +koff*IPS -k5*IP -k6*IP +k3*BP*Vp/Ve;
dIS = -kon*IP*IS +koff*IPS -k4*IS;

dEPS = kon*EP*ES -koff*EPS -k8*EPS -k9*EPS +k6*IPS*Ve/Vi;
dEP = -kon*EP*ES +koff*EPS -k8*EP -k9*EP +k6*IP*Ve/Vi;
dES = -kon*EP*ES +koff*EPS -k7*ES -k9*ES;

dCPS = kon*CP*CS -koff*CPS -k11*CPS +k9*EPS; 
dCP = -kon*CP*CS +koff*CPS -k11*CP +k9*EP;
dCS = -kon*CP*CS +koff*CPS -k10*CS +k9*ES -kcleave*CS;

dCSC = kcleave*CS;

dX=[dBPS; dBP; dBS; dIPS; dIP; dIS; dEPS; dEP; dES; dCPS; dCP; dCS; dCSC];

