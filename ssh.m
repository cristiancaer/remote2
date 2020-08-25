clc
clear all
%% 
% Variables de estado
TS=491;% main steam temperature INITIAL
TM=501;% SSH temperature INITIAL
%%
% Punto de operación
TSi=450;
TGi=1028;
%%
% Entradas
tsi=0;% SSH inlet steam temperature
tg=0;% furnance gas temperature
%%
% Constantes
ai=2210.25;% inside heat transfer area of SSH
ao=2454;% outside heat transfer area
alfams=3616;% metal to steam heat transfer coefficient
alfagm=58.9;% gas to metal heat transfer coefficient
Cv=0.4634;% specific heat of SSH steam at constan volume
Cs=0.4634;% specific heat of main steam 
Csi=0.4634;% especific heat of main steam inlet
Cm=0.15;% especific heat of SSH tube
Fs=690e3;% mass flow rate of stemam in SSH
Mm=252e3;% mass of SSH metal
rs=50.5;% specific wight of steam in SSH
Vs=23;% control volumen of SSH
%%
% Matrices de estados
A1=[-(ai*alfams+Fs*Cs)/(Vs*rs*Cv),ai*alfams/Vs*rs*Cv;ai*alfams/Mm*Cm,-(ao*alfagm+ai*alfams)/(Mm*Cm)]
B1=[Fs*Csi/(Vs*rs*Cv), 0;0, ao*alfagm/(Mm*Cm)]
C=[1,0;0,1];
D=[0,0;0,0];
A=[-3.94751 3.782493;0.053859 -0.054981];
B=[0.165017 0; 0 0.001122];
%%
% Modelo
modelo=ss(A,B,C,D)

%%
% funcion de transferencia con entrada TSi
[num, den]=ss2tf(A,B,C,D,1);
Ts=tf(num(1,:),den);
Tm=tf(num(2,:),den);
%step(Ts)
%%
% funcion de transferencia con entrada Tg
[num, den]=ss2tf(A,B,C,D,2);
Ts=tf(num(1,:),den);
Tm=tf(num(2,:),den);
%step(Ts)
%%
% TIEMPO de muestreo
T=20;
t=[0:T:1e3];
%% 
% codiciones iniciales
X0=[TS,TM];
u0=[TSi,TGi];
t=0:20:4000;
sysSSH=ss(A,B,C,[]);
Tsi=step_signal(t,22, TSi, -2000)';
Tg=step_signal(t,-0*51,TGi,-2000)';
u=[Tsi Tg];
[y1, t1, x1]=lsim(sysSSH,u,t,X0);
figure(2)
subplot(2,1,1)
plot(t1/20,y1);ylim([490 680]);%se divide en 20 para ver en # sample
grid on; grid minor;
subplot(2,1,2)
plotyy(t/20,u(:,1),t/20,u(:,2))
grid on; grid minor;