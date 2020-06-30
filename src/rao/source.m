clear all
close all
clc

%% FLAGS



%% Define necessary parameters


Rt=0.025;
k=1.23;%ratio of specific heats
Mm=28;%kg/kmol, molar mass
Pc=1500;%psi, chamber pressure
Tc=1100;%K, guesstamate
Pa=14.7;%psi, SL pressure, initial design point
PR=Pc/Pa;%pressure ratio


TR=PR^((k-1)/k);%temperature ratio overall
TRt=1+(k-1)/2;%TR at throat
PRt=(1+(k-1)/2)^(k/(k-1));%PRt
M=sqrt(2/(k-1)*(TR-1));%Mach number


AR=((k+1)/2)^(-(k+1)/(2*(k-1)))*(TR)^((k+1)/(2*(k-1)))/M;%Area ratio;


%% Calculate throat area

Tt=Tc*TRt;%K, temp at throat
Pt=Pc*PRt;%Pa, throat pressure
At=pi*(Rt^2);

Re=sqrt(AR)*Rt;                                       %m, exit radius

%Calculate length of nozzle
K=0.8;                                                % percent length of same AR conical nozzle
Ln=K*(sqrt(AR)-1)*Rt/tan(15*pi/180);                  %m,length of nozzle


%% begin defining points

th_1=15*pi/180;%rad, angle at first circle
th_N=40*pi/180;%rad, angle b/w throat exit and parabola

x1=-1.5*Rt*sin(th_1);%x point of circle entering the nozzle
xN=0.382*Rt*sin(th_N);%X point at transition to parabola

%% parabola coefficients

RN=-sqrt((0.382*Rt)^2-xN^2)+1.382*Rt;%radius at xN


tmp1=[2*RN, 1, 0;
RN^2, RN, 1;
Re^2, Re, 1];

tmp2=[1/tan(th_N); xN; Ln];
dd=tmp1\tmp2;

a=dd(1);
b=dd(2);
c=dd(3);


x_c1=linspace(x1,0,50);%note, in cm
x_c2=linspace(0,xN,50);%circle exiting nozzle
y3=linspace(RN,Re,100);

th_e=1/(2*a*Re+b)*180/pi%deg, flow exit angle

ufo=size(x_c1,2); ufa=size(x_c2,2); %ufe=size(x_c3,2);

y1=-((1.5*Rt)^2*ones(1,ufo)-x_c1.^2).^.5+2.5*Rt*ones(1,ufo);
y2=-((0.382*Rt)^2*ones(1,ufa)-x_c2.^2).^.5+1.382*Rt*ones(1,ufa);

x_c3=a*y3.^2+b*y3+c*ones(1,size(y3,2));

x=[x_c1,x_c2,x_c3];
y=[y1,y2,y3];

z=ones(length(x));
plot(x,y,'r')
hold on
plot(x,-y,'r')

%axes('square');
% hold on
% plot(x_c2,y2);
% plot(x_c3,y3);
zz=[x;y;zeros(1,length(x))]';
% fid=fopen('nozzle.txt','w');
% fprintf(fid,'%6.2f %12.8f\n',zz);
% fclose(fid);
% xlswrite('nozzle.xls',zz);
