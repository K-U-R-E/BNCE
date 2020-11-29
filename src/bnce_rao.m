% BNCE RAO

function [xpoints, ypoints] = bnce_rao(G, RT,s,PR_Ratio)

disp("Running Rao Solver");

PR=PR_Ratio^-1;%pressure ratio


TR=PR^((G-1)/G);%temperature ratio overall

M=sqrt(2/(G-1)*(TR-1));%Mach number


AR=((G+1)/2)^(-(G+1)/(2*(G-1)))*(TR)^((G+1)/(2*(G-1)))/M;%Area ratio;


Re=sqrt(AR)*RT                                       %m, exit radius

%Calculate length of nozzle
G=s;                                                % percent length of same AR conical nozzle
Ln=G*(sqrt(AR)-1)*RT/tan(15*pi/180);                  %m,length of nozzle


th_1=15*pi/180;%rad, angle at first circle
th_N=40*pi/180;%rad, angle b/w throat exit and parabola

x1=-1.5*RT*sin(th_1);%x point of circle entering the nozzle
xN=0.382*RT*sin(th_N);%X point at transition to parabola


RN=-sqrt((0.382*RT)^2-xN^2)+1.382*RT;%radius at xN


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

ufo=size(x_c1,2); ufa=size(x_c2,2); %ufe=size(x_c3,2);

y1=-((1.5*RT)^2*ones(1,ufo)-x_c1.^2).^.5+2.5*RT*ones(1,ufo);
y2=-((0.382*RT)^2*ones(1,ufa)-x_c2.^2).^.5+1.382*RT*ones(1,ufa);

x_c3=a*y3.^2+b*y3+c*ones(1,size(y3,2));

x=[x_c1,x_c2,x_c3];
y=[y1,y2,y3];

xpoints = x;
ypoints = y;

end
