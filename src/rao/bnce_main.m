% MIT License
%
% Copyright (c) 2020 Mohammed Zweiri, Pritom Chowdhury and Vinay Williams
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

%% References

% Rao, G.V.R., 1958. Exhaust nozzle contour for optimum thrust. *Journal of Jet Propulsion*, *28*(6), pp.377-382.

% Davis, K., Fortner, E., Heard, M., McCallum, H. and Putzke, H., 2015. Experimental and
% computational investigation of a dual-bell nozzle. In *53rd AIAA Aerospace Sciences Meeting* (p. 0377).



%%
clear all
close all
clc
disp("Bell Nozzle Contour Evaluator")
disp("Written by Mohammed Zweiri, Pritom Chowdhury and Vinay Williams")
disp("at Kingston University, 2020")

%% Input variables:
mdot = 1;
G = 1.475;                                                                                           % Specific heat ratio
                                                                                    % Angle that seperates the parabola and throat exit(Degree)
RT = 0.025;                                                                                          % Diameter of the throat(m)
TC = 3500;                                                                                           % Chamber temperature(K)
PA = 101325;                                                                                         % Ambient Pressure(Pa)
PC = 4882399.71812414;                                                                               % Chamber pressure(Pa)
k = 1;                                                                                             % 80% length has expansion ratio as the 100%, hence the efficiency between them is only 0.2% difference.
PR_Ratio = PA/PC;                                                                                    % Pressure ratio
G1 = G + 1;                                                                                          % Constant
G2 = G - 1;                                                                                          % Constant
                                                                           % Temperature Ratio

%% Start of calculations

TRatio = PR_Ratio^(G2/G);
PR_throat = 1+(G2/2)^(1/G2);                                                                          % Pressure ratio at the throat
TR_throat = 1+(G2/2);                                                                                 % Temperature ratio at the throat
T = TC * TR_throat;                                                                                   % Temperature at the throat(K)
P = PC * PR_throat;                                                                          % Pressure at the throat(Pa)
M = sqrt((TRatio-1)*(2/G2));

AT = pi*(RT^2);                                                                                       % Area of the throat(m^2)
e = ((G2/2)^(1/2))*((2/G1)^((G1)/(2*G2)))*((TRatio)^(-1/G))*(1-((TRatio)^(G2/G)))^(-0.5);         % Expansion Ratio

                                                                          % Exit Mach Number
Ae = e*AT;                                                                                            % Exit Area (m^2)
Re = sqrt(e*AT);                                                                                      % Exit Radius(m)

L  = (k*(sqrt(e)-1)*RT)/tan(15*pi/180);                                                        % Length of the nozzle(m)



Angle_1 = 15*pi/180;                                                                                 % Initial Angle of Converging Region wrt to throat (rad)
Angle_n = 40*pi/180;                                                                                 % Initial Angle of Diverging Region wrt to throat (rad)

x_1 = -1.5*RT*sin(Angle_1);                                                                         % x position of first point in converging region
x_nth = 0.382*RT*sin(Angle_n);                                                                      % x position of first point in diverging region

y_n = -sqrt((0.382*RT)^2-(x_nth)^2)+1.382*RT


                                                    % The radius at x_nth (m)

% Full system of equations for the parabolic coefficients for the first parabola
Matrix_1 = [2*y_n,1,0;2*Re,1,0;(y_n)^2,y_n,1];
Matrix_2 = [(tan(Angle_n))^-1;(tan(Angle_1))^-1;x_nth];
Matrix_3 = Matrix_1\Matrix_2;

a = Matrix_3(1);
b = Matrix_3(2);
c = Matrix_3(3);

N1 = linspace(x_1,0,100);
N2 = linspace(0,x_nth,100);
W1 = linspace(y_n,Re,400);

Angle_Exit = tan((1/((2*a*Re)+b))*180/pi);

UWO = size(N1,2);UXR = size(N2,2);                                                                   % Creating array that return each size of the dimension

y_1 = -sqrt((1.5*RT)^2*(ones(1,UWO))-N1.^2)+2.5*RT*ones(1,UWO);
y_2 = -sqrt((0.382*RT)^2*ones(1,UXR)-N2.^2)+1.382*RT*ones(1,UXR);

N3 = a*W1.^2+b*W1+c*ones(1,size(W1,2));

%% Plot the curve

x=[N1,N2,N3];
y=[y_1,y_2,W1];

z=zeros(length(x));
hold on
plot(x,y,'r')
plot(x, -y, 'r')
hold on
%plot(y,'b')
xlabel('x/m')
ylabel('y/m')
grid on
title('Nozzle Contour')

%-----------------------------------------------------------------------------------
%% FLAGS



%% Define necessary parameters


Rt=0.014;
k=1.475;%ratio of specific heats


Pc=4882399.71812414;%psi, chamber pressure
Tc=3500;%K, guesstamate
Pa=101325;%psi, SL pressure, initial design point
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
%plot(x,y,'g')
plot(x,-y,'g')

%axes('square');
% hold on
% plot(x_c2,y2);
% plot(x_c3,y3);
zz=[x;y;zeros(1,length(x))]';
% fid=fopen('nozzle.txt','w');
% fprintf(fid,'%6.2f %12.8f\n',zz);
% fclose(fid);
% xlswrite('nozzle.xls',zz);

%INPUT VALUES
p_1 = 4882399.71812414;     %CHAMBER PRESSURE
T_1 = 3500;     %CHAMBER TEMP
FT = 0;      %DESIRED THRUST OR....
m_dot = 1.875;   %DESIRED MASS FLOW RATE....
g = 1.475;       %GAMMA
R = 0.25;       %GAS CONSTANT
p_o = 101325;

%%  begin calculation
PR = p_o/p_1;
PR2 = (p_o/p_1)^((g-1)/g);
TT = (2*g*R*T_1)/(g-1);
p_t = ((2/(g+1))^(g/(g-1)))*2.068;
v_t = sqrt((2*g*R*T_1)/(g+1));
v_e = sqrt(TT*(1-PR2));

if m_dot==0
    m_dot=FT/v_e;
elseif FT==0
    FT = m_dot/v_e;
else
    fprintf('You can either set desired thrust OR mass flow rate')
end

T_e = T_1*(p_o/p_1)^((g-1)/g);
a_e = sqrt(g*R*T_e);

Me = v_e/a_e;

% MOC
TR = 0.025; %throat radius (cm)
RTOD = 180/pi;
DTOR = pi/180;
P = []; %x axis points
%% PM FUNCTION
v_PM = @(x) (sqrt((g+1)/(g-1)))*atan(sqrt(((g-1)/(g+1))*(x^2-1))) - atan(sqrt(x^2-1));


%% CALCULATE T_MAX, BREAK UP INTO DIVISIONS
T_max = 0.5*v_PM(Me)*RTOD;
DT = (90-T_max) - fix(90-T_max);
T(1) = DT*DTOR;
n = T_max*2;

for m = 2:n+1
    T(m) = (DT + (m-1))*DTOR;
    %Mach from T(i) using T(i) = v_PM (FALSE POSITION)
    x_int = [1 1.01*Me];
    func = @(x) T(m) - v_PM(x);
    M(m) = fzero(func,x_int);
    P(m) = 0 + TR*tan(T(m)); %X-AXIS POINTS
    %RRSLOPES
    RR(m) = -TR/P(m);
    %LR slopes
    LR(m) = tan(T(m)+asin(1/M(m)));
    SL(m) = -RR(m);
end
%% PLOTTING
P(1) = [];
l = length(P);

for j = 1:l
    P1 = [0 TR];
    P2 = [P(j) 0];
    plot(P2,P1,'k')
    hold on
    xlabel('CENTERLINE')
    ylabel('RADIUS')
end
hold on;
LR(1) = []; RR(1) = [];
SL(1) = [];
F = RR(m-1);

for c = 1:length(P)-1
    x(c) = (TR+SL(c)*P(c))/(SL(c)-F);
    y(c) = F*x(c)+TR;
    X_P = [P(c) x(c)];
    Y_P = [0 y(c)];
    plot(X_P,Y_P,'b');
end
hold on

%% FIRST WALL SECTION
TM = T_max*DTOR;
xw(1) = (TR+SL(1)*P(1))/(SL(1)-tan(TM));
yw(1) = tan(TM)*xw(1)+TR;
X_P2 = [P(1) xw];
Y_P2 = [P(2) yw];
plot(X_P2,Y_P2,'g');
%DIVIDE (delta slopes)
DTW = tan(TM)/(length(P)-1);
s(1) = tan(TM);
b(1) = TR;

    for k = 2:length(P)-1
        s(k) = tan(TM)-(k-1)*DTW; %slope
        b(k) = yw(k-1)-s(k)*xw(k-1); %y-int
        xw(k) = (b(k)+SL(k)*P(k))/(SL(k)-s(k));
        yw(k) = s(k)*xw(k)+b(k);
        X_P3 = [x(k) xw(k)];
        Y_P3 = [y(k) yw(k)];
        plot(X_P3,Y_P3,'r');
    end
    hold on

   % LAST POINT
    xf = (b(length(b))+SL(length(SL))*P(length(P)))/SL(length(SL));
    yf = b(length(b));
    X_F = [P(length(P)) xf];
    Y_F = [0 yf];
    plot(X_F,Y_F,'r');

    xw = [0 xw];
    yw = [TR yw];
    RTHROAT = TR;
    REXIT = yw(length(yw));

    AR = (RTHROAT/REXIT)^2

    %xw
    %yw
