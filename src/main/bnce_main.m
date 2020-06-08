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
fprintf("Bell Nozzle Contour Evaluator")
fprintf("Written by Mohammed Zweiri, Pritom Chowdhury and Vinay Williams")
fprintf("at Kingston University, 2020")

%% Input variables:
mdot = 1;
G = 1.475;                                                                                           % Specific heat ratio
theta_1st = 13;                                                                                      % half cone angle, 15 degrees is often used a reference in comparing lengths and thrust performances(Degree)
theta_nth = 40;                                                                                      % Angle that seperates the parabola and throat exit(Degree)
RT = 0.025;                                                                                          % Diameter of the throat(m)
TC = 3500;                                                                                           % Chamber temperature(K)
PA = 101325;                                                                                         % Ambient Pressure(Pa)
PC = 4882399.71812414;                                                                               % Chamber pressure(Pa)
k = 0.8;                                                                                             % 80% length has expansion ratio as the 100%, hence the efficiency between them is only 0.2% difference.
PR_Ratio = PA/PC                                                                                     % Pressure ratio
G1 = G + 1;                                                                                          % Constant
G2 = G - 1;                                                                                          % Constant
TRatio = PR_Ratio^(G2/G)                                                                             % Temperature Ratio

%% Start of calculations


PR_throat = (1+G2/2)^(1/G2)                                                                          % Pressure ratio at the throat
TR_throat = 1+(G2/2)                                                                                 % Temperature ratio at the throat
T = TC * TR_throat                                                                                   % Temperature at the throat(K)
P = PC * PR_throat                                                                                   % Pressure at the throat(Pa)
AT = pi*(RT^2)                                                                                       % Area of the throat(m^2)

e = ((G2/2)^(1/2))*((2/G1)^((G1)/(2*G2)))*((PR_Ratio)^(-1/G))*(1-((PR_Ratio)^(G2/G)))^(-0.5)         % Expansion Ratio

M = sqrt((TRatio-1)*(2/G2))                                                                          % Exit Mach Number
Ae = e*AT                                                                                            % Exit Area (m^2)
Re = sqrt(e*AT)                                                                                      % Exit Radius(m)

L  = (k*(sqrt(e)-1)*RT)/tan(theta_1st*pi/180)                                                        % Length of the nozzle(m)


%A TOP nozzle, using Rao coefficients to define the circular curves entering and exiting the throat,
%equal to 1.5RT and 0.382RT



Angle_1 = 15*pi/180;                                                                                 % Initial Angle of Converging Region wrt to throat (rad)
Angle_n = 40*pi/180;                                                                                 % Initial Angle of Diverging Region wrt to throat (rad)

x_1 = -1.5*RT*sin(theta_1st)                                                                         % x position of first point in converging region
x_nth = 0.382*RT*sin(theta_nth)                                                                      % x position of first point in diverging region

y_n = -sqrt((0.382*RT)^2-(x_nth)^2)+1.382*RT                                                         % The radius at x_nth (m)

% Full system of equations for the parabolic coefficients for the first parabola
Matrix_1 = [2*y_n,1,0;2*Re,1,0;(y_n)^2,y_n,1]
Matrix_2 = [(tan(Angle_n))^-1;(tan(Angle_1))^-1;x_nth]
Matrix_3 = Matrix_1\Matrix_2

a = Matrix_3(1)
b = Matrix_3(2)
c = Matrix_3(3)

N1 = linspace(x_1,0,100);
N2 = linspace(0,x_nth,100);
W1 = linspace(y_n,Re,400);

Angle_Exit = tan((1/((2*a*Re)+b))*180/pi)

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
plot(x, -y, 'b')
hold on
%plot(y,'b')
xlabel('x/m')
ylabel('y/m')
grid on
title('Nozzle Contour')

%% Performance Characteristics 
te = T* (1 + (G2/2)* M^2);
ve = M * sqrt(G*8.314*te);


if mdot == 0
    f = (PA-PC)*Ae;
else
    f = mdot*ve +(PA-PC)*Ae;
end
disp("Exit Mach Number: "M)
disp(f)
