% BNCE MAIN

standalone = 1;

clc;
close all;
clear;

disp("Getting values and doing preliminary calculations")

solver = xlsread('bnce_inputs.xlsx','bnce_inputs','B3');

if solver == 1 || solver == 0
    % Read the percentage length variable for the Rao solver
    if solver == 0
        fraction_length = xlsread('bnce_inputs.xlsx','bnce_inputs','F3');
    end
    disp("Solver: Good")
else
    disp("Check solver choice")
    return
end


gamma = xlsread('bnce_inputs.xlsx','bnce_inputs','B5');
R = xlsread('bnce_inputs.xlsx','bnce_inputs','B6');

% Check the flag for defined ambient pressure or defined altitude and does
% the relevant calculations
if xlsread('bnce_inputs.xlsx','bnce_inputs','B8') == 1
    pressure_ambient = xlsread('bnce_inputs.xlsx','bnce_inputs','B9');
    disp("Ambient Pressure v Altitude: Good")
elseif xlsread('bnce_inputs.xlsx','bnce_inputs','B8') == 0
    altitude = xlsread('bnce_inputs.xlsx','bnce_inputs','B10');
    
    % Following if elseif else is NASA's atmospheric model
    if (11000>altitude) && (altitude<25000)
        temperature_ambient = -56.46; 
        pressure_ambient = 1000*(22.65*exp(1.73-0.000157*altitude));  
    elseif altitude >=25000
        temperature_ambient = -131.21 + 0.00299*altitude ;
        pressure_ambient = 1000*(2.488*((temperature_ambient+273.1)/216.6)^-11.388);   
    else 
        temperature_ambient = 15.04 - 0.00649*altitude;
        pressure_ambient = 1000*(101.29*((temperauture_ambient+273.1)/288.08)^5.256);
    end
    disp("Ambient Pressure v Altitude: Good")
else
    disp("Check the altitude v ambient pressure flag")
    return
end


pressure_chamber = xlsread('bnce_inputs.xlsx','bnce_inputs','B12');
temperature_chamber = xlsread('bnce_inputs.xlsx','bnce_inputs','B13');

% % Checks the flag for the explicit defintion of throat radius or its
% % calculation based on other parameters
% if xlsread('bnce_inputs.xlsx','bnce_inputs','B2') == 1
%     radius_throat = xlsread('bnce_inputs.xlsx','bnce_inputs','B2');
% else
%     % pass
% end

pressure_ratio = pressure_ambient/pressure_chamber;
pressure_ratio2 = pressure_ratio ^ ((gamma -1)/gamma); 
temperature_throat = (2*temperature_chamber)/(gamma+1); 

pressure_throat = ((2/(gamma+1))^(gamma/(gamma-1)))*2.068;
velocity_throat = sqrt((2*gamma*R*temperature_chamber)/(gamma+1));
velocity_exit = sqrt(temperature_throat*(1-pressure_ratio2));


% Checks the flag for calculating force from mass flow rate or mass flow
% rate from force and then performs the relevant math
if xlsread('bnce_inputs.xlsx','bnce_inputs','B15') == 1
    force = xlsread('bnce_inputs.xlsx','bnce_inputs','B16');
    mass_flow_rate = force / velocity_exit;
    disp("Force v Mass Flow Rate: Good")
elseif xlsread('bnce_inputs.xlsx','bnce_inputs','B15') == 0
    mass_flow_rate = xlsread('bnce_inputs.xlsx','bnce_inputs','B17');
    force = mass_flow_rate * velocity_exit;
    disp("Force v Mass Flow Rate: Good")
else
    disp("Check force v mass flow rate flag")
    return
end

temperature_exit = temperature_chamber*(pressure_ambient/pressure_chamber)^((gamma-1)/gamma);
area_exit = sqrt(gamma*R*temperature_exit);
mach_exit = velocity_exit / area_exit;

disp("Let the show begin");
radius_throat = 0.01;
if solver == 0
    cd src
    [xpoints, ypoints] = bnce_rao(gamma, radius_throat,fraction_length,pressure_ratio);
    plot(xpoints, ypoints);
    hold on 
    cd ..
else
    cd src
    [xpoints, ypoints] = bnce_moc(radius_throat, pressure_ratio, temperature_chamber, gamma, R, mach_exit);
    cd ..
    %pass
end

% @New Shrike Programmers: Unifinished below, I think the math for rotatating the curve has a problem, not quite sure what it is ~ Vinay

%cd src; surfpoints = bnce_spin(xpoints,ypoints);

%     hold on
%     surf2solid(surfpoints(1, :), surfpoints(2, :), surfpoints(3, :));
%     % Plot contour
%     title('3D Contour')
%     ylabel('y')
%     xlabel('Width (mm)')
%     zlabel('Legth (mm)')
%     axis equal
%     grid on
%     
%     x,y,z = size(surfpoints);
%     solid = surf2solid(surfpoints(1, 2:x-5), surfpoints(2, 2:y-5), surfpoints(3, 2:z-5));

    
cd ..;
disp("Let's get some stats")

% @ New Shrike Programmers: Unfinished here, you need to implement a stats page ~ Vinay









