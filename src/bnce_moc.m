% BNCE Method of Characteristics

function [xpoints, ypoints] = bnce_moc(TR, PR, T_1, g, R, Me)

PR2 = (PR)^((g-1)/g);
TT = (2*g*R*T_1)/(g-1);

v_e = sqrt(TT*(1-PR2));

T_e = T_1*(PR)^((g-1)/g);
a_e = sqrt(g*R*T_e);

Me = v_e/a_e;
    

RTOD = 180/pi;
DTOR = pi/180;
P = []; %x axis points
%% PM FUNCTION
A = sqrt((g+1)/(g-1));
B = (g-1)/(g+1);
v_PM = @(x) A*atan(sqrt(B*(x^2-1))) - atan(sqrt(x^2-1));


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
    %plot(P2,P1,'k')
    %hold on
    %xlabel('CENTERLINE')
    %ylabel('RADIUS')
end
%hold on;
LR(1) = []; RR(1) = [];
SL(1) = [];
F = RR(m-1);

for c = 1:length(P)-1
    x(c) = (TR+SL(c)*P(c))/(SL(c)-F);
    y(c) = F*x(c)+TR;
    X_P = [P(c) x(c)];
    Y_P = [0 y(c)];

end


TM = T_max*DTOR;
xw(1) = (TR+SL(1)*P(1))/(SL(1)-tan(TM));
yw(1) = tan(TM)*xw(1)+TR;
X_P2 = [P(1) xw];
Y_P2 = [P(2) yw];

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
        %plot(X_P3,Y_P3,'r');
    end
    %hold on
    
%    LAST POINT
    xf = (b(length(b))+SL(length(SL))*P(length(P)))/SL(length(SL));
    yf = b(length(b));
   % xw(length(xw)) = xf;
   % yw(length(yw)) = yf;
    
%     X_F = [P(length(P)) xf];
%     Y_F = [0 yf];
    %plot(X_F,Y_F,'r');
    
     xw = [0 xw];
     yw = [TR yw];
    
    xpoints = xw;
    ypoints = yw;
    
end

