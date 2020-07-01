function [surfpoints] = bnce_spin(xpoints, ypoints)
surfpoints = [3, length(ypoints), 360]; 
    for i = 1:length(ypoints)
        x = xpoints(i);
        Y = ypoints(i);
        disp(x);
        r = Y;
        for ii = 1:360
            Y = r * cos(ii*(pi/180));
            z = r * sin(ii*(pi/180));
            surfpoints(1, i, ii) = z;
            surfpoints(2, i, ii) = Y;
            surfpoints(3, i, ii) = x;
        end
    end  
end
