function im = writeAttenuationVector( w, h, str, fwhm,CyclesPerMicron)

im = zeros( h, w);
for x = 1: 1: w
    for y = 1: 1: h
        if x <= w/2
            xh = x-1;
        else
            xh = x-w-1;
        end
        if y <= h/2
            yh = y-1;
        else
            yh = y-h-1;
        end
                
        rad = sqrt(xh^2+yh^2);
        cycl = rad * CyclesPerMicron;
        v = (1 - str * exp( -cycl^2)) / (2* (fwhm/2.355)^2);
        im(y,x) = v;    
    end
end

end
