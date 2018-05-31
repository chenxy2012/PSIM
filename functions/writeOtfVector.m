function vec_out = writeOtfVector(otf, vec, band, kx, ky)
% false true
vec_out = zeros(size(vec));

w = size(vec,2);
h = size(vec,1);

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
        
        rad  = sqrt(xh^2 + yh^2);
        cycl = rad * otf.cyclesPerMicron;
        
        if cycl > otf.cutOff
            vec_out(y,x) = 0;
        end
        
        if cycl < otf.cutOff
            val = getOtfVal(otf, band, cycl);
            vec_out(y,x) = val;
        end
    end
end


end