function [band0_out,band1_out] = commonRegion( band0, band1, bn0, bn1, otf, kx, ky, dist, weightLimit, divideByOtf ) 
w = size(band0,2);
h = size(band0,1);
band0_out = zeros(h,w);

weight0 = zeros(h,w);
weight1 = zeros(h,w);
wt0 = zeros(h,w);
wt1 = zeros(h,w);

writeOtfVector( otf, weight0, bn0, 0, 0);
writeOtfVector( otf, weight1, bn1, 0, 0);
writeOtfVector( otf, wt0, bn0, kx, ky);
writeOtfVector( otf, wt1, bn1, -kx, -ky);

cutCount =0;

for y = 1: 1: h
    for x = 1: 1: w
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
        maxv = Math.sqrt( kx*kx+ ky*ky ); 

	    ratio = rad / maxv;
        
        if  abs(weight0(y,x)) < weightLimit || abs(wt0(y,x)) < weightLimit
            cutCount = cutCount + 1;
            band0_out(y,x)=0;
        else
            if divideByOtf
                band0_out(y,x) = band0_out(y,x) / weight0(y,x);
            end
        end
        
        if  abs(weight1(y,x)) < weightLimit || abs(wt1(y,x)) < weightLimit
            band1_out(y,x)=0;
        else
            if divideByOtf
                band1_out(y,x) = band1_out(y,x) / weight1(y,x);
            end
        end
        
        if ratio<dist || ratio>(1-dist) 
            band0_out(y,x) = 0;
            band1_out(y+round(ky)+h,x-round(kx)+w) = 0;
            cutCount = cutCount + 1;
        end

    end
end

end