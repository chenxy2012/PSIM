function v_out = getOtfVal_m( otf, band, r_in, att)
if ~otf.isMultiBand
    band = 1;
end

v_out = zeros(size(r_in));
% v_out(r_in>otf.cutOff) = 0;

pos = r_in / otf.cyclesPerMicron + 1;
lPos = floor( pos );
hPos = ceil( pos );
f = pos - lPos;


for i = 1: 1: size(r_in,1)
    for j = 1: 1: size(r_in,2)
        if r_in(i,j) <= otf.cutOff &&  ceil(pos(i,j)) <= otf.samplesLateral
            l = otf.vals(band,lPos(i,j))*(1-f(i,j));
            h = otf.vals(band,lPos(i,j))*f(i,j);
            
            attl = 	otf.valsOnlyAtt(band,lPos(i,j)) * (1-f(i,j));
            atth = 	otf.valsOnlyAtt(band,lPos(i,j)) * f(i,j);
            if (att)
                v_out(i,j) = (l + h)^2*(attl + atth);
            else
                v_out(i,j) = (l + h)^2;
            end
        end
    end
end

end