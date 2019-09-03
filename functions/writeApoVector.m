% Generating apo matrix from the given system parameters
% Matlab code of a fairSIM function. https://github.com/fairSIM
function vec = writeApoVector( otf, bend, cutOff, h, w)

xx = 1: 1: w;
yy = 1: 1: h;
[x,y] = meshgrid(xx,yy);

x(:,1:w/2) = x(:,1:w/2)-1;
x(:,w/2+1:w) = x(:,w/2+1:w) - w - 1;

y(1:h/2,:) = -(y(1:h/2,:)-1);
y(h/2+1:h,:) = h - (y(h/2+1:h,:)-1);

cycl = sqrt(x.^2 + y.^2) * otf.vecCyclesPerMicron;
frac = cycl / (otf.cutOff*cutOff);

frac = (2/pi)* (acos(frac)-frac.*sqrt(1-frac.^2));

vec = frac.^bend;

end

