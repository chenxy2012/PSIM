function vec_out = maskOtf( otf, vec, kx, ky)
w = size(vec,2);
h = size(vec,1);
mask = zeros(h,w);
xx = 1: 1: w;
yy = 1: 1: h;
[x,y] = meshgrid(xx,yy);

x(:,1:w/2) = x(:,1:w/2)-1;
x(:,w/2+1:w) = x(:,w/2+1:w) - w - 1;

y(1:h/2,:) = -(y(1:h/2,:)-1);
y(h/2+1:h,:) = h - (y(h/2+1:h,:)-1);

rad = sqrt( (x-kx).^2 + (y-ky).^2 ) * otf.vecCyclesPerMicron;
mask(rad<=otf.cutOff) = 1; 
vec_out = vec.*mask;
end