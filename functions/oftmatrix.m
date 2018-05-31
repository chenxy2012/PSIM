function v_out = oftmatrix( otf, w, h)

val = zeros(h,w);

xx = 1: 1: w;
yy = 1: 1: h;
[x,y] = meshgrid(xx,yy);

x(:,1:w/2) = x(:,1:w/2)-1;
x(:,w/2+1:w) = x(:,w/2+1:w) - w - 1;

y(1:h/2,:) = -(y(1:h/2,:)-1);
y(h/2+1:h,:) = h - (y(h/2+1:h,:)-1);

rad = sqrt( (x).^2 + (y).^2 );
cycl = rad * otf.vecCyclesPerMicron;


pos = cycl / otf.cyclesPerMicron + 1;
lPos = floor( pos );
hPos = ceil( pos );
f = pos - lPos;


for i = 1: 1: size(cycl,1)
    for j = 1: 1: size(cycl,2)
        if cycl(i,j) <= otf.cutOff &&  ceil(pos(i,j)) <= otf.samplesLateral
            
            l = 	otf.vals(1,lPos(i,j)) * (1-f(i,j));
            h = 	otf.vals(1,lPos(i,j)) * f(i,j);
            val(i,j) = l + h;
        end
    end
end
v_out = val;
end