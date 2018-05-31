function outV = fourierShift( inV, kx, ky )
% val = inV;
N = size(inV,1);

xx = 1: 1: N;
yy = 1: 1: N;
[y,x] = meshgrid(xx,yy);
phaVal = 2*pi*(x*ky/N+y*kx/N);


outV = inV .* exp(1i*phaVal);

% 
% for x = 1: 1: N
%     for y = 1: 1: N
%         phaVal = 2*pi*(kx*x+ky*y) / N;
% 
%         co = cos( phaVal );
%         si = sin( phaVal );
%                 
%         val(y,x) = val(y,x) * (co+1i*si); 
%     end 
% end
% outV = val;

end