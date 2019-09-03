% Fade the border of the input images
% Matlab code of a fairSIM function. https://github.com/fairSIM
function img_out = fadeBorderCos(img, px)
w = size(img,2);
h = size(img,1);

img_out = img;

fac = (1/px) *( pi/2);

for y = 1: 1: px
    for x = 1: 1: w
        img_out(y,x) = img_out(y,x) * sin( (y-1) * fac )^2;
    end
end;

for y = h-px+1: 1: h
    for x = 1: 1: w
        img_out(y,x) = img_out(y,x) * sin( (h-y) * fac )^2;
    end
end

for y=1:1:h
    for x = 1: 1: px
        img_out(y,x) = img_out(y,x) * sin( (x-1) * fac )^2;
    end
end

for y=1:1:h
    for x = w-px+1: 1: w
        img_out(y,x) = img_out(y,x) * sin( (w-x) * fac )^2;
    end
end

end