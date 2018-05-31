function outV = pasteFreq(inV)
wi = size(inV,2);
hi = size(inV,1);

wo = 2*wi;
ho = 2*hi;

outV  = zeros(ho,wo);


outV( 1:wi/2, 1:wi/2) = inV( 1:wi/2, 1:wi/2);
outV( 1:wi/2, wo/2+(wi/2+1:wi)) = inV( 1:wi/2, wi/2+1:wi);
outV( wo/2+(wi/2+1:wi), 1:wi/2) = inV( wi/2+1:wi, 1:wi/2);
outV( wo/2+(wi/2+1:wi), wo/2+(wi/2+1:wi)) = inV( wi/2+1:wi, wi/2+1:wi);