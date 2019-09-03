% Generating OTF information based on SIM parameters
% Matlab code of a fairSIM function. https://github.com/fairSIM
function  ret =  otfGenerator( NA, lambda, a, attStrength, attFWHM)
ret.samplesLateral = 512;
ret.NA = NA; 
ret.lambda = lambda;
ret.cutOff = 1000 / ( lambda / NA / 2 ); 
ret.cyclesPerMicron = ret.cutOff / ret.samplesLateral ;

ret.vals	= zeros(1,ret.samplesLateral);
ret.valsAtt	= zeros(1,ret.samplesLateral);
ret.valsOnlyAtt = zeros(1,ret.samplesLateral);
ret.isMultiBand = false;
ret.isEstimate  = true;
ret.estimateAValue = a;

for i = 1: 1: ret.samplesLateral
    v = (i-1)/ ret.samplesLateral ;
    valIdealOTF = (2/pi)*(acos(v) - v*sqrt(1-v*v));
    r = valIdealOTF * ( a^v );
    ret.vals(1,i) = r;
end


ret.attStrength = attStrength;
ret.attFWHM = attFWHM;

for v = 1: 1: length(ret.vals)
    dist = (v-1) *  ret.cyclesPerMicron;
    vv = 1 - attStrength * exp( -dist^2  / (2* ( attFWHM/2.355)^2) );
    ret.valsOnlyAtt(1,v) = vv ;
    ret.valsAtt(1,v) = ret.vals(1,v) * ret.valsOnlyAtt(1,v);
end

end