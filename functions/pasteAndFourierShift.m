% Doubleing the Fourier domain and paste the original part in the shifted
% region
% Matlab code of a fairSIM function. https://github.com/fairSIM
function outV = pasteAndFourierShift( inV, kx, ky)
outV0 = pasteFreq( inV);
if_outV = ifft2(outV0);

shitf_if_outV = fourierShift(if_outV, kx, -ky );
outV = fft2(shitf_if_outV);
end