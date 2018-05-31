function outV = pasteAndFourierShift( inV, kx, ky)
outV0 = pasteFreq( inV);
if_outV = ifft2(outV0);

shitf_if_outV = fourierShift(if_outV, kx, -ky );
outV = fft2(shitf_if_outV);
end