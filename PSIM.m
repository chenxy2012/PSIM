%% The main program
% This program is developed for pSIM reconstruction based on the open source fairSIM plugin
% More functions are under development
% 
% For some reference:
% pSIM: Zhanghao, Karl, et al. "Super-resolution Imaging of Fluorescent Dipoles by Polarized Structured Illumination Microscopy." arXiv preprint arXiv:1712.05092 (2017).
% fairSIM: M¨¹ller, Marcel, et al. "Open-source image reconstruction of super-resolution structured illumination microscopy data in ImageJ." Nature communications 7 (2016): 10980.
% 
% website:
% pSIM:
% fairSIM: https://github.com/fairSIM
% 
% For any question, please contact:
% chenxy16@mails.tsinghua.edu.cn


%% Input your SIM parameters here
% For the SIM parameters, if the images is captured by commercial SIM syste there should be an accompanied parameter file. 
% It can also be estimated with the help of fairSIM plugin, which is
% available at: https://github.com/fairSIM


% We apply fairSIM as our SIM reconstruction algorithm. However, sometimes we do notice that the quality of fairSIM result
% is not as good as the HR images provided by commerical microscope (e.g. GE OMX SR). In this situation, we can skip the fairSIM part 
% and directly use the HR images for pSIM reconstruction.

%% Parameter k
% Variable k should be a matrix of m*2, where m denotes the number of illumination patterns; 
k(1,:) = [142.656,-152.367];
k(2,:) = [-58.567,-200.1];
k(3,:) = [-201.833,-51.678];

%% Parameter phase
% Variable phase should be a vector of 1*m, where m (3 for 2D SIM) denotes the number of illumination patterns;
% phase(i) denotes the referred phase in the ith pattern, the three
% phases in this pattern should be phase(i), phase(i)+2/3*pi and phase(i)+4/3*pi
phase = [1.523,2.404,-2.033];

%% Directory of the raw images
% The number of the raw images should be m*n, where m (3 for 2D SIM) denotes the number of
% illumination patterns and n (3 for 2D SIM) denotes the number of the phases in each
% pattern.
% The input images in the folder should named as '1.tif', '2.tif', ... 'm*n.tif'  
% The test image is the alexa 488 label actin in BAPE cell
read_dir = 'Input\';

%% Directory of the output images
saveDir = 'Output\';

%% Calibration
% Calibration is to compensate uneven illumination intensity among
% different directions. Different setup has different performance.
% calibrationFlag: whether or not use system calibration
% if use calibration, give the directory of the calibamp
calibrationFlag = false; 
if calibrationFlag
    calib1path = 'Calib\calib1.tif'; 
    calib2path = 'Calib\calib2.tif'; 
end

%% whether display the result
displayFlag = true;


%% whether use otf attentuation or not
attFlag = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The PSIM begins
addpath('functions\');
close all;

%% read images
for i = 1: 1: 9
    rawImage(:,:,i) = double(imread([read_dir num2str(i) '.tif'])); 
end

w = size(rawImage,2);
h = size(rawImage,1);

imgs = zeros(size(rawImage));
for i = 1: 1: 9
    imgs(:,:,i) = fadeBorderCos(rawImage(:,:,i),10); 
end

inFFT = zeros(size(rawImage));
for i = 1: 1: 9
    inFFT(:,:,i) = fft2(imgs(:,:,i)); 
end


%% global parameters
showDialog =  false;

% the band in SIM setup
nrBands  = 2;   % beam number
nrDirs   = 3;   % illumination direction number
nrPhases = 3;   % phase number


emWavelen = 528;    % emission wavelength
otfNA     = 1.4;    % NA
otfCorr   = 0.31;
pxSize    = 0.080;  % pixel size

wienParam   = 0.05; % Wiener filter parameter

attStrength = 0.995;
attFWHM     = 1.2;
doAttenuation = true; % otf attentuation

apoB=0.9;   % apo parameter
apoF=2;

%% OTF
otfPr = otfGenerator( otfNA, emWavelen, otfCorr, attStrength, attFWHM);
otfPr.doAttenuation = doAttenuation;

%% SIM
param = SimParamCreate(nrBands, nrDirs, nrPhases, w, pxSize);

%%
otfPr.vecCyclesPerMicron = param.cyclesPerMicron;

%% OTFs
OTF_raw = oftmatrix( otfPr, w, h);
% OTF_general = oftmatrix( otfPr, 2*w, 2*h);

%% calibration 
%% if the calibrationFlag is set as "false", this module will not execute
if calibrationFlag
    % calibration map is saved as ratio*10000 in pixel
    calib1 = double(imread(calib1path))/10000;
    calib2 = double(imread(calib2path))/10000;
    
    img_calib = zeros(size(inFFT));
    mask = OTF_raw>0;
    for i = 1: 1: 9
        tmpinFFT = inFFT(:,:,i);
        tmpoutFFT = zeros(size(tmpinFFT));
%         tmpoutFFT(mask) = tmpinFFT(mask);
        tmpoutFFT(mask) = tmpinFFT(mask) ./  OTF_raw(mask);
        tmpout = ifft2(tmpoutFFT);
        img_calib(:,:,i) = tmpout;
    end
    img_calib(:,:,4) = img_calib(:,:,4)./ calib1;
    img_calib(:,:,5) = img_calib(:,:,5)./ calib1;
    img_calib(:,:,6) = img_calib(:,:,6)./ calib1;
    
    img_calib(:,:,7) = img_calib(:,:,7)./ calib2;
    img_calib(:,:,8) = img_calib(:,:,8)./ calib2;
    img_calib(:,:,9) = img_calib(:,:,9)./ calib2;
    
    for i = 1: 1: 9
         tmpF = fft2(img_calib(:,:,i));
         tmpF(mask) = tmpF(mask) .* OTF_raw(mask);
         tmpF(~mask) = 0;
         inFFT(:,:,i) = tmpF;
    end
end

%% Wienar Filter
wFilter = zeros(2*h,2*w);
xx = 1: 1: 2*w;
yy = 1: 1: 2*h;
[x,y] = meshgrid(xx,yy);
x(:,1:w) = x(:,1:w)-1;
x(:,w+1:2*w) = x(:,w+1:2*w)-2*w - 1;

y(1:h,:) = -(y(1:h,:)-1);
y(h+1:2*h,:) = 2*h - (y(h+1:2*h,:)-1);

for d = 1: 1: 3
    for b = 1: 1: 2
        rad1 = sqrt( (x-k(d,1)*(b-1)).^2 + (y-k(d,2)*(b-1)).^2 ) * otfPr.vecCyclesPerMicron;
        rad2 = sqrt( (x+k(d,1)*(b-1)).^2 + (y+k(d,2)*(b-1)).^2 ) * otfPr.vecCyclesPerMicron;
        v_out1 = getOtfVal_m( otfPr, b, rad1, attFlag);
        v_out2 = getOtfVal_m( otfPr, b, rad2, attFlag);
        wFilter = wFilter + v_out1 + v_out2;
    end
end


%% SIM Reconstruction (matlab version of fairSIM)
fullResult = zeros(2*h,2*w);
p_zero_component = zeros(h,w,3);

for angIdx = 1: 1: 3
    kx = k(angIdx,1);
    ky = k(angIdx,2);

    M = [1, 0.5*exp( 1i * (phase(angIdx))), 0.5*exp( -1i * phase(angIdx));
        1, 0.5*exp( 1i * (phase(angIdx)+pi*2/3)), 0.5*exp( -1i * (phase(angIdx)+pi*2/3));
        1, 0.5*exp( 1i * (phase(angIdx)+pi*4/3)), 0.5*exp( -1i * (phase(angIdx)+pi*4/3))];
    invM = inv(M);
    
    separate = zeros(size(inFFT,1),size(inFFT,2),3);
    separate(:,:,1) = invM(1,1) * inFFT(:,:,(angIdx-1)*3+1)+invM(1,2) * inFFT(:,:,(angIdx-1)*3+2)+invM(1,3) * inFFT(:,:,(angIdx-1)*3+3);
    separate(:,:,2) = invM(2,1) * inFFT(:,:,(angIdx-1)*3+1)+invM(2,2) * inFFT(:,:,(angIdx-1)*3+2)+invM(2,3) * inFFT(:,:,(angIdx-1)*3+3);
    separate(:,:,3) = invM(3,1) * inFFT(:,:,(angIdx-1)*3+1)+invM(3,2) * inFFT(:,:,(angIdx-1)*3+2)+invM(3,3) * inFFT(:,:,(angIdx-1)*3+3);
    
    separate_OTF = zeros(size(inFFT,1),size(inFFT,2),3);
    for i=1: 1: 3  
        separate_OTF(:,:,i) = applyOtf( otfPr, separate(:,:,i), 0, 0, attFlag);
    end
    
    p_zero_component(:,:,angIdx) = separate_OTF(:,:,1);
    
    shifted = zeros(2*w, 2*h,5);
    shifted(:,:,1) = pasteFreq( separate_OTF(:,:,1));
    
    pos = 3;
    neg = 2;
    
    shifted(:,:,pos) = pasteAndFourierShift( separate_OTF(:,:,pos), kx, ky );
    shifted(:,:,neg) = pasteAndFourierShift( separate_OTF(:,:,neg), -kx, -ky );
	
    shifted_mask = zeros(size(shifted));
    shifted_mask(:,:,1) = maskOtf( otfPr, shifted(:,:,1),  0,  0);
    shifted_mask(:,:,pos) = maskOtf( otfPr, shifted(:,:,pos),  kx,  ky);
    shifted_mask(:,:,neg) = maskOtf( otfPr, shifted(:,:,neg),  -kx,  -ky);
    
    thisD = shifted_mask(:,:,1)+shifted_mask(:,:,pos)+shifted_mask(:,:,neg);
    
    
    for i = 1: 1: 3
        fullResult = fullResult + shifted_mask(:,:,i);
    end
end

denom = 1./(wFilter+wienParam^2);
fullResult_filtered = fullResult .* denom;

apo = writeApoVector( otfPr, apoB, apoF, 2*h, 2*w);
fullResult_filtered = fullResult_filtered .* apo;


%% pSIM reconstruction
theta = atan2(k(:,2),k(:,1));
M_ld = 0.5*[1 0.5*exp(2i*theta(1)) 0.5*exp(-2i*theta(1));
    1 0.5*exp(2i*theta(2)) 0.5*exp(-2i*theta(2));
    1 0.5*exp(2i*theta(3)) 0.5*exp(-2i*theta(3))];
invM_ld = inv(M_ld);

p_component(:,:,1) = invM_ld(1,1) * p_zero_component(:,:,1) + invM_ld(1,2) * p_zero_component(:,:,2) + invM_ld(1,3) * p_zero_component(:,:,3);
p_component(:,:,2) = invM_ld(2,1) * p_zero_component(:,:,1) + invM_ld(2,2) * p_zero_component(:,:,2) + invM_ld(2,3) * p_zero_component(:,:,3);
p_component(:,:,4) = invM_ld(3,1) * p_zero_component(:,:,1) + invM_ld(3,2) * p_zero_component(:,:,2) + invM_ld(3,3) * p_zero_component(:,:,3);

psim_f(:,:,1) = zeros(size(fullResult)); 
psim_f(:,:,2) = fftshift(pasteFreq(p_component(:,:,4))); psim_f(:,:,4) = fftshift(pasteFreq(p_component(:,:,2)));


psim_f(:,:,3) = fftshift(fullResult);
for i = 1: 1: 4
    psim_f(:,:,i) = psim_f(:,:,i).*fftshift(denom);
    psim_f(:,:,i) = psim_f(:,:,i).*fftshift(apo);
end

% % % We can also directly use the SR images provided by the software accompanied 
% % % by the commerical microscope and skip the multiplying with denom and apo as:
% % psim_f(:,:,3) = fftshift(fft2(SIMSRImg))*scaleFactor;
% % for i = 1: 1: 4
% %     if i ~= 3
% %         psim_f(:,:,i) = psim_f(:,:,i).*fftshift(denom);
% %         psim_f(:,:,i) = psim_f(:,:,i).*fftshift(apo);
% %     end
% % end
% % % "SIMSRImg is the HR reconstructed from the software of your commerial
% % % SIM scope
% % % The scaleFactor is applied since the reconstrucved SR image and the raw
% % % LR images may not in the same scale.


% inverse Fourier Transform
psim = abs(ifft(ifft(ifft(ifftshift(psim_f),[],1),[],2),[],3));

% display the psim image as dipole orientation indicited by colorwheel
[sim, psim_om,cm, ~] = PSIM_display(psim,min(psim(:)),max(psim(:)),false);

%% Display the result 
if displayFlag
    figure;
    imshow(sum(rawImage,3),[]);
    title('WF');
    
    figure;
    imshow(sim,[]);
    title('SIM');
    
    figure;
    imshow(psim_om);
    title('PSIM');
end

%% output result
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

imwrite(uint16(sim), [saveDir 'sim','.tif']);

if calibrationFlag
    imwrite(uint8(psim_om*255), [saveDir 'pSIM_cal', '.png']);
else
    imwrite(uint8(psim_om*255), [saveDir 'pSIM_noncal', '.png']);
end

% Color wheel
imwrite(uint8(cm*255), [saveDir 'cm.png']);

