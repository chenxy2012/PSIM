# PSIM
Matlab source code for the reconstruction of Polarized Structured Illumination Microscopy (PSIM)

This program is developed based on fairSIM.

In the reconstruction process, the polarization of the sample can be resolved, and displayed in pseudo color mode.

PSIM reconstruction required precise SIM parameters, which is available by the accompanyed parameter file of the 
commerical SIM system or estimated from fairSIM. 


A testfile is available in the folder "input".

This program starts with PSIM.m at the root directory. The UI interface and more functions are under development.

References: 

The principle of pSIM has been published at https://arxiv.org/abs/1712.05092 .

FairSIM plugin (for ImageJ) is available at https://github.com/fairSIM

The publication accompanying the fairSIM plugin has been published in Nature Communications (open access):
Marcel Müller, Viola Mönkemöller, Simon Hennig, Wolfgang Hübner, Thomas Huser (2016).
"Open-source image reconstruction of super-resolution structured illumination microscopy data in ImageJ",
Nature Communications, doi: 10.1038/ncomms10980

For any questions, please contact: chenxy16@mails.tsinghua.edu.cn