% Project: Water fat separation at 50mT MRI
% Filename: mainpriorGRE_T2x.m
%   Template for phantom fat/water separation using
%   the two-point dixon algorithm incorporating the R2* effect and priori information.
%
% 
% Modified:
%
% Last Modified: 2023/12/05
% Created:       2023/12/05
% Written by Cai Wan and Zheng Xu, 2023
% (c) School of Electrical Engineering, Chongqing University

% Please send comments and suggestions to
% 20181113069t@cqu.edu.cn or xuzheng@cqu.edu.cn 
clc; clear all;
%% Add to matlab path
BASEPATH = '../';
addpath([BASEPATH 'phantom data/']);

%% Add phantom image k-space data
[coil1,seqParam1] = getKSpace1('ScanGRE2-20ms1.mrd');
[coil2,seqParam2] = getKSpace1('ScanGRE2-20ms2.mrd');

coe=quandratureCoilCoe(coil1,coil2);
newspace=coil1+coe*coil2;
newspace1=newspace;

img1=ifftshift(ifftn(coil1));
img2=ifftshift(ifftn(coil2));
imgCo=ifftshift(ifftn(newspace));
imgCo1=ifftshift(ifftn((newspace1)));

img1=rot90(img1,-2);
img2=rot90(img2,-2);
imgCo1=(rot90(imgCo1,-2));

[coil3,seqParam3] = getKSpace1('ScanGRE2-40ms1.mrd');
[coil4,seqParam4] = getKSpace1('ScanGRE2-40ms2.mrd');

coe=quandratureCoilCoe(coil3,coil4);
newspace3=coil3+coe*coil4;
newspace4=newspace3;

img3=ifftshift(ifftn(coil3));
img4=ifftshift(ifftn(coil4));
imgCo3=ifftshift(ifftn(newspace3));
imgCo4=ifftshift(ifftn((newspace4)));

img3=rot90(img3,-2);
img4=rot90(img4,-2);
imgCo4=(rot90(imgCo4,-2));

[coil5,seqParam5] = getKSpace1('ScanGRE2-60ms1.mrd');
[coil6,seqParam6] = getKSpace1('ScanGRE2-60ms2.mrd');

coe=quandratureCoilCoe(coil5,coil6);
newspace5=coil5+coe*coil6;
newspace6=newspace5;

img5=ifftshift(ifftn(coil5));
img6=ifftshift(ifftn(coil6));
imgCo5=ifftshift(ifftn(newspace5));
imgCo6=ifftshift(ifftn((newspace6)));

img5=rot90(img5,-2);
img6=rot90(img6,-2);
imgCo6=(rot90(imgCo6,-2));

%% Add priori phantom image k-space data
[coil7,seqParam7] = getKSpace1('ScanGRE-corrent-20ms1.mrd');
[coil8,seqParam8] = getKSpace1('ScanGRE-corrent-20ms2.mrd');

coe=quandratureCoilCoe(coil7,coil8);
newspace7=coil7+coe*coil8;
newspace8=newspace7;

img7=ifftshift(ifftn(coil7));
img8=ifftshift(ifftn(coil8));
imgCo7=ifftshift(ifftn(newspace7));
imgCo8=ifftshift(ifftn((newspace8)));


img7=rot90(img7,-2);
img8=rot90(img8,-2);
imgCo8=(rot90(imgCo8,-2));

[coil9,seqParam9] = getKSpace1('ScanGRE-corrent-40ms1.mrd');
[coil10,seqParam10] = getKSpace1('ScanGRE-corrent-40ms2.mrd');

coe=quandratureCoilCoe(coil9,coil10);
newspace9=coil9+coe*coil10;
newspace10=newspace9;

img9=ifftshift(ifftn(coil9));
img10=ifftshift(ifftn(coil10));
imgCo9=ifftshift(ifftn(newspace9));
imgCo10=ifftshift(ifftn((newspace10)));


img9=rot90(img9,-2);
img10=rot90(img10,-2);
imgCo10=(rot90(imgCo10,-2));

[coil11,seqParam11] = getKSpace1('ScanGRE-corrent-60ms1.mrd');
[coil12,seqParam12] = getKSpace1('ScanGRE-corrent-60ms2.mrd');

coe=quandratureCoilCoe(coil11,coil12);
newspace11=coil11+coe*coil12;
newspace12=newspace11;

img11=ifftshift(ifftn(coil11));
img12=ifftshift(ifftn(coil12));
imgCo11=ifftshift(ifftn(newspace11));
imgCo12=ifftshift(ifftn((newspace12)));

img11=rot90(img11,-2);
img12=rot90(img12,-2);
imgCo12=(rot90(imgCo12,-2));

S1=(imgCo1(:,:,4));
S2=(imgCo4(:,:,4));
S3=(imgCo6(:,:,4));

Sp1=(imgCo8(:,:,4));
Sp2=(imgCo10(:,:,4));
Sp3=(imgCo12(:,:,4));

DATA=cat(3,S1,S2,S3);
TE1=str2num(seqParam1.te)*1e-3;
TE2=str2num(seqParam3.te)*1e-3;
TE3=str2num(seqParam5.te)*1e-3;
TE=[TE1,TE2,TE3];
[r2s,t2s, s0] = R2star_ARLO_mag(DATA,TE); % ARLO calculates R2*map


SS1=(abs(S1).*exp((angle(Sp1)-angle(S1))*1i)).*exp((TE1-TE1)*r2s);
SS3=(abs(S3).*exp((angle(Sp3)-angle(S3))*1i)).*exp((TE3-TE1)*r2s); % Eq. 2

[W,F]=dixonprior(abs(SS1),abs(SS3),angle(SS1),angle(SS3), TE3-TE1);  % two-point dixon 

FF= abs(F)./(abs(W)+abs(F)); % Calculating Fat Fraction

ha = tight_subplot(1,5,[.001 .001],[.001 .001],[.001 .001]); 
axes(ha(1)); 
imshow(abs(SS1),[]);
axes(ha(2)); 
imshow(abs(SS3),[]);
axes(ha(3)); 
imshow(abs(W),[]);
axes(ha(4)); 
imshow(abs(F),[]);
axes(ha(5)); 
imshow(abs(FF),[]);
