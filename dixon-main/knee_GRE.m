% Project: Water fat separation at 50mT MRI
% Filename: knee_GRE.m
%   Template for knee fat/water separation using
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
addpath([BASEPATH 'knee data/']);

%% Add knee image k-space data
[coil1,seqParam1] = getKSpace1('Scanwxcge3dte241.mrd');
[coil2,seqParam2] = getKSpace1('Scanwxcge3dte242.mrd');

coe=quandratureCoilCoe(coil1,coil2); % quandratureCoil
newspace=coil1+coe*coil2;
newspace1=newspace;

img1=ifftshift(ifftn(coil1));
img2=ifftshift(ifftn(coil2));
imgCo=ifftshift(ifftn(newspace));
imgCo1=ifftshift(ifftn((newspace1)));

img1=rot90(img1,-3);
img2=rot90(img2,-3);
imgCo1=fliplr(rot90(imgCo1,-3));

[coil3,seqParam3] = getKSpace1('Scanwxcge3dte441.mrd');
[coil4,seqParam4] = getKSpace1('Scanwxcge3dte442.mrd');

coe=quandratureCoilCoe(coil3,coil4);
newspace3=coil3+coe*coil4;
newspace4=newspace3;

img3=ifftshift(ifftn(coil3));
img4=ifftshift(ifftn(coil4));
imgCo3=ifftshift(ifftn(newspace3));
imgCo4=ifftshift(ifftn((newspace4)));

img3=rot90(img3,-3);
img4=rot90(img4,-3);
imgCo4=fliplr(rot90(imgCo4,-3));

[coil5,seqParam5] = getKSpace1('Scanwxcge3dte641.mrd');
[coil6,seqParam6] = getKSpace1('Scanwxcge3dte642.mrd');

coe=quandratureCoilCoe(coil5,coil6);
newspace5=coil5+coe*coil6;
newspace6=newspace5;

img5=ifftshift(ifftn(coil5));
img6=ifftshift(ifftn(coil6));
imgCo5=ifftshift(ifftn(newspace5));
imgCo6=ifftshift(ifftn((newspace6)));

img5=rot90(img5,-3);
img6=rot90(img6,-3);
imgCo6=fliplr(rot90(imgCo6,-3));

%% Add priori phantom image k-space data
[coil7,seqParam7] = getKSpace1('Scanwxcge3dte24sagphan1.mrd');
[coil8,seqParam8] = getKSpace1('Scanwxcge3dte24sagphan2.mrd');

coe=quandratureCoilCoe(coil7,coil8);
newspace7=coil7+coe*coil8;
newspace8=newspace7;

img7=ifftshift(ifftn(coil7));
img8=ifftshift(ifftn(coil8));
imgCo7=ifftshift(ifftn(newspace7));
imgCo8=ifftshift(ifftn((newspace8)));

img7=rot90(img7,-3);
img8=rot90(img8,-3);
imgCo8=fliplr(rot90(imgCo8,-3));

[coil9,seqParam9] = getKSpace1('Scanwxcge3dte44sagphan1.mrd');
[coil10,seqParam10] = getKSpace1('Scanwxcge3dte44sagphan2.mrd');

coe=quandratureCoilCoe(coil9,coil10);
newspace9=coil9+coe*coil10;
newspace10=newspace9;

img9=ifftshift(ifftn(coil9));
img10=ifftshift(ifftn(coil10));
imgCo9=ifftshift(ifftn(newspace9));
imgCo10=ifftshift(ifftn((newspace10)));

img9=rot90(img9,-3);
img10=rot90(img10,-3);
imgCo10=fliplr(rot90(imgCo10,-3));

[coil11,seqParam11] = getKSpace1('Scanwxcge3dte64sagphan1.mrd');
[coil12,seqParam12] = getKSpace1('Scanwxcge3dte64sagphan2.mrd');

coe=quandratureCoilCoe(coil11,coil12);
newspace11=coil11+coe*coil12;
newspace12=newspace11;

img11=ifftshift(ifftn(coil11));
img12=ifftshift(ifftn(coil12));
imgCo11=ifftshift(ifftn(newspace11));
imgCo12=ifftshift(ifftn((newspace12)));

img11=rot90(img11,-3);
img12=rot90(img12,-3);
imgCo12=fliplr(rot90(imgCo12,-3));

S14=(imgCo1(:,:,14));
S24=(imgCo4(:,:,14));
S34=(imgCo6(:,:,14));

Sp14=(imgCo8(:,:,14));
Sp24=(imgCo10(:,:,14));
Sp34=(imgCo12(:,:,14));

DATA4=cat(3,S14,S24,S34);
TE1=str2num(seqParam1.te)*1e-3;
TE2=str2num(seqParam3.te)*1e-3;
TE3=str2num(seqParam5.te)*1e-3;
TE=[TE1,TE2,TE3];
[r2s4,t2s, s0] = R2star_ARLO_mag(DATA4,TE);


SS14=(abs(S14).*exp((angle(Sp14)-angle(S14))*1i)).*exp((TE1-TE1)*r2s4);

SS24=(abs(S24).*exp((angle(Sp24)-angle(S24))*1i)).*exp((TE2-TE1)*r2s4);

[W,F]=dixonprior((SS14),(SS24),angle(SS14),angle(SS24),TE2);  % two-point dixon 

FF= abs(F)./(abs(W)+abs(F));

ha = tight_subplot(1,5,[.001 .01],[.001 .001],[.001 .025]);
axes(ha(1)); 
imshow(abs(S14(55:432,:)),[]);
axes(ha(2));
imshow(abs(S24(55:432,:)),[]);
axes(ha(3)); 
imshow((abs(W(55:432,:))),[]);
axes(ha(4)); 
imshow((abs(F(55:432,:))),[]);
axes(ha(5)); 
imshow((abs(FF(55:432,:))),[]);
