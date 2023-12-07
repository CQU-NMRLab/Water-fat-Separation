
function [W,F]=dixonprior(in_m,out_m,in_pp,out_pp,TE)

% Author:  Cai Wan
%          Modified from hzr0071
% Created: December 2023
% Last Modified: December 5th, 2023

% Inputs:
%   in_m,out_m,in_pp,out_pp:  The amplitude and phase of magnetic resonance imaging
%   TE:                       Timing of the offset

% Outputs:
%   W, F                   :  Water-only and fat-only images

% Reference:  Two-Point Water-Fat Imaging With Partially-Opposed-Phase (POP) Acquisition: An Asymmetric Dixon Method

algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 4.70;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29];
algoParams.species(2).relAmps = [88 642 58 62 58 6 39 10 37];
gyro = 42.577481;
FieldStrength=0.05;
flextime=0.015;
for m = 1:length(algoParams.species)
algoParams.species(m).relAmps=algoParams.species(m).relAmps/sum(algoParams.species(m).relAmps);
end
omega_p = 2*pi*2.323*(algoParams.species(1).frequency-algoParams.species(2).frequency);
alpha=exp(complex(0,(TE+flextime)'*omega_p))*algoParams.species(2).relAmps';

that=angle((alpha));
if that>=0
    that=that;
else
    that=pi+(pi-abs(that));
end
    
mb=40;
in=(in_m);
out=(out_m);

in_p=(in_pp);
out_p=(out_pp);
[sizem,sizen]=size(out);

J1=abs(in);
J2=out.*exp((out_p-in_p)*1i); % Eq.8 from reference
 
M1=abs(in);
M2=abs(out);
ang=that;
t1=clock;

[S,B]=getSB(M1,M2,ang); % Eq.1 and 2
Pu=J2./(B+S.*exp(ang*1i));
Pv=J2./(S+B.*exp(ang*1i)); % Eq.9 and 10

 Pu=Pu./abs(Pu);
 Pv=Pv./abs(Pv);

 for m=1:sizem
   for n=1:sizen
       if isnan(Pu(m,n))
           Pu(m,n)=0;
       end
       if isnan(Pv(m,n))
           Pv(m,n)=0;
       end
       
   end
end
P0=get_init_phase(Pu,Pv);

Pn=smoth(P0,mb); % Eq.12

for m=1:200
    Pn=get_new_phase(Pn,Pu,Pv); % Eq. 13
    Pn=smoth(Pn,mb);
end


Ps=Pn./abs(Pn);
Jc=J2.*conj(Ps);

A=[1,1;1,cos(ang);0,sin(ang)];
A=((A'*A)\A');
[m,n]=size(in);
J=zeros(3,1,m,n);
Xls=zeros(2,1,m,n);
W=zeros(m,n);
F=zeros(m,n);
for a=1:m
    for b=1:n
   J(:,:,a,b)=[J1(a,b),real(Jc(a,b)),imag(Jc(a,b))];
   Xls(:,:,a,b)=A*J(:,:,a,b); % Eq. 20
   W(a,b)=Xls(1,1,a,b);
   F(a,b)=Xls(2,1,a,b);
    end
end
t2=clock;
etime(t2,t1)
end
