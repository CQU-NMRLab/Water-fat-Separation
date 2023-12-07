function P0=get_init_phase(Pu,Pv)
PS=zeros(size(Pu));
PC=zeros(size(Pu));
[ms,ns]=size(Pu);
% for m=0:ms/8-1
%     for n=0:ns/8-1
%         bm=m*8;
%         bn=n*8;
%         [PS(bm+1:bm+8,bn+1:bn+8),PC(bm+1:bm+8,bn+1:bn+8)] =getweight(Pu(bm+1:bm+8,bn+1:bn+8),Pv(bm+1:bm+8,bn+1:bn+8),8)  ;     
%     end
% end
% MS=(sum(sum(abs(PS))));
% MC=(sum(sum(abs(PC))));
% C=abs(MS.^2-MC.^2)/abs(MS.^2+MC.^2);
% P0=PS.*C;

 P0=(Pv+Pu)/2;

 
end