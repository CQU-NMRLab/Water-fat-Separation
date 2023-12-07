function [PS,PC]=getweight(Pu,Pv,num)
[m,n]=size(Pu);

PS=zeros([m,n]);
PC=zeros([m,n]);
for i=1:num
    for j=1:num
        if abs(Pu(i,j)) >= abs(Pv(i,j))
            PS(i,j)=(Pu(i,j));
            PC(i,j)=(Pv(i,j));
        elseif abs(Pu(i,j)) < abs(Pv(i,j))
            PS(i,j)=(Pv(i,j));
            PC(i,j)=(Pu(i,j));
        else 
            PS(i,j)=0;
            PC(i,j)=0;
        end  
    end
end

end