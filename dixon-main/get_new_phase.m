function NPn=get_new_phase(Pn,Pu,Pv)
[m,n]=size(Pn);
NPn=zeros([m,n]);
P1=abs((Pu)-(Pn));
P2=abs((Pv)-(Pn));
for a=1:m
    for b=1:n
    if P1(a,b)<P2(a,b)
        NPn(a,b)=Pu(a,b);
    elseif P1(a,b)>P2(a,b)
        NPn(a,b)=Pv(a,b);
    elseif P1(a,b)==P2(a,b)
        NPn(a,b)=0;
    end
    end
end

end