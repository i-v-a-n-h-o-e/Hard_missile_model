function [A1,B1]=refzw(A,B);
% Преобразование z-ПФ, заданной полиномами A (числитель), 
% B(знаменатель), в w-ПФ на основе соотношения z=(1+w)/(1-w).

n=[length(A),length(B)];  nm=max(n);
Z=zeros(2,nm);
for s=1:2
        if s==1
                V=A;
        else
                V=B;
        end
        for s1=1:n(s)
                Pzn=1;
                Pzp=1;
                for k=1:n(s)-s1
                        Pzn=conv(Pzn,[1,1]);
                end
                for k=1:nm-n(s)+s1-1
                        Pzp=conv(Pzp,[-1,1]);
                end
                Z(s,:)=Z(s,:)+V(s1)*conv(Pzn,Pzp);
        end
end
if abs(Z(2,1))>1e-6
        Z=Z/Z(2,1);
end
B1=Z(2,:);
A1=Z(1,:);
