function objd=corr_c2d(objc,T);
%Ïğåîáğàçîâàíèå àíàëîãîâîãî êîğğåêòîğà (ñèñòåìû) â äèñêğåòíóş ôîğìó.
%objd=corr_c2d(objc,T);
objc=tf(objc)
[numd,dend]=refwz(T,objc.num{1},objc.den{1});
objd=tf(numd,dend,T);

function [A1,B1]=refwz(T,A,B);
% Ïğåîáğàçîâàíèå ÏÔ îò ïñåâäî÷àñòîû w, çàäàííîé ïîëèíîìàìè A (÷èñëèòåëü), 
% B(çíàìåíàòåëü), â z-ÏÔ íà îñíîâå ñîîòíîøåíèÿ w=(2/Ò)*(z-1)/(z+1).

n=[length(A),length(B)];  nm=max(n);
S={A,B};
for c=1:2
   M=S{c};
	rab=2/T;
	for b=length(M)-1:-1:1
   	M(b)=M(b)*rab;
   	rab=rab*2/T;
   end
   S{c}=M;
end
Z=zeros(2,nm);
for s=1:2
        V=S{s};
        for s1=1:n(s)
                Pzn=1;
                Pzp=1;
                for k=1:n(s)-s1
                        Pzn=conv(Pzn,[1,-1]);
                end
                for k=1:nm-n(s)+s1-1
                        Pzp=conv(Pzp,[1,1]);
                end
                Z(s,:)=Z(s,:)+V(s1)*conv(Pzn,Pzp);
        end
end
if abs(Z(2,1))>1e-6
        Z=Z/Z(2,1);
end
B1=Z(2,:);
A1=Z(1,:);  
