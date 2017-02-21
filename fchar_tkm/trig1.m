function G=trig1(a,b,c)
%TRIG1 ���������� ���� 
%������� ������������-
%������� ��������� 
%acosgamma ����
%bsingamma ����� c
%��������� trig1
zn=sqrt(a^2+b^2);
kc=a/zn;
ks=b/zn;
kn=c/zn;
if ks<0
    kc=-kc;
    ks=-ks;
    kn=-kn;
end
if abs(kc)>1    %����������� �������������� ������
    kc=sign(kc);
end
if abs(kn)>1
    kn=sign(kn);
end
gmn=acos(kc);
G(1)=gmn-acos(kn);
G(2)=gmn+acos(kn);
if G(2)>pi
   G(2)=-2*pi+G(2);
end



