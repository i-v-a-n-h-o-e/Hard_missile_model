function G=trig2m(a1,b1,c1,a2,b2,c2)
%TRIG2 Вычисление решения системы двух тригон. уравнений 
%acosgamma + bsingamma = c
%Вычислить G=trig2(a1,b1,c1,a2,b2,c2)
det=a1*b2-a2*b1;
C=(c1*b2-c2*b1)/det;
S=(c2*a1-c1*a2)/det;
if abs(C)>1
    C=sign(C);    %Парирование вычислительных ошибок
end
if S>=0
   G=acos(C);
else
   G=-acos(C);
end






    





