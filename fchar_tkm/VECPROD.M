% Векторное произведение трехм. векторов vecprod(a,b)
function y=vecprod(a,b);
   y(1)=a(2)*b(3)-a(3)*b(2);
   y(2)=-a(1)*b(3)+a(3)*b(1);
   y(3)=a(1)*b(2)-a(2)*b(1);
   [m,n]=size(a);
   if m>1 y=y(:); end
