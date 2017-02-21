function y=lint(x,tabx,taby);
% Линейная интерполяция. x, tabx - векторы, taby - матрица.
[m,l]=size(taby);
for k=1:length(x)
        for i=2:l
                if (tabx(i-1)-x(k))*(tabx(i)-x(k))<=0
                        I=i;
                        break
                end
        end
        for i=1:m
           y(i,k)=(taby(i,I-1)*(tabx(I)-x(k))+taby(i,I)*(x(k)-tabx(I-1)))...
           /(tabx(I)-tabx(I-1));
        end
end

