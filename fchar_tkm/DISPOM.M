function [Y]=dispom(Mp,kf);
% Преобразует коэфф. полиномов в строках матрицы Mp для изменения
% характерных частот в kf раз.

[n,m]=size(Mp);   
for s=1:m-1
        Mp(:,s)=Mp(:,s)/(kf^(m-s));
end
Y=Mp;
