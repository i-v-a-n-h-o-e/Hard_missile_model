function [Y]=dispom(Mp,kf);
% �८�ࠧ�� �����. ��������� � ��ப�� ������ Mp ��� ���������
% �ࠪ���� ���� � kf ࠧ.

[n,m]=size(Mp);   
for s=1:m-1
        Mp(:,s)=Mp(:,s)/(kf^(m-s));
end
Y=Mp;
