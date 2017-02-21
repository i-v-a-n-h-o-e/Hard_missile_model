%%*********************************************
%
% ���������� ������� ������� ������� ��������� (���) � �������� �������
% �� ����� .mat
% ������ ������� ������� � ������ ������� t. ������� ������� �� ������� ����������
% (������� A,B,C,G  dm. ), �������� ������������(as), �������� ������� 
% ������ �������������� ������ ������������� 
% ��� ��������� � ����� ���� �� ����� ������� dm.t
%
% ���������� ������ ������� t
% ������� �, B, C, G � ���� *.mat
% ��������������� ������ matrixABGinit.m
%
%***********************************************

load([prj.path, 'DUS&rp.mat']); %�������� ������� ��� � ��
sl = ceil(t/dm.T0) + 1; %����� ���� ������ A,B,C,G
Wobj = ss(dm.A(:,:,sl),dm.B(:,:,sl),dm.C(:,:,sl),0);
l = as.l(sl);
    
Cnt = append(Wrp,Wobj,Wdus*eye(prj.ndus));
Q = [2 1;
     (3:2+prj.ndus)' (3:2+prj.ndus)'];
inputs = 1;
if prj.ndus == 1
    a = 5;
else
   a = (6:7)';
end;   
outputs = [2; a; 3+prj.ndus];
W1 = connect(Cnt,Q,inputs, outputs);
clearvars Cnt Q inputs outputs
set(W1,'OutputDelay',[0.0214; 0.0114*ones(prj.ndus,1); 0.0214]);
%W1 = c2d(W1,dm.T0);
i = ss(tf(1,[1 0]));
if prj.ndus == 2
    alpha = oscil.dfx(2,1,t)/(oscil.dfx(2,1,t)-oscil.dfx(1,1,t));
    complex = ss([1-alpha alpha]);
    Cnt = append(W1,complex,as.a0,as.a1,as.n0*i,as.n1,l);
    Q = [3 3 0 0;
        2 2 0 0 ;
        5 5 0 0;
        7 4 -10 8;
        6 4 -10 0;
        8 5 0 0;
        4 1 7 9];
        outputs = 6;
else    
    Cnt = append(W1,as.a0,as.a1,as.n0*i,as.n1,l);
    Q = [3 2 0 0;
        6 2 0 0;
        5 3 -8 6;
        4 3 -8 0;
        2 1 5 7];
        outputs = 4;
end;
inputs = 1;
W = connect(Cnt,Q,inputs, outputs);
clearvars Cnt Q inputs outputs sl