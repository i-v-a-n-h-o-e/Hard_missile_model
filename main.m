%%*********************************************
%
% ���������� ���� ��� ����������� ������� ������� � 
% ����������� �������� ����������� � ���� ������ ������� 
%
% ��������� ImportFromMatAndExcel, matrixABGinit, contour
% ��������� boderus
%
%***********************************************
tic
%clear;
%clc;

%% ������ � ��
prj.name = '����-2.1�';
prj.pn = 2175;
prj.ndus = 1;
prj.ntonsall = 3; % ����� ��������� ����� ������� ��������� 
%% ����������� ������������ ������
prj.ntons = 3; %{0,1,2,3} ����� ����������� ����� ������� ��������� ������� ������ 
prj.ntanks = 2; % {0,2} ����� ����� � ������ �����������
%����� �����
prj.spaceport = '��������';%{'��������','�������'}
prj.season = '����';%{'����','����'}
prj.windtype = '���������';%{'���������','��������'}

%prj.mode = '�����������';
%prj.mode = '����������';
prj.mode = '����������';

%prj.path = 'C:\Users\admin\Documents\MATLAB\';%���� �� ������
prj.path = 'D:\YandexDisk\Documents\Ivan\!!Tasks&Projects\RocketStabilizationSystem\';%���� �� �����  
prj.id = [prj.path,prj.name...
    ' �� = ' num2str(prj.pn)];

dm.T0 = 0.03; %������ �������������,� 

t = 80; %��������������� ������ �������
as.a0 = 19; % ����� ����������������
as.a1 = 0.8; % ���������������� �� �������
as.n0 = 0; % ������������ �� ��
as.n1 = 0.02; %�������������q �� ��

run([prj.path,'ImportFromMatAndExcel.m']);
run([prj.path,'matrixABGinit.m']);
run([prj.path,'contour.m']);

boderus_A_v2(W,3);
title({[prj.name, ' �����: ',prj.mode, ' t = ',num2str(t),' c'],'����� �������'});
figure;
step(feedback(W,1,1));
grid on;
toc
