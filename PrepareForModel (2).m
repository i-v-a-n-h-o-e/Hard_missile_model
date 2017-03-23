%%*********************************************
%
% ������� ������ ��� ������ � ���������. ��������� ��������� ��������
% ������������
%
% ��������� ImportFromMatAndExcel, matrixABGinit
%
%***********************************************
tic
%clear;
%clc;
%% ������ � ��
prj.name = '����-2.1�';
prj.pn = 2175;
prj.ndus = 1; % Xbckj
prj.ntonsall = 3; % ����� ��������� ����� ������� ��������� 
%% ����������� ������������ ������
prj.ntons = 3; %{0,1,2,3} ����� ����������� ����� ������� ��������� ������� ������ 
prj.ntanks = 2; % {0,2} ����� ����� � ������ �����������
%����� �����
prj.spaceport = '�������';%{'��������','�������'}
prj.season = '����';%{'����','����'}
prj.windtype = '���������';%{'���������','��������'}
%�����
prj.mode = '����������';

%prj.path = 'C:\Users\admin\Documents\MATLAB\';%���� �� ������
prj.path = 'D:\YandexDisk\Documents\Ivan\!!Tasks&Projects\RocketStabilizationSystem\';%���� �� �����
prj.id = [prj.path,prj.name...
    ' �� = ' num2str(prj.pn)];
dm.T0 = 0.03; %������ �������������,� 
run([prj.path,'ImportFromMatAndExcel.m']);
run([prj.path,'matrixABGinit.m']);
as.data = (xlsread([prj.id,' ��.xlsx'],...
    '��������� �������� ������������','C1:L5'))';
as.t = as.data(:,1);
as.a0 = as.data(:,2);
as.a1 = as.data(:,3);
as.n0 = as.data(:,4);
as.n1 = as.data(:,5);
as.l = interp1(dm.t,as.l,as.t);
as = rmfield(as,'data');
clearvars options
toc;
% disp('�������������...')
% tic;
% [dm.Psi,dm.Phi] = CntrlTrMtrx(dm.A,dm.B,dm.T0);
% options = odeset('RelTol',1e-9,'AbsTol',1e-12);
% [dm.Gamma,~] = CntrlTrMtrx(dm.A,dm.G,dm.T0,options);
% toc;
disp('������ �������������...')
sim('model.slx');
result.t = ymeasure.time;
result.y = ymeasure.signals.values;
result.Vy = Vymeasure1.signals.values;
result.thet = thetmeasure.signals.values;
result.dthet = dthetmeasure.signals.values;
clearvars ymeasure Vymeasure1 thetmeasure dthetmeasure;
%% graphics
subplot(4,1,1);
plot(result.t,result.thet); grid on;
xlabel('t,c');
ylabel('\vartheta, ����');
subplot(4,1,2);
plot(result.t,result.y); grid on;
xlabel('t,c');
ylabel('y, �');
subplot(4,1,3);
plot(result.t,result.dthet); grid on;
xlabel('t,c');
ylabel('d\vartheta/dt, ����/c');
subplot(4,1,4);
plot(result.t,result.Vy); grid on;
xlabel('t,c');
ylabel('V_y, �/�');