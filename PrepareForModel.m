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
prj.ntons = 0; %{0,1,2,3} ����� ����������� ����� ������� ��������� ������� ������ 
prj.ntanks = 0; % {0,2} ����� ����� � ������ �����������
%����� �����
prj.spaceport = '��������';%{'��������','�������'}
prj.season = '����';%{'����','����'}
prj.windtype = '���������';%{'���������','��������'}
%�����
prj.mode = '����������';

%prj.path = 'C:\Users\admin\Documents\MATLAB\';%���� �� ������
prj.path = 'D:\YandexDisk\Documents\Ivan\!!Tasks&Projects\hard rocket\Hard_missile_model\';%���� �� �����
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
Q = [1/(3*pi/180)^2 0 0 0;
    0   1/500^2 0 0;
    0 0 1/(1*pi/180)^2 0;
    0 0 0 1/(10^2)];
for tmp = 1:length(dm.t)
    as.k(:,:,tmp) =lqr(dm.A(:,:,tmp),dm.B(:,:,tmp),eye(4),1);
end
as = rmfield(as,'data');
clearvars options
toc;
lqr.a0(1:length(dm.t)) = as.k(1,1,:);
lqr.a0a1(1:length(dm.t)) = as.k(1,3,:);
lqr.a0n1(1:length(dm.t)) = as.k(1,4,:);
lqr.a0n1n0(1:length(dm.t)) = as.k(1,2,:);
figure;
subplot(2,2,1); plot(dm.t,-lqr.a0,as.t,as.a0); grid on;
subplot(2,2,2); plot(dm.t,-lqr.a0n1n0,as.t,as.a0.*as.n1.*as.n0); grid on;
subplot(2,2,3); plot(dm.t,-lqr.a0a1,as.t,as.a0.*as.a1); grid on;
subplot(2,2,4); plot(dm.t,-lqr.a0n1,as.t,as.a0.*as.n1); grid on;
% disp('�������������...')
% tic;
% [dm.Psi,dm.Phi] = CntrlTrMtrx(dm.A,dm.B,dm.T0);
% options = odeset('RelTol',1e-9,'AbsTol',1e-12);
% [dm.Gamma,~] = CntrlTrMtrx(dm.A,dm.G,dm.T0,options);
% toc;
%% 
disp('������ �������������...')
sim('model.slx');
result.t = ymeasure.time;
result.y = ymeasure.signals.values;
result.Vy = Vymeasure1.signals.values;
result.thet = thetmeasure.signals.values;
result.dthet = dthetmeasure.signals.values;
clearvars ymeasure Vymeasure1 thetmeasure dthetmeasure;
%% graphics
figure;
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