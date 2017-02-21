%%*********************************************
%
% ������ ������������� �������� ���� � ���������� �� ��������
% �� ������ ��.
% ��������� ���� *.xls � ����� �������
%
%***********************************************
tic
clear;
clc;

prj.name = '����-2.1�';
prj.pn = 2175;

prj.mode = '�����������';

%'C:\Users\admin\Documents\MATLAB\'
prj.id = ['D:\YandexDisk\Documents\Ivan\!!Tasks&Projects\RocketStabilizationSystem\',prj.name...
    ' �� = ' num2str(prj.pn)];
dm.T0 = 0.03; %������ �������������,� 
%% ��������� ������ �� EXCEL
trj.data = xlsread([prj.id,' ��.xlsx'],'����������');
% ��������� ����������
% 1: t,c - �����
% 2: V, m/s - ��������
% 3: m, � - �����
% 4: P, �� - ��������� ����
% 5: H, �� - ������
% 6: q, ���/�2 - ���������� �����
% 7: �, �/� - ����� ����
% 2: theta, ���� - ����������� ������

trj.massdata = xlsread([prj.id,' ��.xlsx'],'�������� ������������ ���-��');
% �������� ������������ ��������������
% 1: m, �� - �����
% 2: Jzz, �� �2 - ������ ������ ������� 
% 3: Cm, % - ������������� ��������� ������ ����
trj.aerodata = xlsread([prj.id,' ��.xlsx'],'���������������� ���-��');
% �������� ������������ ��������������
% 1: �, �/� - ����� ����
% 2: cya, 1/����
% 3: �yaalpha, 1/����
% 4: Cd, % - ������������� ��������� ������ ��������

%% �������� ������ � �������� ��������
switch prj.name
    case '����-2.1�'    
        par.Jupr = 298.3 * 9.81; % � * � / ��  
        par.Sm = 6.82; % �^2
        par.mo = 62.306; % kg/s
        par.mf = 28.398; % kg/s
        par.Supr = pi * (0.32 ^ 2) / 4; % m^2
        par.Xupr = 0.175; % m
        par.nuo = 2;% ����� ����������� ������� � ������ �������(��������)
        par.xgp = 32.4215;
    otherwise 
        clear;
        clc;
        error('������������ ������');
end;

% switch prj.pn
%     case 2175    
     par.L = 43.14;% m ������ ������
     par.t12 = 196.94;% c, ����� �� ���������� ������ � ������ �������
%     otherwise 
%         clear;
%         clc;
%         error('������������ �������� ��������');
% end;
t = [0:dm.T0:par.t12]'; %����� �������
%% ��������� ������ �� data

trj.H = 1000 * interp1(trj.data(:,1),trj.data(:,5),t);
trj.t = trj.data(:,1);
trj.V = interp1(trj.data(:,1),trj.data(:,2),t);
trj.m = 1000 .* interp1(trj.data(:,1),trj.data(:,3),t);
trj.P = interp1(trj.data(:,1),trj.data(:,4),t) .* 1000 * 9.81;
trj.M = interp1(trj.data(:,1),trj.data(:,7),t);
trj.q = interp1(trj.data(:,1),trj.data(:,6),t) * 9.81;
trj.cya = interp1(trj.aerodata(:,1),trj.aerodata(:,2),trj.M) * 180 / pi;
trj.cyaalpha = interp1(trj.aerodata(:,1),trj.aerodata(:,3),trj.M) * 180 / pi;
trj.Jz = interp1(trj.massdata(:,1),trj.massdata(:,2),trj.m);
trj.Cm = interp1(trj.massdata(:,1),trj.massdata(:,3),trj.m) / 100;
trj.Cd = interp1(trj.aerodata(:,1),trj.aerodata(:,4),trj.M) / 100;
trj.l = par.xgp - par.L * trj.Cm; 
%% ������������ ��������� ���������
[~,~,trj.Pa,~] = atmoscoesa(trj.H);
trj.Pupr = par.Jupr * (par.mo + par.mf) / 4 - trj.Pa * par.Supr;

%% ������������ ������ ������ ������� ������������
switch prj.mode
    case '����������'
        trj.cya = trj.cya * (1 - 0.12);
        trj.Cm = trj.Cm + 0.12;
        trj.Cd = trj.Cd - 0.03;
        trj.Jz = trj.Jz * (1 + 0.15);
        trj.P = trj.P * (1 + 0.015);
        trj.Pupr = trj.Pupr * (1 + 0.015);
        trj.cyaalpha = trj.cyaalpha * (1 - 0.12);
    case '����������'
        trj.cya = trj.cya * (1 + 0.12);
        trj.Cm = trj.Cm - 0.12;
        trj.Cd = trj.Cd + 0.03;
        trj.Jz = trj.Jz * (1 - 0.15);
        trj.P = trj.P * (1 - 0.015);
        trj.Pupr = trj.Pupr * (1 - 0.015);
        trj.cyaalpha = trj.cyaalpha * (1 + 0.12);
end
%% ������������ ������������

trj.Cthetatheta = trj.cya .* trj.q .* par.Sm .* par.L ...
    .* (trj.Cm - trj.Cd) ./ trj.Jz;
trj.CthetaVy = -trj.Cthetatheta ./ trj.V;
trj.Cthetadelta = par.nuo * trj.Pupr .* (trj.Cm .* par.L - par.Xupr) ./ trj.Jz;
trj.CVytheta = - (trj.P + trj.cyaalpha .* trj.q .* par.Sm) ./ trj.m;
trj.CVyVy = trj.cya .* trj.q .* par.Sm ./ (trj.m .* trj.V);
trj.CVydelta = - par.nuo * trj.Pupr ./ trj.m;
%% ������ �������

subplot(2,3,1); hold on;
plot(t,trj.Cthetatheta);
xlabel('t,c');
ylabel('C_{{\vartheta}\vartheta}, 1/c^2');
grid on;

subplot(2,3,2); hold on;
plot(t,trj.CthetaVy);
xlabel('t,c');
ylabel('C_{\varthetaV_y}, (�*�)^{-1}');
grid on;

subplot(2,3,3); hold on;
plot(t,trj.Cthetadelta);
xlabel('t,c');
ylabel('C_{\vartheta\delta},1/c^2');
grid on;

subplot(2,3,4); hold on;
plot(t,trj.CVytheta);
xlabel('t,c');
ylabel('C_{V_y\vartheta}, 1/c^2');
grid on;

subplot(2,3,5); hold on;
plot(t,trj.CVyVy);
xlabel('t,c');
ylabel('C_{V_yV_y}, (�*�)^{-1}');
grid on;

subplot(2,3,6); hold on;
plot(t,trj.CVydelta);
xlabel('t,c');
ylabel('C_{V_y\delta}, 1/c^2');
grid on;
%title(gca,'������� ��������� ������������� ��������� ������������ ��������')
toc;