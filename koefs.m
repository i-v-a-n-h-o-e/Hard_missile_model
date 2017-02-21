%%*********************************************
%
% Расчет коэффициентов твердого тела и построение их графиков
% на основе ИД.
% НЕОБХОДИМ ФАЙЛ *.xls в папке проекта
%
%***********************************************
tic
clear;
clc;

prj.name = 'Союз-2.1в';
prj.pn = 2175;

prj.mode = 'Номинальный';

%'C:\Users\admin\Documents\MATLAB\'
prj.id = ['D:\YandexDisk\Documents\Ivan\!!Tasks&Projects\RocketStabilizationSystem\',prj.name...
    ' ПН = ' num2str(prj.pn)];
dm.T0 = 0.03; %период дискретизации,с 
%% Считываем данные из EXCEL
trj.data = xlsread([prj.id,' ИД.xlsx'],'Траектория');
% Параметры траектории
% 1: t,c - время
% 2: V, m/s - скорость
% 3: m, т - масса
% 4: P, тс - суммарная тяга
% 5: H, км - высота
% 6: q, кгс/м2 - скоростной напор
% 7: М, б/р - число Маха
% 2: theta, град - программный тангаж

trj.massdata = xlsread([prj.id,' ИД.xlsx'],'Массовые центровочные хар-ки');
% Массовые центровочные характеристики
% 1: m, кг - масса
% 2: Jzz, кг м2 - осевой момент инерции 
% 3: Cm, % - относительное полоэение центра масс
trj.aerodata = xlsread([prj.id,' ИД.xlsx'],'Аэродинамические хар-ки');
% Массовые центровочные характеристики
% 1: М, б/р - число Маха
% 2: cya, 1/град
% 3: сyaalpha, 1/град
% 4: Cd, % - относительное положение центра давления

%% Выбираем ракету и полезную нагрузку
switch prj.name
    case 'Союз-2.1в'    
        par.Jupr = 298.3 * 9.81; % Н * с / кг  
        par.Sm = 6.82; % м^2
        par.mo = 62.306; % kg/s
        par.mf = 28.398; % kg/s
        par.Supr = pi * (0.32 ^ 2) / 4; % m^2
        par.Xupr = 0.175; % m
        par.nuo = 2;% Число управляющих органов в канале тангажа(рыскания)
        par.xgp = 32.4215;
    otherwise 
        clear;
        clc;
        error('Неправильная ракета');
end;

% switch prj.pn
%     case 2175    
     par.L = 43.14;% m высота ракеты
     par.t12 = 196.94;% c, время до разделения первой и второй ступени
%     otherwise 
%         clear;
%         clc;
%         error('Неправильная полезная нагрузка');
% end;
t = [0:dm.T0:par.t12]'; %сетка времени
%% Считываем данные из data

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
%% Рассчитываем некоторые параметры
[~,~,trj.Pa,~] = atmoscoesa(trj.H);
trj.Pupr = par.Jupr * (par.mo + par.mf) / 4 - trj.Pa * par.Supr;

%% Формирование режима работы системы стабилизации
switch prj.mode
    case 'Повышенный'
        trj.cya = trj.cya * (1 - 0.12);
        trj.Cm = trj.Cm + 0.12;
        trj.Cd = trj.Cd - 0.03;
        trj.Jz = trj.Jz * (1 + 0.15);
        trj.P = trj.P * (1 + 0.015);
        trj.Pupr = trj.Pupr * (1 + 0.015);
        trj.cyaalpha = trj.cyaalpha * (1 - 0.12);
    case 'Пониженный'
        trj.cya = trj.cya * (1 + 0.12);
        trj.Cm = trj.Cm - 0.12;
        trj.Cd = trj.Cd + 0.03;
        trj.Jz = trj.Jz * (1 - 0.15);
        trj.P = trj.P * (1 - 0.015);
        trj.Pupr = trj.Pupr * (1 - 0.015);
        trj.cyaalpha = trj.cyaalpha * (1 + 0.12);
end
%% Рассчитываем коэффициенты

trj.Cthetatheta = trj.cya .* trj.q .* par.Sm .* par.L ...
    .* (trj.Cm - trj.Cd) ./ trj.Jz;
trj.CthetaVy = -trj.Cthetatheta ./ trj.V;
trj.Cthetadelta = par.nuo * trj.Pupr .* (trj.Cm .* par.L - par.Xupr) ./ trj.Jz;
trj.CVytheta = - (trj.P + trj.cyaalpha .* trj.q .* par.Sm) ./ trj.m;
trj.CVyVy = trj.cya .* trj.q .* par.Sm ./ (trj.m .* trj.V);
trj.CVydelta = - par.nuo * trj.Pupr ./ trj.m;
%% Строим графики

subplot(2,3,1); hold on;
plot(t,trj.Cthetatheta);
xlabel('t,c');
ylabel('C_{{\vartheta}\vartheta}, 1/c^2');
grid on;

subplot(2,3,2); hold on;
plot(t,trj.CthetaVy);
xlabel('t,c');
ylabel('C_{\varthetaV_y}, (м*с)^{-1}');
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
ylabel('C_{V_yV_y}, (м*с)^{-1}');
grid on;

subplot(2,3,6); hold on;
plot(t,trj.CVydelta);
xlabel('t,c');
ylabel('C_{V_y\delta}, 1/c^2');
grid on;
%title(gca,'Графики изменения коэффициентов уравнений возмущенного движения')
toc;