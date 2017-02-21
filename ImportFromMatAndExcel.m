%% *********************************************
%
% Расчет коэффициентов твердого тела на основе ИД(данных о траектории).
% Расчет структуры TRJ по сетке времени dm.t
% НЕОБХОДИМЫ ФАЙЛ *.xls  в папке проекта
%
%***********************************************

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
% 8: theta, град - программный тангаж

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
trj.winddata = xlsread([prj.id,' ИД.xlsx'],[prj.spaceport,' ',prj.season,' ',prj.windtype]);
% Характеристики ветра
% Если огибающая, то
% 1: Высота H, км
% 2: Wсист, м/с
% 3: sigma, m/s
% 4: Wогиб, m/s
% 
% Если градиент, то 
% 1: Высота H, км
% 2: sigma, m/s
% 4: sigma sloy = 1km, m/s
% 6: sigma sloy = 2km, m/s
% 8: sigma sloy = 3km, m/s
if (prj.ntons > 0) && (prj.ntons <= prj.ntonsall)
    oscil.par = (xlsread([prj.id,' ИД.xlsx'],'УКК параметры колебаний'))';
    oscil.par(1,:) = [];
    oscil.form = (xlsread([prj.id,' ИД.xlsx'],'УКК формы колебаний'))';
    oscil.form(1,:) = [];
end
if (prj.ntanks > 0)
    tanks.data1 = xlsread([prj.id,' ИД.xlsx'],'ЖНБ бак №1');
    tanks.data2 = xlsread([prj.id,' ИД.xlsx'],'ЖНБ бак №2');
end
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
        par.PS = 0;
    otherwise 
        clear;
        clc;
        error('Неправильная ракета');
end;

% switch prj.pn
%     case 2175    
     par.L = 43.14;% m 
     par.t12 = 196.94;% c 
     Az = pi;
%     otherwise 
%         clear;
%         clc;
%         error('Неправильная полезная нагрузка');
% end;
dm.t = [0:dm.T0:par.t12]';
    
%% Считываем данные из data

trj.H = 1000 * interp1(trj.data(:,1),trj.data(:,5),dm.t);
trj.V = interp1(trj.data(:,1),trj.data(:,2),dm.t);
trj.m = 1000 .* interp1(trj.data(:,1),trj.data(:,3),dm.t);
trj.P = interp1(trj.data(:,1),trj.data(:,4),dm.t) .* 1000 * 9.81;
trj.M = interp1(trj.data(:,1),trj.data(:,7),dm.t);
trj.q = interp1(trj.data(:,1),trj.data(:,6),dm.t) * 9.81;
trj.cya = interp1(trj.aerodata(:,1),trj.aerodata(:,2),trj.M) * 180 / pi;
trj.cyaalpha = interp1(trj.aerodata(:,1),trj.aerodata(:,3),trj.M) * 180 / pi;
trj.Jz = interp1(trj.massdata(:,1),trj.massdata(:,2),trj.m);
trj.Cm = interp1(trj.massdata(:,1),trj.massdata(:,3),trj.m) / 100;
trj.Cd = interp1(trj.aerodata(:,1),trj.aerodata(:,4),trj.M) / 100;
trj.thetpr = interp1(trj.data(:,1),trj.data(:,8),dm.t); %degrees
trj.l = par.xgp - par.L * trj.Cm;
as.l = trj.l;
%% Считываем ветер
switch prj.windtype
    case 'Градиент'
         ;
    case 'Огибающая'
        trj.wind = interp1(trj.winddata(:,1) .* 1000,trj.winddata(:,4)-trj.winddata(:,2),min(trj.H,20000));
    otherwise
        clear;
        clc;
        error('Неправильный ветер');
end;
%% Считываем параметры упругих колебаний корпуса
if (prj.ntons > 0) && (prj.ntons <= prj.ntonsall)
    for tmp = 1:prj.ntons
        oscil.om2(:,tmp) = interp1(oscil.par(:,1),oscil.par(:,1+tmp),dm.t);
        oscil.k(:,tmp) = interp1(oscil.par(:,1),oscil.par(:,4+tmp),dm.t);
        oscil.mu(:,tmp) = interp1(oscil.par(:,1),oscil.par(:,7+tmp),dm.t);
        oscil.nu(:,tmp) = interp1(oscil.par(:,1),oscil.par(:,10+tmp),dm.t);
        
        oscil.fx(1,tmp,:) = interp1(oscil.form(:,1),oscil.form(:,1+tmp),dm.t);
        oscil.fx(2,tmp,:) = interp1(oscil.form(:,1),oscil.form(:,4+tmp),dm.t);
        oscil.dfx(1,tmp,:) = interp1(oscil.form(:,1),oscil.form(:,7+tmp),dm.t);
        oscil.dfx(2,tmp,:) = interp1(oscil.form(:,1),oscil.form(:,10+tmp),dm.t);
    end
end  
if (prj.ntanks > 0)
    tanks.r(:,1) = interp1(tanks.data1(:,1),tanks.data1(:,2),dm.t);
    tanks.r(:,2) = interp1(tanks.data2(:,1),tanks.data2(:,2),dm.t);
    
    tanks.u(:,1) = interp1(tanks.data1(:,1),tanks.data1(:,3),dm.t);
    tanks.u(:,2) = interp1(tanks.data2(:,1),tanks.data2(:,3),dm.t);
    
    tanks.t(:,1) = interp1(tanks.data1(:,1),tanks.data1(:,4),dm.t);
    tanks.t(:,2) = interp1(tanks.data2(:,1),tanks.data2(:,4),dm.t);
    
    tanks.T(:,1) = interp1(tanks.data1(:,1),tanks.data1(:,5),dm.t);
    tanks.T(:,2) = interp1(tanks.data2(:,1),tanks.data2(:,5),dm.t);
    
    tanks.P(:,1) = interp1(tanks.data1(:,1),tanks.data1(:,6),dm.t);
    tanks.P(:,2) = interp1(tanks.data2(:,1),tanks.data2(:,6),dm.t);
    
    tanks.U(:,1) = interp1(tanks.data1(:,1),tanks.data1(:,7),dm.t);
    tanks.U(:,2) = interp1(tanks.data2(:,1),tanks.data2(:,7),dm.t);
    
    tanks.e(:,1) = interp1(tanks.data1(:,1),tanks.data1(:,8),dm.t);
    tanks.e(:,2) = interp1(tanks.data2(:,1),tanks.data2(:,8),dm.t);
end
%% Рассчитываем некоторые параметры
[~,~,trj.Pa,~] = atmoscoesa(trj.H);
trj.Pupr = par.Jupr * (par.mo + par.mf) / 4 - trj.Pa * par.Supr;
%% Формирование режима работы системы стабилизации
switch prj.mode
    case 'Повышенный'
        trj.cya = trj.cya * (1 - 0.12);
        trj.Cm = trj.Cm + 0.012;
        trj.Cd = trj.Cd - 0.03;
        trj.Jz = trj.Jz * (1 + 0.15);
        trj.P = trj.P * (1 + 0.015); 
        trj.Pupr = trj.Pupr * (1 + 0.015);
        trj.cyaalpha = trj.cyaalpha * (1 - 0.12);
    case 'Пониженный'
        trj.cya = trj.cya * (1 + 0.12);
        trj.Cm = trj.Cm - 0.012;
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
clearvars n tmp
trj = rmfield(trj,'data');
trj = rmfield(trj,'winddata');
trj = rmfield(trj,'massdata');
trj = rmfield(trj,'aerodata');
if (prj.ntons > 0)
     oscil = rmfield(oscil,'par');
     oscil = rmfield(oscil,'form');
end;
if (prj.ntanks > 0)
    tanks = rmfield(tanks,'data1');
    tanks = rmfield(tanks,'data2');
end;   