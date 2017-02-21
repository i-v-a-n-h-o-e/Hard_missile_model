%%*********************************************
%
% Построение АФЧХ для конкретного момента времени с 
% настройками автомата стаблизации в этот момент времени 
%
% Запускает ImportFromMatAndExcel, matrixABGinit, contour
% Необходим boderus
%
%***********************************************
tic
%clear;
%clc;

%% Ракета и ПН
prj.name = 'Союз-2.1в';
prj.pn = 2175;
prj.ndus = 1;
prj.ntonsall = 3; % Число доступных тонов упргуих колебаний 
%% Подробность рассмотрения модели
prj.ntons = 3; %{0,1,2,3} Число учитываемых тонов упругих колебаний корпуса ракеты 
prj.ntanks = 2; % {0,2} Число баков с жидким наполнением
%Выбор ветра
prj.spaceport = 'Байконур';%{'Байконур','Плесецк'}
prj.season = 'Лето';%{'Лето','Зима'}
prj.windtype = 'Огибающая';%{'Огибающая','Градиент'}

%prj.mode = 'Номинальный';
%prj.mode = 'Повышенный';
prj.mode = 'Пониженный';

%prj.path = 'C:\Users\admin\Documents\MATLAB\';%Путь на работе
prj.path = 'D:\YandexDisk\Documents\Ivan\!!Tasks&Projects\RocketStabilizationSystem\';%Путь на ноуте  
prj.id = [prj.path,prj.name...
    ' ПН = ' num2str(prj.pn)];

dm.T0 = 0.03; %период дискретизации,с 

t = 80; %рассматриваемый момент времени
as.a0 = 19; % Общий пропорциоанльный
as.a1 = 0.8; % Дифференциальный по тангажу
as.n0 = 0; % Интегральный по ЦМ
as.n1 = 0.02; %Пропорциальныq по ЦМ

run([prj.path,'ImportFromMatAndExcel.m']);
run([prj.path,'matrixABGinit.m']);
run([prj.path,'contour.m']);

boderus_A_v2(W,3);
title({[prj.name, ' Режим: ',prj.mode, ' t = ',num2str(t),' c'],'Канал тангажа'});
figure;
step(feedback(W,1,1));
grid on;
toc
