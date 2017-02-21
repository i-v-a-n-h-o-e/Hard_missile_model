%%*************************************************************************
%
% Считывание моделей Датчика Угловых Скоростей (ДУС) и рулевого привода
% из файла .mat
% Сборка контура системы в момент времени t. Система состоит из Объекта управления
% (Матрицы A,B,C,G  dm. ), Автомата стабилизации(as), рулевого привода 
% Расчет нестационарных матриц коэффициентов 
% для уравнений в форме Коши по сетке времени dm.t
%
% НЕОБХОДИМЫ момент времени t
% матрицы А, B, C, G и файл *.mat
% Предварительный запуск matrixABGinit.m
%
%**************************************************************************

load([prj.path, 'DUS&rp.mat']); %Загрузка моделей ДУС и РП
sl = ceil(t/dm.T0) + 1; %Номер слоя матриц A,B,C,G
dWobj = ss(dm.Phi(:,:,sl),dm.Psi(:,:,sl),dm.C(:,:,sl),0,dm.T0);
l = as.l(sl);
dWdus = c2d(Wdus,dm.T0);
dWrp = c2d(Wrp,dm.T0);
Cnt = append(dWrp,dWobj,dWdus*eye(prj.ndus));
Q = [2 1;
     (3:2+prj.ndus)' (3:2+prj.ndus)'];
inputs = 1;
if prj.ndus == 1
    a = 5;
else
   a = (6:7)';
end;   
outputs = [2; a; 3+prj.ndus];
dW1 = connect(Cnt,Q,inputs, outputs);
clearvars Cnt Q inputs outputs
%set(W1,'OutputDelay',[0.0214; 0.0114*ones(prj.ndus,1); 0.0214]);
%W1 = c2d(W1,dm.T0);
di = filt([0 dm.T0],[1 -1],dm.T0);
if prj.ndus == 2
    alpha = oscil.dfx(2,1,t)/(oscil.dfx(2,1,t)-oscil.dfx(1,1,t));
    dcomplex = c2d(ss([1-alpha alpha]));
    Cnt = append(dW1,dcomplex,as.a0,as.a1,as.n0*di,as.n1,l);
    Q = [3 3 0 0;
        2 2 0 0 ;
        5 5 0 0;
        7 4 -10 8;
        6 4 -10 0;
        8 5 0 0;
        4 1 7 9];
        outputs = 6;
else    
    Cnt = append(dW1,as.a0,as.a1,as.n0*di,as.n1,l);
    Q = [3 2 0 0;
        6 2 0 0;
        5 3 -8 6;
        4 3 -8 0;
        2 1 5 7];
        outputs = 4;
end;
inputs = 1;
dW = connect(Cnt,Q,inputs, outputs);
clearvars Cnt Q inputs outputs sl