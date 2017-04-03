%%*********************************************
%
% Расчет нестационарных матриц коэффициентов 
% для уравнений в форме Коши по сетке времени dm.t
%
% НЕОБХОДИМЫ СЧИТАННЫЕ ДАННЫЕ В СТРУКТУРУ TRJ!!!
% Предварительный запуск
% koefs или ImportFromMatAndExcel 
%
%***********************************************
n = 2 + prj.ntanks + prj.ntons;
koef = -0.06 / pi;
dm.A = zeros(2*n,2*n,length(dm.t));
dm.B = zeros(2*n,1,length(dm.t));
dm.G = zeros(2*n,3,length(dm.t));
dm.C = zeros(2+prj.ndus,2*n,length(dm.t));
for tmp = 1:length(dm.t)
        dm.A(:,:,tmp) = [0                0 1       0;
                         0                0 0       1;
                    -trj.Cthetatheta(tmp) 0 0 -trj.CthetaVy(tmp);
                    -trj.CVytheta(tmp)    0 0 -trj.CVyVy(tmp)];
        dm.B(:,:,tmp) = [0;
                         0;
                     -trj.Cthetadelta(tmp);
                     -trj.CVydelta(tmp)];
         dm.C(:,:,tmp) = [1 0 0 0;
                 0 0 1 0;
                 0 0 trj.l(tmp) 1];
         dm.G(:,:,tmp) = [0 0 0;
             0 0 0;
             -trj.CthetaVy(tmp).*sind(trj.thetpr(tmp)).*cos(Az-pi-par.PS) 1./trj.Jz(tmp) 0;
             -trj.CVyVy(tmp).*sind(trj.thetpr(tmp)).*cos(Az-pi-par.PS)  0 1./trj.m(tmp)];
end; 
clearvars tmp Az koef