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
I = eye(n);
O = zeros(n);
dm.A = zeros(2*n,2*n,length(dm.t));
dm.B = zeros(2*n,1,length(dm.t));
dm.G = zeros(2*n,3,length(dm.t));
dm.C = zeros(2+prj.ndus,2*n,length(dm.t));
if (prj.ntons == 0) && (prj.ntanks == 0)
    for tmp = 1:length(dm.t)
        dm.A(:,:,tmp) = [0 0 1 0;
             0 0 0 1;
            -trj.Cthetatheta(tmp) 0 0 -trj.CthetaVy(tmp);
            -trj.CVytheta(tmp) 0 0 -trj.CVyVy(tmp)];
        dm.B(:,:,tmp) = [0;
                   0;
             -trj.Cthetadelta(tmp);
             -trj.CVydelta(tmp)];
         dm.C = [1 0 0 0;
             0 0 1 0;
             0 0 trj.l(tmp) 1];
         dm.G = [
             0 0 0;
             0 0 0;
             -trj.CthetaVy(tmp).*sind(trj.thetpr(tmp)).*cos(Az-pi-par.PS) 1./trj.Jz(tmp) 0;
             -trj.CVyVy(tmp).*sind(trj.thetpr(tmp)).*cos(Az-pi-par.PS)  0 1./trj.m(tmp)];
    end; 
end
for tmp = 1:length(dm.t)
    R = I;
    if (prj.ntanks > 0)
        R(3:(prj.ntanks+2),1) = (tanks.t(tmp,:))';
        R(3:(prj.ntanks+2),2) = -1.*(tanks.u(tmp,:))';
        R(1,3:(prj.ntanks+2)) = (tanks.T(tmp,:))';
        R(2,3:(prj.ntanks+2)) = -1.*(tanks.U(tmp,:))';
    end;    
    M = O;
    N = O;
    M(1,1) = -trj.Cthetatheta(tmp);
    M(2,1) = -trj.CVytheta(tmp);
    if (prj.ntanks > 0)
        M(3:(prj.ntanks+2),1) = -1*(tanks.r(tmp,:))';
        M(1,3:(prj.ntanks+2)) = -1*(tanks.P(tmp,:))';
        M(1,(prj.ntanks+3):n) = -1*(oscil.nu(tmp,:))';
        M(2,(prj.ntanks+3):n) = -1*(oscil.mu(tmp,:))';
    end;   
    N(1,1) = -trj.CthetaVy(tmp);
    N(2,1) = -trj.CVyVy(tmp);
    if prj.ntanks > 0
        for tmp2 = 3:(prj.ntanks+2)
            M(tmp2,tmp2) = -tanks.r(tmp,tmp2-2);
            N(tmp2,tmp2) = -tanks.e(tmp,tmp2-2);
        end;
        for tmp2 = (prj.ntanks+3):n
            M(tmp2,tmp2) = -oscil.om2(tmp,tmp2-(prj.ntanks+2));
            N(tmp2,tmp2) = sqrt(oscil.om2(tmp,tmp2-(prj.ntanks+2))) * koef;
        end;
         dm.G(:,:,tmp) = [zeros(n,3);
        -trj.CthetaVy(tmp).*sind(trj.thetpr(tmp)).*cos(Az-pi-par.PS) 1./trj.Jz(tmp) 0;
        -trj.CVyVy(tmp).*sind(trj.thetpr(tmp)).*cos(Az-pi-par.PS)  0 1./trj.m(tmp);
        zeros(n-2-prj.ntons,3);
        (oscil.k(tmp,:))' zeros(prj.ntons,2)];
    else 
         dm.G(:,:,tmp) = [zeros(2,3);
        -trj.CthetaVy(tmp).*sind(trj.thetpr(tmp)).*cos(Az-pi-par.PS) 1./trj.Jz(tmp) 0;
        -trj.CVyVy(tmp).*sind(trj.thetpr(tmp)).*cos(Az-pi-par.PS)  0 1./trj.m(tmp)];
    end;     
    dm.B(:,:,tmp) = zeros(2*n,1);
    if prj.ntons ~= 0
        H = [-trj.Cthetadelta(tmp);
            -trj.CVydelta(tmp);
            zeros(prj.ntanks,1);
            (oscil.k(tmp,:))'];
        C1 = [1 zeros(1,n-1-prj.ntons) (oscil.dfx(2,:,tmp));
        zeros(1+prj.ndus,n)];
        C2 = [zeros(1,n);
        ones(prj.ndus,1) zeros(prj.ndus,n-1-prj.ntons) (oscil.dfx(prj.ndus:-1:1,:,tmp));
        trj.l(tmp) 1 zeros(1,n-2-prj.ntons) (oscil.fx(2,:,tmp))];
    else     
        H = [-trj.Cthetadelta(tmp);
            -trj.CVydelta(tmp);
            zeros(prj.ntanks,1);];
        C1 = [1 zeros(1,n-1-prj.ntons);
                zeros(2,n)];
        C2 = [zeros(1,n);
            1 zeros(1,n-1-prj.ntons);
            trj.l(tmp) 1 zeros(1,n-2-prj.ntons)];
    end;     
    dm.A(:,:,tmp) = [O I;
        R\M R\N];
    dm.B(:,:,tmp) = [zeros(n,1);
                         R\H];
    dm.C(:,:,tmp) = [C1 C2];  
end
clearvars tmp Az tmp2 R M N O I H C1 C2 koef