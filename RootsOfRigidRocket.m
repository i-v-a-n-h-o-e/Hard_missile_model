%% *********************************************
%
% ѕќстроение действительных частей полюсов системы дл€ €вной проверки
% устойчивости
%
%***********************************************
p = zeros(length(t),3);
for tmp=1:length(t)
    a=[1;
    trj.CVyVy(tmp);
    trj.Cthetatheta(tmp);
    trj.Cthetatheta(tmp)*trj.CVyVy(tmp)-trj.CVytheta(tmp)*trj.CthetaVy(tmp)];
    p(tmp,:) = roots(a);
end   
figure;
plot(t, real(p));
grid on;
xlabel('t,c');
ylabel('real(\lambda)');
legend('\lambda_1','\lambda_2','\lambda_3');
title(' арта вещественных частей полюсов в зависимости от времени');