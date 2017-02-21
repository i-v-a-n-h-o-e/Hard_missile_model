figure;
plot(trj.t,trj.P / 9810);
grid on;
hold on;
[~,~,Pa,~] = atmoscoesa(trj.H);
P = 172 + 0.2983 * (63.98 + 29.24) - Pa * pi * 0.32^2 / 9810 - Pa * pi * 1.47777^2 /4 /9810; 
plot(trj.t,P);