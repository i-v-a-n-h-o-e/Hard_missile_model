%trj.thetpr = interp1(dm.t,trj.thetpr,result.t);
plot(result.t,trj.thetpr*pi/180+result.thet)
hold on;plot(result.t,trj.thetpr*pi/180);grid on;
