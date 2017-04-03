tic;
for i = 1:4
    for j=1:4
        a(:) = dm.A(i,j,:);
        subplot(4,4,4*(i-1)+j);
        plot(dm.t,a);
    end;
end;
toc;