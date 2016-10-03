N=1:100;
fig = figure;

sigma_w1=0.01;
Y1=2*qfunc(pi/4./sqrt(N*sigma_w1^2));
plot(N,Y1);
xlabel('N (#verstuurde symbolen)');
ylabel('SER');
title('Symboolfoutprobabiliteit voor \sigma_w = 0.01');
print(fig, '-dpng', 'SER1.png');

sigma_w2=0.05;
Y2=2*qfunc(pi/4./sqrt(N*sigma_w2^2));
plot(N,Y2);
xlabel('N (#verstuurde symbolen)');
ylabel('SER');
title('Symboolfoutprobabiliteit voor \sigma_w = 0.05');
print(fig, '-dpng', 'SER2.png');

sigma_w3=0.1;
Y3=2*qfunc(pi/4./sqrt(N*sigma_w3^2));
plot(N,Y3);
xlabel('N (#verstuurde symbolen)');
ylabel('SER');
title('Symboolfoutprobabiliteit voor \sigma_w = 0.1');
print(fig, '-dpng', 'SER3.png');

sigma_w4=0.5;
Y4=2*qfunc(pi/4./sqrt(N*sigma_w4^2));
plot(N,Y4);
xlabel('N (#verstuurde symbolen)');
ylabel('SER');
title('Symboolfoutprobabiliteit voor \sigma_w = 0.5');
print(fig, '-dpng', 'SER4.png');

plot(N,Y1,N,Y2,N,Y3,N,Y4);
xlabel('N (#verstuurde symbolen)');
ylabel('SER');
title('Symboolfoutprobabiliteit');
legend('\sigma_w=0.01','\sigma_w=0.05','\sigma_w=0.1','\sigma_w=0.5','Location','Northeast');
print(fig, '-dpng', 'SER.png');
